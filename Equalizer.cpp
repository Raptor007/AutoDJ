#include "Equalizer.h"

#include <cmath>
extern "C" {
#include <libavutil/opt.h>
}


// --------------------------------------------------------------------------------------


EqualizerParam::EqualizerParam( void )
{
}


EqualizerParam::~EqualizerParam()
{
}


float EqualizerParam::LookupScale( float freq ) const
{
	std::map<float,float>::const_iterator found = FreqScale.lower_bound( freq );
	
	// If it's higher than our highest specified, use the highest.
	if( found == FreqScale.end() )
		return FreqScale.size() ? FreqScale.rbegin()->second : 1.f;
	
	// If it's less than or equal to our lowest specified, use the lowest.
	if( found == FreqScale.begin() )
		return found->second;
	
	std::map<float,float>::const_iterator prev = found;
	prev --;
	
	// Interpolate linearly between nearest frequency volume scales.
	float b_part = (freq - prev->first) / (found->first - prev->first);
	return prev->second * (1. - b_part) + found->second * b_part;
}


float EqualizerParam::GetScale( float freq ) const
{
	#define TRANSITION_AT (20094.)
	#define TRANSITION_TO (20627.)
	
	// Use EQ specifications for audible frequencies.
	if( freq <= TRANSITION_AT )
		return LookupScale(freq);
	
	// Above the transition band, use 0Hz scale.
	if( freq >= TRANSITION_TO )
		return LookupScale(0.f);
	
	// Within the transition band, interpolate from top EQ scale to 0Hz scale.
	float end_part = (freq - TRANSITION_AT) / (TRANSITION_TO - TRANSITION_AT);
	return LookupScale(TRANSITION_AT) * (1.f - end_part) + LookupScale(0.f) * end_part;
}


float EqualizerParam::GetScale( size_t index, size_t frames, unsigned int rate ) const
{
	return GetScale( ((float) index) * rate / frames );
}


float EqualizerParam::MaxScale( void ) const
{
	if( ! FreqScale.size() )
		return 1.f;
	
	std::map<float,float>::const_iterator iter = FreqScale.begin();
	float max = iter->second;
	for( iter ++; iter != FreqScale.end(); iter ++ )
	{
		if( iter->second > max )
			max = iter->second;
	}
	return max;
}


float EqualizerParam::MinScale( void ) const
{
	if( ! FreqScale.size() )
		return 1.f;
	
	std::map<float,float>::const_iterator iter = FreqScale.begin();
	float min = iter->second;
	for( iter ++; iter != FreqScale.end(); iter ++ )
	{
		if( iter->second < min )
			min = iter->second;
	}
	return min;
}


float EqualizerParam::AvgScale( void ) const
{
	if( ! FreqScale.size() )
		return 1.f;
	
	float total = 0.f;
	for( std::map<float,float>::const_iterator iter = FreqScale.begin(); iter != FreqScale.end(); iter ++ )
		total += iter->second;
	return total / FreqScale.size();
}


// --------------------------------------------------------------------------------------


EqualizerFFT::EqualizerFFT( unsigned int channels, unsigned int rate, size_t frames )
{
	Channels = channels;
	Rate = rate;
	Frames = frames;
	Context1 = av_fft_init( log2(Frames), false );
	Context2 = av_fft_init( log2(Frames), true );
	Complex = (FFTComplex*) av_mallocz( Frames * sizeof(FFTComplex) );
}


EqualizerFFT::~EqualizerFFT()
{
	av_fft_end( Context1 );
	av_fft_end( Context2 );
	av_free( Complex );
}


void EqualizerFFT::Process( float *buffer, EqualizerParam *param )
{
	for( size_t ch = 0; ch < Channels; ch ++ )
	{
		memset( Complex, 0, Frames * sizeof(FFTComplex) );
		for( size_t i = 0; i < Frames; i ++ )
			Complex[ i ].re = buffer[ (i * Channels) + ch ];
		av_fft_permute( Context1, Complex );
		av_fft_calc( Context1, Complex );
		
		// Apply equalizer to all frequency bins, and their negative bins.
		for( size_t i = 1; i < Frames / 2; i ++ )
		{
			float scale = param->GetScale( i, Frames, Rate );
			Complex[ i ].re *= scale;
			Complex[ i ].im *= scale;
			Complex[ Frames - i ].re *= scale;
			Complex[ Frames - i ].im *= scale;
		}
		
		float scale0 = param->GetScale(0.f);
		Complex[ 0 ].re *= scale0;
		Complex[ 0 ].im *= scale0;
		float scaleH = param->GetScale( Frames / 2, Frames, Rate );
		Complex[ Frames / 2 ].re *= scaleH;
		Complex[ Frames / 2 ].im *= scaleH;
		
		av_fft_permute( Context2, Complex );
		av_fft_calc( Context2, Complex );
		
		#define EQ_ANTIPOP (Rate / 500)  // 2ms each end = 88 samples each end at 44.1KHz
		float new_scale = 1.f / Frames;
		
		for( size_t i = EQ_ANTIPOP; i < Frames - EQ_ANTIPOP; i ++ )
			buffer[ (i * Channels) + ch ] = Complex[ i ].re * new_scale;
		
		size_t end_antipop = std::min<size_t>( Frames / 2, EQ_ANTIPOP );
		for( size_t i = 0; i < end_antipop; i ++ )
		{
			float new_part = i / (float) EQ_ANTIPOP;
			buffer[ (i * Channels) + ch ] = new_part * Complex[ i ].re * new_scale + (1.f - new_part) * scale0 * buffer[ (i * Channels) + ch ];
			size_t j = Frames - 1 - i;
			buffer[ (j * Channels) + ch ] = new_part * Complex[ j ].re * new_scale + (1.f - new_part) * scale0 * buffer[ (j * Channels) + ch ];
		}
	}
}


// --------------------------------------------------------------------------------------


Equalizer::Equalizer( void )
{
}


Equalizer::~Equalizer()
{
}


void Equalizer::Process( float *buffer, unsigned int channels, unsigned int rate, size_t frames, EqualizerParam *param )
{
#ifdef EQ_FRAMES
	if( frames > EQ_FRAMES )
	{
		for( size_t chunk = 0; (chunk + 1) * EQ_FRAMES <= frames; chunk ++ )
			Process( (float*) (((char*) buffer) + sizeof(*buffer) * channels * chunk * EQ_FRAMES), channels, rate, EQ_FRAMES, param );
		return;
	}
#endif
	
	if( param )
	{
		// NOTE: We can do this because our audio output never changes rate/channels.
		if( EqualizerFFTs.find(frames) == EqualizerFFTs.end() )
			EqualizerFFTs[ frames ] = new EqualizerFFT( channels, rate, frames );
		
		EqualizerFFTs[ frames ]->Process( buffer, param );
	}
}
