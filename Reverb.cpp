#include "Reverb.h"

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <stdint.h>

#ifndef likely
#ifdef __GNUC__
#define likely(x) (__builtin_expect( !!(x), 1 ))
#else
#define likely(x) (x)
#endif
#endif

#ifndef unlikely
#ifdef __GNUC__
#define unlikely(x) (__builtin_expect( !!(x), 0 ))
#else
#define unlikely(x) (x)
#endif
#endif


// --------------------------------------------------------------------------------------


ReverbBounce::ReverbBounce( void )
{
	FramesBack = 0;
	AmpScale = 0;
}


ReverbBounce::ReverbBounce( size_t f, float a )
{
	FramesBack = f;
	AmpScale = a;
}


ReverbBounce::~ReverbBounce()
{
}


// --------------------------------------------------------------------------------------


ReverbParam::ReverbParam( void )
{
	SpeakerSide = 1.0;
	SpeakerFront = 3.2;
	SideWall = 2.5;
	FrontWall = 4.2;
	BackWall = 0.8;
	Ceiling = 1.4;
	Floor = 0.9;
	HeadWidth = 0.15;
	BehindScale = 0.5;
	BounceEnergy = 0.3;
}


ReverbParam::~ReverbParam()
{
}


float ReverbParam::SpeakerDist( void ) const
{
	return sqrt( SpeakerSide * SpeakerSide + SpeakerFront * SpeakerFront );
}

	
float ReverbParam::BouncedDist( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite ) const
{
	float x = SpeakerSide;
	if( x_bounces % 2 )
		x += 2 * (opposite ? SideWall : SideWall - SpeakerSide);
	x += (x_bounces / 2) * SideWall * 4;
	x += HeadWidth * (opposite ? 0.5f : -0.5f);
	
	float y = SpeakerFront;
	if( y_bounces % 2 )
		y += 2 * BackWall;
	y += (y_bounces / 2) * (FrontWall + BackWall) * 2;
	
	float z = 0.f;
	if( z_bounces % 2 )
		z += 2 * (up ? Ceiling : Floor);
	z += (z_bounces / 2) * (Ceiling + Floor) * 2;
	
	return sqrt( x * x + y * y + z * z );
}


float ReverbParam::AmpScale( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite ) const
{
	if( !(x_bounces || y_bounces || z_bounces) )
		return 0.f;
	if( opposite && ! x_bounces )
		return 0.f;
	float bounced = BouncedDist( x_bounces, y_bounces, z_bounces, up, opposite );
	float speaker = SpeakerDist();
	float scale = pow( 0.5, bounced / speaker ) * pow( BounceEnergy, x_bounces + y_bounces + z_bounces );
	if( y_bounces % 2 ) // Arriving from behind.
		scale *= pow( BehindScale, y_bounces / (float)(x_bounces + y_bounces) );
	return scale;
}


void ReverbParam::AddBounce( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite, float speaker, int rate )
{
	float amp_scale = AmpScale( x_bounces, y_bounces, z_bounces, up, opposite );
	if( amp_scale )
	{
		size_t frames_back = (BouncedDist( x_bounces, y_bounces, z_bounces, up, opposite ) - speaker) * rate / 344.;
		if( opposite )
			OppoSide.push_back( ReverbBounce( frames_back, amp_scale ) );
		else
			SameSide.push_back( ReverbBounce( frames_back, amp_scale ) );
	}
}


void ReverbParam::Setup( unsigned int rate )
{
	SameSide.clear();
	OppoSide.clear();
	
	float speaker = SpeakerDist();
	
	for( int xb = 0; xb <= 4; xb ++ )
	{
		for( int yb = 0; yb <= 3; yb ++ )
		{
			for( int zb = 0; zb <= 4; zb ++ )
			{
				AddBounce( xb, yb, zb, false, false, speaker, rate );
				AddBounce( xb, yb, zb, true,  false, speaker, rate );
				AddBounce( xb, yb, zb, false, true,  speaker, rate );
				AddBounce( xb, yb, zb, true,  true,  speaker, rate );
			}
		}
	}
	
	TotalScale = 1.;
	size_t bounces_same = SameSide.size();
	size_t bounces_oppo = OppoSide.size();
	for( size_t i = 0; i < bounces_same; i ++ )
		TotalScale += SameSide[ i ].AmpScale;
	for( size_t i = 0; i < bounces_oppo; i ++ )
		TotalScale += OppoSide[ i ].AmpScale;
}


// --------------------------------------------------------------------------------------


Reverb::Reverb( void )
{
	History = NULL;
	Rate = 44100;
}


Reverb::~Reverb()
{
}

	
#define REVERB_HISTORY_FRAMES 262144
	
void Reverb::Process( float *buffer, unsigned int channels, unsigned int rate, size_t frames, ReverbParam *param )
{
	bool changed_rate = (Rate != rate);
	Rate = rate;
	
	size_t buff_bytes = sizeof(*History) * channels * frames;
	size_t hist_bytes = sizeof(*History) * channels * REVERB_HISTORY_FRAMES;
	uint8_t *raw_history = (uint8_t*) History;
	
	if( ! History )
	{
		History = (float*) malloc( hist_bytes );
		memset( History, 0, hist_bytes );
		raw_history = (uint8_t*) History;
	}
	else
		memmove( raw_history, raw_history + buff_bytes, hist_bytes - buff_bytes );
	memcpy( raw_history + hist_bytes - buff_bytes, buffer, buff_bytes );
	
	if( ! param )
		return;
	
	// If the reverb hasn't been calculated yet or the rate has changed, calculate now.
	if( changed_rate || (param->SameSide.size() + param->OppoSide.size() == 0) )
		param->Setup( rate );
	
	size_t bounces_same = param->SameSide.size();
	size_t bounces_oppo = param->OppoSide.size();
	
	for( unsigned int ch = 0; ch < channels; ch ++ )
	{
		for( size_t frame = 0; frame < frames; frame ++ )
		{
			size_t index = channels * frame + ch;
			float val = buffer[ index ];
			size_t frames_back = frames - 1 - frame;
			for( size_t i = 0; i < bounces_same; i ++ )
			{
				int from_frame = REVERB_HISTORY_FRAMES - 1 - frames_back - param->SameSide[ i ].FramesBack;
				if(likely( from_frame >= 0 ))
					val += History[ channels * from_frame + ch ] * param->SameSide[ i ].AmpScale;
			}
			for( size_t i = 0; i < bounces_oppo; i ++ )
			{
				int from_frame = REVERB_HISTORY_FRAMES - 1 - frames_back - param->OppoSide[ i ].FramesBack;
				if(likely( from_frame >= 0 ))
					val += History[ channels * from_frame + (ch+1)%channels ] * param->OppoSide[ i ].AmpScale;
			}
			buffer[ index ] = val;
		}
	}
}
