#define __STDC_CONSTANT_MACROS
#define __STDC_FORMAT_MACROS
#include <stdint.h>
#include <dirent.h>
#include <sys/stat.h>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <deque>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <SDL/SDL.h>
#include <SDL/SDL_audio.h>
#include "AudioFile.h"
#include "PlaybackBuffer.h"
#include "FontDraw.h"
#include "FontBin.h"
#ifdef WAVEOUT
#include <windows.h>
#include <mmreg.h>
#endif
extern "C" {
#include <libavcodec/avfft.h>
}

class Song;
class SongLoader;
class EqualizerParam;
class EqualizerFFT;
class Equalizer;
class ReverbBounce;
class ReverbParam;
class Reverb;
class UserData;

#define VISUALIZERS            5
#define VISUALIZER_COLORS     10
#define VISUALIZER_BACKGROUNDS 4

int SongLoaderThread( void *data_ptr );
float LinearCrossfade( float a, float b, float b_percent );
float EqualPowerCrossfade( float a, float b, float b_percent );
Sint16 FloatTo16( float value );
bool CalculateCrossfade( const Song *current_song, const Song *next_song, double *crossfade_for_beats, double *crossfade_at_beat );
void AudioCallback( void *userdata, Uint8* stream, int len );
void BufferedAudioCallback( void *userdata, Uint8* stream, int len );
void UnbufferedAudioCallback( void *userdata, Uint8* stream, int len );
void UpdateVisualizerBuffer( UserData *ud, Uint8 *stream, int len );
std::deque<std::string> DirSongs( std::string path );

#ifdef __GNUC__
#define ALLOW_VLA
#endif

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

#ifdef WIN32
const char *PATH_SEPARATOR = "\\";
#else
const char *PATH_SEPARATOR = "/";
#endif


// --------------------------------------------------------------------------------------


bool AlwaysLoadFloat = false;


class Song
{
public:
	AudioFile Audio;
	size_t TotalFrames;
	double BPM;
	double CurrentFrame, FirstBeatFrame, FirstOutroFrame;
	int IntroBeats, OutroBeats;
	std::map<size_t,double> Beats;
	float VolumeMax, VolumeAverage;
	volatile bool *RunningPtr;
	
	Song( const char *filename, volatile bool *running_ptr )
	{
		TotalFrames = 0;
		BPM = 138.;
		CurrentFrame = 0.;
		FirstBeatFrame = 0.;
		FirstOutroFrame = 0.;
		IntroBeats = 96;
		OutroBeats = 128;
		VolumeMax = 1.;
		VolumeAverage = 0.25;
		RunningPtr = running_ptr;
		
		if( AlwaysLoadFloat )
			Audio.SampleFormat = AV_SAMPLE_FMT_FLT;
		
		if( Audio.Load( filename, RunningPtr ) )
		{
			TotalFrames = Audio.Size / (Audio.BytesPerSample * Audio.Channels);
			
			if( TotalFrames / (double) Audio.SampleRate < 30. )
			{
				Audio.Clear();
				TotalFrames = 0;
			}
		}
		else
			TotalFrames = 0;
		
		if( Audio.Tags.find("title") == Audio.Tags.end() )
		{
			const char *last_slash = filename, *this_slash = NULL;
			while(( this_slash = strstr( last_slash + 1, PATH_SEPARATOR ) ))
				last_slash = this_slash;
			std::string title;
			if( strncmp( last_slash, PATH_SEPARATOR, strlen(PATH_SEPARATOR) ) == 0 )
				title = std::string( last_slash + strlen(PATH_SEPARATOR) );
			else
				title = std::string( last_slash );
			
			std::string::size_type last_dot = title.rfind('.');
			if( last_dot != std::string::npos )
				Audio.Tags["title"] = title.substr( 0, last_dot );
		}
	}
	
	const char *GetTag( const std::string &name ) const
	{
		std::map<std::string,std::string>::const_iterator tag = Audio.Tags.find(name);
		if( tag != Audio.Tags.end() )
			return tag->second.c_str();
		return NULL;
	}
	
	~Song()
	{
	}
	
	double TotalSeconds( void ) const
	{
		return TotalFrames / (double) Audio.SampleRate;
	}
	
	void SetFirstBeat( int minutes, double seconds )
	{
		FirstBeatFrame = Audio.SampleRate * ((60. * minutes) + seconds);
	}
	
	float SampleAtIndex( size_t index ) const
	{
		if(likely( index < Audio.Size / Audio.BytesPerSample ))
		{
			switch( Audio.BytesPerSample )
			{
				case 8:
				{
					double *buffer64 = (double*) Audio.Data;
					return buffer64[ index ];
				}
				case 4:
				{
					float *buffer32 = (float*) Audio.Data;
					return buffer32[ index ];
				}
				case 2:
				{
					Sint16 *buffer16 = (Sint16*) Audio.Data;
					return buffer16[ index ] / 32768.;
				}
				case 1:
				{
					Uint8 *buffer8 = (Uint8*) Audio.Data;
					return buffer8[ index ] / 128. - 1.;
				}
			}
		}
		return 0.f;
	}
	
	float SanitizedSample( size_t index ) const
	{
		// Avoid corrupted samples, which ffmpeg/libav sometimes feeds us.
		float sample = SampleAtIndex(index);
		if(likely( fabs(sample) <= 1.f ))
			return sample;
		return 0.f;
	}
	
	float NearestFrame( Uint8 channel ) const
	{
		channel %= Audio.Channels;
		size_t a_index = Audio.Channels * (size_t) CurrentFrame + channel;
		size_t b_index = a_index + Audio.Channels;
		float a = SampleAtIndex( a_index );
		float b = SampleAtIndex( b_index );
		double unused = 0.;
		return (modf( CurrentFrame, &unused ) >= 0.5) ? b : a;
	}
	
	float CubicFrame( Uint8 channel ) const
	{
		channel %= Audio.Channels;
		size_t a_index = Audio.Channels * (size_t) CurrentFrame + channel;
		size_t b_index = a_index + Audio.Channels;
		long prev_index = a_index - Audio.Channels;
		size_t next_index = b_index + Audio.Channels;
		float a = SampleAtIndex( a_index );
		float b = SampleAtIndex( b_index );
		float prev = likely(prev_index >= 0) ? SampleAtIndex( prev_index ) : 0.f;
		float next = SampleAtIndex( next_index );
		double unused = 0.;
		float b_part = modf( CurrentFrame, &unused );
		/*
		float along_a_tangent = a + b_part * (b - prev) / 2.;
		float along_b_tangent = b - (1. - b_part) * (next - a) / 2.;
		float remaining_portion = 1. - b_part * b_part - (1. - b_part) * (1. - b_part);
		float linear = a + (b - a) * b_part;
		return along_a_tangent * (1. - b_part) * (1. - b_part) + along_b_tangent * b_part * b_part + linear * remaining_portion;
		*/
		// Paul Breeuwsma came up with a simpler cubic interpolation than mine, so I'm using it.
		// http://www.paulinternet.nl/?page=bicubic
		return a + 0.5f * b_part * (b - prev + b_part * (2.f * prev - 5.f * a + 4.f * b - next + b_part * (3.f * (a - b) + next - prev)));
	}
	
	void Advance( double playback_bpm, int playback_rate )
	{
		CurrentFrame += playback_bpm * Audio.SampleRate / (BPM * playback_rate);
	}
	
	double Beat( void ) const
	{
		return BeatAtFrame( CurrentFrame );
	}
	
	double TotalBeats( void ) const
	{
		return BeatAtFrame( TotalFrames );
	}
	
	double FirstOutroBeat( void ) const
	{
		return BeatAtFrame( FirstOutroFrame );
	}
	
	double BeatAtFrame( double frame ) const
	{
		return (frame - FirstBeatFrame) * BPM / (60. * Audio.SampleRate);
	}
	
	double FrameAtBeat( double beat ) const
	{
		return beat * (60. * Audio.SampleRate) / BPM + FirstBeatFrame;
	}
	
	double NearestBeatFrameAtFrame( double frame ) const
	{
		if( ! Beats.size() )
			return frame;
		
		std::map<size_t,double>::const_iterator exact_beat_found = Beats.find( (size_t)( frame + 0.5 ) );
		if( exact_beat_found != Beats.end() )
		{
			// We found an exact match.
			
			return exact_beat_found->first;
		}
		else
		{
			// No exact match, check for beats near where we guessed.
			
			std::map<size_t,double>::const_iterator ahead = Beats.lower_bound( frame );
			std::map<size_t,double>::const_reverse_iterator behind( ahead );
			behind ++;
			
			double nearest_beat = 0.;
			double off_by = FLT_MAX;
			
			if( behind != Beats.rend() )
			{
				nearest_beat = behind->first;
				off_by = frame - nearest_beat;
			}
			
			if( (ahead != Beats.end()) && ((ahead->first - frame) < off_by) )
			{
				nearest_beat = ahead->first;
				off_by = ahead->first - frame;
			}
			
			return nearest_beat;
		}
	}
	
	double NearestBeatFrameAtBeat( double beat ) const
	{
		return NearestBeatFrameAtFrame( FrameAtBeat( beat ) );
	}
	
	double NearestBeatAtFrame( double frame ) const
	{
		return BeatAtFrame( NearestBeatFrameAtFrame( frame ) );
	}
	
	double NearestBeatAtBeat( double beat ) const
	{
		return BeatAtFrame( NearestBeatFrameAtBeat( beat ) );
	}
	
	bool Finished( void ) const
	{
		return (CurrentFrame >= TotalFrames);
	}
	
	float VolumeAdjustment( void ) const
	{
		return 0.25f / VolumeAverage;
	}
	
	float VolumeAdjustmentToMax( void ) const
	{
		return 1.f / VolumeMax;
	}
	
	void SetIntroOutroBeats( int crossfade_in = -1, int crossfade_out = -1 )
	{
		if( crossfade_in >= 0 )
			IntroBeats = crossfade_in;
		else
			IntroBeats = 96;
		
		if( crossfade_out >= 0 )
			OutroBeats = crossfade_out;
		else
			OutroBeats = 128;
		
		double total_beats = TotalBeats();
		
		if( total_beats < 8. )
		{
			IntroBeats = 1;
			OutroBeats = 1;
		}
		else if( total_beats < 16. )
		{
			IntroBeats = std::min<int>( 4, IntroBeats );
			OutroBeats = std::min<int>( 4, OutroBeats );
		}
		else if( total_beats < 32. )
		{
			IntroBeats = std::min<int>( 8, IntroBeats );
			OutroBeats = std::min<int>( 8, OutroBeats );
		}
		else if( total_beats < 64. )
		{
			IntroBeats = std::min<int>( 16, IntroBeats );
			OutroBeats = std::min<int>( 16, OutroBeats );
		}
		else if( total_beats < 128. )
		{
			IntroBeats = std::min<int>( 32, IntroBeats );
			OutroBeats = std::min<int>( 32, OutroBeats );
		}
		else if( total_beats < 256. )
		{
			IntroBeats = std::min<int>( 64, IntroBeats );
			OutroBeats = std::min<int>( 64, OutroBeats );
		}
	}
	
	void Analyze( size_t first_frame = 0 )
	{
		if( ! *RunningPtr )
			return;
		if( !( Audio.Data && Audio.Size ) )
			return;
		
		// Don't analyze more than 10 minutes of audio per song.
		size_t max_analysis_frame = std::min<size_t>( Audio.SampleRate * 60 * 10 + first_frame, TotalFrames );
		
		if( first_frame >= max_analysis_frame )
			return;
		
		// Determine song volume.
		VolumeMax = 0.;
		for( size_t frame = first_frame; frame < max_analysis_frame; frame ++ )
		{
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
			{
				float this_amp = fabs( SampleAtIndex( frame * Audio.Channels + channel ) );
				if( likely(this_amp <= 1.f) && (this_amp > VolumeMax) )
					VolumeMax = this_amp;
			}
		}
		double volume_sum = 0.;
		size_t samples_in_average = 0;
		for( size_t frame = first_frame; frame < max_analysis_frame; frame ++ )
		{
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
			{
				double this_amp = fabs( SampleAtIndex( frame * Audio.Channels + channel ) );
				if( likely(this_amp <= 1.) && (this_amp > VolumeMax / 8.) )
				{
					volume_sum += this_amp;
					samples_in_average ++;
				}
			}
		}
		if( samples_in_average )
			VolumeAverage = volume_sum / samples_in_average;
		else
			VolumeAverage = VolumeMax / 4.;
		
		SDL_Delay( 1 );
		if( ! *RunningPtr )
			return;
		
		
		// Apply low-pass filter and keep track of average peak height over a period of samples.
		
		size_t lpf_samples = std::min<size_t>( max_analysis_frame - first_frame, Audio.SampleRate * 2 / 49 ); // 44100 -> 1800
		double prev = 0., prev_abs = 0., prev_prev_abs = 0.;
		std::map<size_t,double> avg_peak;
		int num_nearby_peaks = 1;
		double highest = 0.;
		
		// Get starting block sum for the low-pass filter.
		for( size_t frame = first_frame; frame < lpf_samples + first_frame; frame ++ )
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
				prev += SanitizedSample( frame * Audio.Channels + channel ) / (double) Audio.Channels;
		
		SDL_Delay( 1 );
		if( ! *RunningPtr )
			return;
		
		// Process the rest of the audio.
		for( size_t frame = lpf_samples + first_frame; frame < max_analysis_frame; frame ++ )
		{
			// Subtract the oldest sample of the block and add the new one.
			double point = prev;
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
			{
				point -= SanitizedSample( (frame - lpf_samples) * Audio.Channels + channel ) / (double) Audio.Channels;
				point += SanitizedSample( frame * Audio.Channels + channel ) / (double) Audio.Channels;
			}
			double point_abs = fabs(point);
			
			if( (prev_abs > point_abs) && (prev_abs > prev_prev_abs) )
			{
				// Found a peak on the previous frame.
				size_t peak_frame = frame - (lpf_samples / 2) - 1;
				
				if( avg_peak.size() )
				{
					size_t last_peak = avg_peak.rbegin()->first;
					
					if( last_peak + (lpf_samples / 2) < peak_frame )
					{
						// This peak is really close to one we already detected, so add them and keep count.
						avg_peak[ last_peak ] += prev_abs;
						num_nearby_peaks ++;
					}
					else
					{
						// If the previous peak had some adjacent, now we divide by the count.
						avg_peak[ last_peak ] /= num_nearby_peaks;
						if( avg_peak[ last_peak ] > highest )
							highest = avg_peak[ last_peak ];
						
						// Start a new peak.
						avg_peak[ peak_frame ] = prev_abs;
						num_nearby_peaks = 1;
					}
				}
				else
				{
					// This is the first peak.
					avg_peak[ peak_frame ] = prev_abs;
					num_nearby_peaks = 1;
				}
			}
			
			prev_prev_abs = prev_abs;
			prev_abs = point_abs;
			prev = point;
		}
		
		if( avg_peak.size() )
		{
			size_t last_peak = avg_peak.rbegin()->first;
			
			// If the last peak had some adjacent, now we divide by the count.
			avg_peak[ last_peak ] /= num_nearby_peaks;
			if( avg_peak[ last_peak ] > highest )
				highest = avg_peak[ last_peak ];
		}
		
		SDL_Delay( 1 );
		if( ! *RunningPtr )
			return;
		
		
		// Look for significant increases in peak averages (bass hits).
		
		Beats.clear();
		size_t want_beats = max_analysis_frame / Audio.SampleRate;
		size_t prev_at = first_frame;
		for( size_t attempt = 0; (attempt < 10) && (Beats.size() < want_beats); attempt ++ )
		{
			double min_value = highest * pow( 3. / 4., ((int)(attempt / 2)) );
			double min_inc_factor = pow( 2., 7 - ((attempt + 1) / 2) );
			Beats.clear();
			prev = 0.;
			bool prev_adjacent = false;
			
			for( std::map<size_t,double>::const_iterator avg_iter = avg_peak.begin(); avg_iter != avg_peak.end(); avg_iter ++ )
			{
				bool is_near = (avg_iter->first - prev_at < lpf_samples);
				
				if( prev_adjacent && (avg_iter->second > prev) && is_near )
				{
					// The low-frequency peaks are still getting louder, and we want the crest.
					std::map<size_t,double>::iterator last = Beats.end();
					last --;
					Beats.erase( last );
					Beats[ avg_iter->first ] = avg_iter->second;
					prev_at = avg_iter->first;
				}
				else if( (avg_iter->second >= min_value) && ((! is_near) || (avg_iter->second > prev * min_inc_factor)) )
				{
					// These low-frequency peaks are enough louder than the previous, so it's the start of a bass hit.
					Beats[ avg_iter->first ] = avg_iter->second;
					prev_at = avg_iter->first;
					prev_adjacent = true;
				}
				else
					prev_adjacent = false;
				
				prev = avg_iter->second;
			}
			
			SDL_Delay( 1 );
			if( ! *RunningPtr )
				return;
			
			// No use trying again if we don't have any more low-frequency peaks detected.
			if( Beats.size() >= avg_peak.size() )
				break;
		}
		
		size_t first_beat = Beats.size() ? Beats.begin()->first : avg_peak.begin()->first;
		
		
		// Search a range of likely BPMs and see how well the beats match.
		
		double best_frame_skip = Audio.SampleRate * 60. / 138.;
		size_t best_first_beat = first_beat;
		double best_error = FLT_MAX;
		int best_doff = INT_MAX;
		
		for( int bpm = 120; (bpm <= 150) && (best_error > 0.); bpm ++ )
		{
			// Search at a BPM and see how closely it matches.
			
			double frame_skip = Audio.SampleRate * 60. / (double) bpm;
			size_t first_beat_matched = 0;
			
			// Try a few different starting points.
			for( std::map<size_t,double>::const_iterator beat_iter = Beats.begin(); (beat_iter != Beats.end()) && (beat_iter->first < max_analysis_frame / 2); beat_iter ++ )
			{
				size_t start = beat_iter->first;
				std::vector<int> doff_behind, doff_ahead;
				int off_behind = 0, off_ahead = 0;
				int misses_behind = 0, misses_ahead = 0;
				
				// Look at each place we expect the beat to hit.
				for( double frame_float = start; frame_float + 0.5 < max_analysis_frame; frame_float += frame_skip )
				{
					size_t frame = frame_float + 0.5;
					
					std::map<size_t,double>::const_iterator exact_beat_found = Beats.find( frame );
					if( exact_beat_found != Beats.end() )
					{
						// We found an exact match.
					
						if( ! misses_behind )
							doff_behind.push_back( -off_behind );
						if( ! misses_ahead )
							doff_ahead.push_back( -off_ahead );
						
						off_behind = 0;
						off_ahead = 0;
					
						misses_behind = 0;
						misses_ahead = 0;
					}
					else
					{
						// No exact match, check for beats near where we guessed.
						
						std::map<size_t,double>::const_iterator ahead = Beats.lower_bound( frame );
						std::map<size_t,double>::const_reverse_iterator behind( ahead );
						behind ++;
						
						if( (behind != Beats.rend()) && (frame - behind->first < frame_skip / 2) )
						{
							int off = behind->first - frame;
							
							if( ! misses_behind )
								doff_behind.push_back( off - off_behind );
							
							off_behind = off;
							
							misses_behind = 0;
						}
						else
							misses_behind ++;
						
						if( (ahead != Beats.end()) && (ahead->first - frame < frame_skip / 2) )
						{
							int off = ahead->first - frame;
							
							if( ! misses_ahead )
								doff_ahead.push_back( off - off_ahead );
							
							off_ahead = off;
							
							misses_ahead = 0;
						}
						else
							misses_ahead ++;
					}
					
					// Assume the first beat of the first detected pair is the first real beat.
					// NOTE: This might catch weird intro beats before the regular pattern.
					if( ! first_beat_matched )
					{
						if( doff_behind.size() )
							first_beat_matched = frame + off_behind;
						else if( doff_ahead.size() )
							first_beat_matched = frame + off_ahead;
					}
				}
				
				
				// Get the mean square errors between differences (how consistently it slipped).
				double mean_behind = 0., mean_ahead = 0.;
				double error_behind = FLT_MAX, error_ahead = FLT_MAX;
				
				if( doff_behind.size() >= 2 )
				{
					for( std::vector<int>::const_iterator doff_iter = doff_behind.begin(); doff_iter != doff_behind.end(); doff_iter ++ )
						mean_behind += *doff_iter / (double) doff_behind.size();
					
					error_behind = 0.;
					
					for( std::vector<int>::const_iterator doff_iter = doff_behind.begin(); doff_iter != doff_behind.end(); doff_iter ++ )
					{
						double ddoff = (*doff_iter - mean_behind);
						error_behind += (ddoff * ddoff) / (double) doff_behind.size();
					}
					
					// When it's close, prefer lower BPM to avoid finding every 5th beat.
					error_behind *= bpm;
				}
				
				if( doff_ahead.size() >= 2 )
				{
					for( std::vector<int>::const_iterator doff_iter = doff_ahead.begin(); doff_iter != doff_ahead.end(); doff_iter ++ )
						mean_ahead += *doff_iter / (double) doff_ahead.size();
					
					error_ahead = 0.;
					
					for( std::vector<int>::const_iterator doff_iter = doff_ahead.begin(); doff_iter != doff_ahead.end(); doff_iter ++ )
					{
						double ddoff = (*doff_iter - mean_ahead);
						error_ahead += (ddoff * ddoff) / (double) doff_ahead.size();
					}
					
					// When it's close, prefer lower BPM to avoid finding every 5th beat.
					error_ahead *= bpm;
				}
				
				// Sort differences so we can find median.
				std::sort( doff_behind.begin(), doff_behind.end() );
				std::sort( doff_ahead.begin(), doff_ahead.end() );
				
				// Find the lowest difference between differences for the closest prediction.
				if( doff_behind.size() && (error_behind < best_error) )
				{
					best_doff = doff_behind[ doff_behind.size() / 2 ];
					best_frame_skip = frame_skip + best_doff;
					best_first_beat = first_beat_matched ? first_beat_matched : first_beat;
					best_error = error_behind;
				}
				if( doff_ahead.size() && (error_ahead < best_error) )
				{
					best_doff = doff_ahead[ doff_ahead.size() / 2 ];
					best_frame_skip = frame_skip + best_doff;
					best_first_beat = first_beat_matched ? first_beat_matched : first_beat;
					best_error = error_ahead;
				}
			}
			
			SDL_Delay( 1 );
			if( ! *RunningPtr )
				return;
		}
		
		// If it's close, round to the nearest whole BPM.
		double bpm_rounding = 0.106;
		double bpm_fpart = modf( Audio.SampleRate * 60. / (double) best_frame_skip, &BPM );
		if( bpm_fpart >= (1. - bpm_rounding) )
			BPM += 1.;
		else if( bpm_fpart > bpm_rounding )
			BPM += bpm_fpart;
		
		// Sanity-check the BPM, and scale it to the expected range.
		if( BPM < 90. )
			BPM *= 2.;
		else if( BPM > 180. )
			BPM /= 2.;
		
		FirstBeatFrame = best_first_beat;
		
		SetIntroOutroBeats();
	}
};


// --------------------------------------------------------------------------------------


class SongLoader
{
public:
	volatile bool Finished;
	volatile bool Running;
	std::string Filename;
	Song *LoadingSong;
	SDL_Thread *Thread;
	
	SongLoader( void )
	{
		Finished = true;
		Running = true;
		LoadingSong = NULL;
		Thread = SDL_CreateThread( &SongLoaderThread, this );
	}
	
	~SongLoader()
	{
		Quit();
	}
	
	void StartLoading( const char *filename )
	{
		Filename.clear();
		Filename.append( filename );
		Finished = false;
	}
	
	void Quit( void )
	{
		Running = false;
		if( Thread )
		{
			int unused = 0;
			SDL_WaitThread( Thread, &unused );
			Thread = NULL;
		}
		Finished = true;
	}
};


// --------------------------------------------------------------------------------------


class EqualizerParam
{
public:
	std::map<float,float> FreqScale;
	
	EqualizerParam( void ){}
	~EqualizerParam(){}
	
	float LookupScale( float freq ) const
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
		
		// Interpolate linearly between nearest freqency volume scales.
		float b_part = (freq - prev->first) / (found->first - prev->first);
		return prev->second * (1. - b_part) + found->second * b_part;
	}
	
	float GetScale( float freq ) const
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
	
	float GetScale( size_t index, size_t frames, unsigned int rate ) const
	{
		return GetScale( ((float) index) * rate / frames );
	}
	
	float MaxScale( void ) const
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
	
	float MinScale( void ) const
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
	
	float AvgScale( void ) const
	{
		if( ! FreqScale.size() )
			return 1.f;
		
		float total = 0.f;
		for( std::map<float,float>::const_iterator iter = FreqScale.begin(); iter != FreqScale.end(); iter ++ )
			total += iter->second;
		return total / FreqScale.size();
	}
};


class EqualizerFFT
{
private:
	FFTContext *Context1, *Context2;
	FFTComplex *Complex;

public:
	unsigned int Channels, Rate;
	size_t Frames;
	
	EqualizerFFT( unsigned int channels, unsigned int rate, size_t frames )
	{
		Channels = channels;
		Rate = rate;
		Frames = frames;
		Context1 = av_fft_init( log2(Frames), false );
		Context2 = av_fft_init( log2(Frames), true );
		Complex = (FFTComplex*) av_mallocz( Frames * sizeof(FFTComplex) );
	}
	
	~EqualizerFFT()
	{
		av_fft_end( Context1 );
		av_fft_end( Context2 );
		av_free( Complex );
	}
	
	void Process( float *buffer, EqualizerParam *param, float max = 0.f )
	{
		// If clipping prevention is enabled, scale down if needed.
		float pre_eq = 1.f;
		if( max )
		{
			float highest = param->MaxScale();
			if( highest > max )
				pre_eq = max / highest;
		}
		
		for( size_t ch = 0; ch < Channels; ch ++ )
		{
			memset( Complex, 0, Frames * sizeof(FFTComplex) );
			for( size_t i = 0; i < Frames; i ++ )
				Complex[ i ].re = buffer[ (i * Channels) + ch ] * pre_eq;
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
			float old_scale = pre_eq * scale0;
			
			for( size_t i = EQ_ANTIPOP; i < Frames - EQ_ANTIPOP; i ++ )
				buffer[ (i * Channels) + ch ] = Complex[ i ].re * new_scale;
			
			size_t end_antipop = std::min<size_t>( Frames / 2, EQ_ANTIPOP );
			for( size_t i = 0; i < end_antipop; i ++ )
			{
				float new_part = i / (float) EQ_ANTIPOP;
				buffer[ (i * Channels) + ch ] = new_part * Complex[ i ].re * new_scale + (1.f - new_part) * old_scale * buffer[ (i * Channels) + ch ];
				size_t j = Frames - 1 - i;
				buffer[ (j * Channels) + ch ] = new_part * Complex[ j ].re * new_scale + (1.f - new_part) * old_scale * buffer[ (j * Channels) + ch ];
			}
		}
	}
};


class Equalizer
{
private:
	std::map<size_t,EqualizerFFT*> EqualizerFFTs;

public:
	Equalizer( void ){}
	~Equalizer(){}
	
	void Process( float *buffer, unsigned int channels, unsigned int rate, size_t frames, EqualizerParam *param, float max = 0. )
	{
		#ifdef EQ_FRAMES
		if( frames > EQ_FRAMES )
		{
			for( size_t chunk = 0; (chunk + 1) * EQ_FRAMES <= frames; chunk ++ )
				Process( (float*) (((char*) buffer) + sizeof(*buffer) * channels * chunk * EQ_FRAMES), channels, rate, EQ_FRAMES, param, max );
			return;
		}
		#endif
		
		if( param )
		{
			// NOTE: We can do this because our audio output never changes rate/channels.
			if( EqualizerFFTs.find(frames) == EqualizerFFTs.end() )
				EqualizerFFTs[ frames ] = new EqualizerFFT( channels, rate, frames );
			
			EqualizerFFTs[ frames ]->Process( buffer, param, max );
		}
	}
};


Equalizer GlobalEQ;


// --------------------------------------------------------------------------------------


class ReverbBounce
{
public:
	size_t FramesBack;
	float AmpScale;
	
	ReverbBounce( void ) { FramesBack = 0; AmpScale = 0; }
	ReverbBounce( size_t f, float a ) { FramesBack = f; AmpScale = a; }
	~ReverbBounce(){}
};


class ReverbParam
{
public:
	float SpeakerSide, SpeakerFront, SideWall, FrontWall, BackWall, Ceiling, Floor;
	float HeadWidth, BehindScale;
	float BounceEnergy;
	float TotalScale;
	
	std::vector<ReverbBounce> SameSide, OppoSide;
	
	ReverbParam( void )
	{
		SpeakerSide = 1.5;
		SpeakerFront = 0.5;
		SideWall = 3;
		FrontWall = 3;
		BackWall = 3;
		Ceiling = 2;
		Floor = 1;
		HeadWidth = 0.15;
		BehindScale = 0.5;
		BounceEnergy = 0.375;
	}
	~ReverbParam(){}
	
	float SpeakerDist( void ) const
	{
		return sqrt( SpeakerSide * SpeakerSide + SpeakerFront * SpeakerFront );
	}
	
	float BouncedDist( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite ) const
	{
		float x = SpeakerSide - HeadWidth / 2.f;
		if( x_bounces % 2 )
			x += 2 * (opposite ? SideWall : SideWall - SpeakerSide);
		x += (x_bounces / 2) * SideWall * 4;
		if( opposite && ! x_bounces )
			x += HeadWidth;
		
		float y = SpeakerFront;
		if( y_bounces % 2 )
			y += 2 * BackWall;
		y += (y_bounces / 2) * (FrontWall + BackWall) * 2;
		
		float z = 0.f;
		if( z_bounces % 2 )
			z += 2 * (up ? Ceiling : Floor);
		z += (z_bounces / 2) * (Ceiling + Floor) * 2;
		
		return sqrt( x * x + y * y + z * z ) - SpeakerDist();
	}
	
	float AmpScale( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite ) const
	{
		float bounced = BouncedDist( x_bounces, y_bounces, z_bounces, up, opposite );
		float speaker = SpeakerDist();
		float scale = pow( BounceEnergy, x_bounces + y_bounces + z_bounces ) * speaker / (bounced + speaker);
		scale *= pow( 0.75, up ? (z_bounces / 2) : ((z_bounces + 1) / 2) ); // Carpet on the floor.
		if( y_bounces % 2 ) // Arriving from behind.
			scale *= pow( BehindScale, y_bounces / (float)(x_bounces + y_bounces) );
		if( opposite && ! x_bounces )
			scale *= SpeakerFront / sqrt( SpeakerSide * SpeakerSide + SpeakerFront * SpeakerFront );
		return scale;
	}
	
	void Setup( unsigned int rate )
	{
		SameSide.clear();
		OppoSide.clear();
		
		for( int xb = 0; xb <= 5; xb ++ )
		{
			for( int yb = 0; yb <= 5; yb ++ )
			{
				if( xb || yb )
					SameSide.push_back( ReverbBounce( BouncedDist( xb, yb,  0, false, false ) * rate / 343., AmpScale( xb, yb,  0, false, false ) ) );
				OppoSide.push_back( ReverbBounce( BouncedDist( xb, yb,  0, false, true  ) * rate / 343., AmpScale( xb, yb,  0, false, true  ) ) );
				
				for( int zb = 1; zb <= 2; zb ++ )
				{
					SameSide.push_back( ReverbBounce( BouncedDist( xb, yb, zb, false, false ) * rate / 343., AmpScale( xb, yb, zb, false, false ) ) );
					SameSide.push_back( ReverbBounce( BouncedDist( xb, yb, zb, true,  false ) * rate / 343., AmpScale( xb, yb, zb, true,  false ) ) );
					OppoSide.push_back( ReverbBounce( BouncedDist( xb, yb, zb, false, true  ) * rate / 343., AmpScale( xb, yb, zb, false, true  ) ) );
					OppoSide.push_back( ReverbBounce( BouncedDist( xb, yb, zb, true,  true  ) * rate / 343., AmpScale( xb, yb, zb, true,  true  ) ) );
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
};


class Reverb
{
private:
	float *History;
	unsigned int Rate;

public:
	Reverb( void ) { History = NULL; Rate = 44100; }
	~Reverb(){}
	
	#define REVERB_HISTORY_FRAMES 262144
	
	void Process( float *buffer, unsigned int channels, unsigned int rate, size_t frames, ReverbParam *param, float max = 0.f )
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
		
		// If clipping prevention is enabled, scale down if needed.
		float scale = 1.f;
		if( max && (param->TotalScale > max) )
			scale = max / param->TotalScale;
		
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
					if( from_frame >= 0 )
						val += History[ channels * from_frame + ch ] * param->SameSide[ i ].AmpScale;
				}
				for( size_t i = 0; i < bounces_oppo; i ++ )
				{
					int from_frame = REVERB_HISTORY_FRAMES - 1 - frames_back - param->OppoSide[ i ].FramesBack;
					if( from_frame >= 0 )
						val += History[ channels * from_frame + (ch+1)%channels ] * param->OppoSide[ i ].AmpScale;
				}
				buffer[ index ] = val * scale;
			}
		}
	}
};


Reverb GlobalReverb;


// --------------------------------------------------------------------------------------


class UserData
{
public:
	bool Playing;
	bool Running;
	bool Crossfading;
	SDL_AudioSpec Spec;
	bool HighRes;
	double BPM;
	bool SourceBPM;
	double SourcePitchScale;
	int CrossfadeIn;
	int CrossfadeOut;
	float Volume;
	bool VolumeMatching;
	float VolumeLimit;
	bool Repeat;
	bool Shuffle;
	size_t ShuffleDelay;
	int Metronome;
	PlaybackBuffer Buffer;
	FILE *WriteTo;
	size_t WroteBytes;
	std::deque<std::string> Queue;
	std::deque<Song*> Songs;
	volatile Song *MostRecentToRemove;
	SongLoader Loader;
	Uint8 *VisualizerBuffer;
	int VisualizerBufferSize;
	clock_t VisualizerClock;
	char Title[ 128 ];
	char Artist[ 128 ];
	char Album[ 128 ];
	char Message[ 128 ];
	time_t MessageUntil;
	EqualizerParam *EQ;
	ReverbParam *Reverb;
	
	UserData( void )
	{
		Playing = false;
		Running = true;
		Crossfading = false;
		memset( &Spec, 0, sizeof(Spec) );
		HighRes = false;
		BPM = 140.;
		SourceBPM = true;
		SourcePitchScale = 1.;
		CrossfadeIn = 96;
		CrossfadeOut = 128;
		Volume = 1.;
		VolumeMatching = true;
		VolumeLimit = 0.;
		Repeat = true;
		Shuffle = true;
		ShuffleDelay = 0;
		Metronome = 0;
		Buffer.SetSize( 131072 );
		WriteTo = NULL;
		WroteBytes = 0;
		MostRecentToRemove = NULL;
		VisualizerBuffer = NULL;
		VisualizerBufferSize = 0;
		VisualizerClock = 0;
		memset( Title, 0, 128 );
		strcpy( Title, "Loading..." );
		memset( Artist, 0, 128 );
		memset( Album, 0, 128 );
		memset( Message, 0, 128 );
		MessageUntil = time(NULL) - 1;
		EQ = NULL;
		Reverb = NULL;
	}
	
	size_t BytesPerSample( void ) const
	{
		return (HighRes ? sizeof(float) : sizeof(Sint16));
	}
	
	void SetMessage( const char *message, int seconds )
	{
		snprintf( Message, 128, "%s", message );
		MessageUntil = time(NULL) + seconds;
	}
	
	void QueueSong( const char *filename )
	{
		Queue.push_back( filename );
	}
	
	void CheckSongLoading( void )
	{
		if( Running )
		{
			bool changed_song = false;
			
			// Remove any songs that have finished playback.
			while( MostRecentToRemove && Songs.size() )
			{
				Song *current_song = Songs.front();
				if( current_song == MostRecentToRemove )
					MostRecentToRemove = NULL;
				Songs.pop_front();
				delete current_song;
				
				changed_song = true;
			}
			
			if( Loader.Finished )
			{
				// Get the song pointer and take it away from the loader.
				Song *song = Loader.LoadingSong;
				Loader.LoadingSong = NULL;
				
				// Make sure the song loaded okay.
				if( song )
				{
					song->SetIntroOutroBeats( CrossfadeIn, CrossfadeOut );
					
					double sec = song->TotalSeconds();
					int min = ((int) sec ) / 60;
					sec -= min * 60.;
					
					printf( "%s: %.2f BPM, Length %i:%04.1f, Start %.3f sec, Beats %i, Volume Avg %.2f Max %.2f\n", Loader.Filename.c_str(), song->BPM, min, sec, song->FirstBeatFrame / song->Audio.SampleRate, (int) song->Beats.size(), song->VolumeAverage, song->VolumeMax );
					fflush( stdout );
					
					bool first_song = Songs.empty();
					
					// Only the first song should play pre-beat lead-in.
					if( ! first_song )
						song->CurrentFrame = song->FirstBeatFrame;
					
					// Done loading and analyzing, so put it in the playback queue.
					Songs.push_back( song );
					
					// Start playback after we load the first one successfully.
					if( first_song )
					{
						Playing = true;
						changed_song = true;
					}
				}
				
				// If shuffle is enabled, shuffle the queue once per loop through.
				if( Shuffle && Queue.size() && (Songs.size() < 3) )
				{
					if( ShuffleDelay )
						ShuffleDelay --;
					else
					{
						size_t dont_shuffle_end = std::min<size_t>( Queue.size(), Songs.size() );
						std::random_shuffle( Queue.begin(), Queue.end() - dont_shuffle_end );
						ShuffleDelay = Queue.size() - dont_shuffle_end - 1;
						printf( "Shuffling next %i songs.\n", (int) ShuffleDelay + 1 );
					}
				}
				
				// If we're repeating tracks, put it at the end of the queue.
				if( Repeat && Loader.Filename.size() )
					QueueSong( Loader.Filename.c_str() );
				
				Loader.Filename.clear();
				
				// Start loading the next in the queue.
				if( Queue.size() && (Songs.size() < 3) )
				{
					std::string song_name = Queue.front();
					Queue.pop_front();
					
					Loader.StartLoading( song_name.c_str() );
				}
			}
			
			// Update visualizer text.
			if( changed_song )
			{
				Title[ 0 ] = '\0';
				Artist[ 0 ] = '\0';
				Album[ 0 ] = '\0';
				if( Songs.size() )
				{
					const Song *current_song = Songs.front();
					const char *title = current_song->GetTag("title");
					if( title )
						snprintf( Title, 128, "%s", title );
					const char *artist = current_song->GetTag("artist");
					if( artist )
						snprintf( Artist, 128, "%s", artist );
					const char *album = current_song->GetTag("album");
					if( album )
						snprintf( Album, 128, "%s", album );
				}
			}
		}
		else
			Loader.Running = false;
	}
	
	float VisualizerAmpScale( void ) const
	{
		// Make the visualizer look good regardless of volume/EQ/reverb/clipping settings.
		
		float volume = Volume;
		float excess = 1.f;
		
		if( EQ )
		{
			excess *= std::max<float>( EQ->MaxScale(), EQ->GetScale(0.) );
			float eq_avg = EQ->AvgScale();
			if( eq_avg )
			{
				volume *= eq_avg;
				excess /= eq_avg;
			}
		}
		
		if( Reverb )
		{
			excess *= Reverb->TotalScale;
			float reverb_avg = pow( (float) Reverb->TotalScale, 0.3f );
			if( reverb_avg )
			{
				volume *= reverb_avg;
				excess /= reverb_avg;
			}
		}
		
		if( VolumeLimit && (volume * excess > VolumeLimit) )
			volume = VolumeLimit / excess;
		
		if( (volume <= 0.f) || (volume > 1.f) )
			return 1.f;
		
		return 1.f / volume;
	}
	
	float VisualizerSample( size_t index ) const
	{
		if( ! HighRes )
			return (((Sint16*)( VisualizerBuffer ))[ index % (VisualizerBufferSize / sizeof(Sint16)) ] / 32768.f) * VisualizerAmpScale();
		return ((float*)( VisualizerBuffer ))[ index % (VisualizerBufferSize / sizeof(float)) ] * VisualizerAmpScale();
	}
	
	~UserData(){}
};


// --------------------------------------------------------------------------------------


int SongLoaderThread( void *loader_ptr )
{
	SongLoader *loader = (SongLoader*) loader_ptr;
	
	while( loader->Running )
	{
		if( ! loader->Finished )
		{
			loader->LoadingSong = new Song( loader->Filename.c_str(), &(loader->Running) );
			
			// Make sure it loaded okay.
			if( loader->LoadingSong->Audio.Data && loader->LoadingSong->Audio.Size )
			{
				SDL_Delay( 1 );
				loader->LoadingSong->Analyze();
				
				// Retry if junk got into the start of the audio and messed up analysis.
				for( size_t retry = 0; (retry < 4) && (loader->LoadingSong->Beats.size() < 10); retry ++ )
				{
					SDL_Delay( 1 );
					loader->LoadingSong->Analyze( loader->LoadingSong->FirstBeatFrame + 400 );
				}
			}
			else
			{
				// It didn't load correctly, so we won't queue it.  Memory cleanup.
				delete loader->LoadingSong;
				loader->LoadingSong = NULL;
			}
			
			// Notify main thread that it should now take posession of LoadingSong.
			loader->Finished = true;
		}
		
		// The loader can wait 5sec each time.
		SDL_Delay( loader->Running ? 5000 : 1 );
	}
	
	return 0;
}


// --------------------------------------------------------------------------------------


float LinearCrossfade( float a, float b, float b_percent )
{
	return a * (1.f - b_percent) + b * b_percent;
}


float EqualPowerCrossfade( float a, float b, float b_percent )
{
	return a * cos( b_percent * M_PI * 0.5f ) + b * cos( (1.f - b_percent) * M_PI * 0.5f );
}


Sint16 FloatTo16( float value )
{
	value *= 32768.f;
	if( value >= 0.f )
	{
		if(likely( value < 32767.f ))
			return value + 0.5f;
		else
			return 32767;
	}
	else
	{
		if(likely( value > -32768.f ))
			return value - 0.5f;
		else
			return -32768;
	}
}


// --------------------------------------------------------------------------------------


bool CalculateCrossfade( const Song *current_song, const Song *next_song, double *crossfade_for_beats, double *crossfade_at_beat )
{
	if( !( current_song && next_song ) )
		return false;
	
	double current_total_beats = current_song->TotalBeats();
	
	// Use the minimum of the intro and outro beats for crossfade length.
	*crossfade_for_beats = std::min<int>( current_song->OutroBeats, next_song->IntroBeats );
	
	if( current_song->FirstOutroFrame )
		*crossfade_at_beat = current_song->BeatAtFrame( current_song->FirstOutroFrame );
	else if( *crossfade_for_beats < 1. )
	{
		// No crossfade means don't bother lining up beats.
		*crossfade_at_beat = current_total_beats - *crossfade_for_beats;
		return true;
	}
	else
	{
		// If FirstOutroFrame is not specified, guess a good crossfade point.
		double expected_beat = current_total_beats - fmod( current_total_beats, 16. ) - current_song->OutroBeats;
		double actual_beat = current_song->NearestBeatAtBeat( expected_beat );
		*crossfade_at_beat = (fabs( actual_beat - expected_beat ) < 32.) ? actual_beat : expected_beat;
	}
	
	// Sanity check to avoid something strange like negative crossfade.
	if( *crossfade_at_beat > current_total_beats )
		*crossfade_at_beat = current_total_beats;
	
	// Make sure the current song doesn't end before the cross-fade is complete.
	if( *crossfade_for_beats > current_total_beats - *crossfade_at_beat )
		*crossfade_for_beats = current_total_beats - *crossfade_at_beat;
	
	return true;
}


// --------------------------------------------------------------------------------------


void AudioCallback( void *userdata, Uint8* stream, int len )
{
	if( ! len )
		return;
	memset( stream, 0, len );
	
	UserData *ud = (UserData*) userdata;
	
	// If we're not filling an intermediate buffer and playback is paused, leave it silent.
	if( !( ud->Playing || ud->Buffer.BufferSize ) )
		return;
	
	size_t songs_size = ud->Songs.size();
	if( ! songs_size )
		return;
	
	size_t samples = len / sizeof(Sint16);
#ifdef ALLOW_VLA
	float streamTemp[ ud->HighRes ? 1 : samples ];
#endif

	float *streamF = NULL;
	if( ! ud->HighRes )
	{
		// Create an intermediate floating-point buffer for audio processing.
#ifdef ALLOW_VLA
		streamF = streamTemp;
#else
		streamF = (float*) malloc( samples * sizeof(float) );
#endif
		memset( streamF, 0, samples * sizeof(float) );
	}
	else
	{
		// Write directly to high-resolution floating-point audio stream.
		streamF = (float*) stream;
		samples = len / sizeof(float);
	}
	
	bool calculated_crossfade = false;
	double crossfade_for_beats = 96.;
	double crossfade_at_beat = 0.;
	bool crossfade_now = false;
	float crossfade = 0.;
	float volume = 1.;
	float volume_next = 1.;
	
	Song *current_song = ud->Songs.front();
	Song *next_song = (songs_size >= 2) ? ud->Songs.at( 1 ) : NULL;
	
	double bpm = ud->SourceBPM ? (current_song->BPM * ud->SourcePitchScale) : ud->BPM;
	
	volume = ud->VolumeMatching ? ud->Volume * current_song->VolumeAdjustment() : ud->Volume;
	volume_next = (ud->VolumeMatching && next_song) ? ud->Volume * next_song->VolumeAdjustment() : ud->Volume;
	if( ud->VolumeLimit )
	{
		volume = std::min<float>( current_song->VolumeAdjustmentToMax() * ud->VolumeLimit, volume );
		if( next_song )
			volume_next = std::min<float>( next_song->VolumeAdjustmentToMax() * ud->VolumeLimit, volume_next );
	}
	
	for( size_t i = 0; i < samples; i += ud->Spec.channels )
	{
		// Check for crossfade into second song.
		
		crossfade_now = false;
		crossfade = 0.;
		
		if( next_song )
		{
			// Calculate when to do the crossfade and for how long.
			if( ! calculated_crossfade )
				calculated_crossfade = CalculateCrossfade( current_song, next_song, &crossfade_for_beats, &crossfade_at_beat );
			
			// Now that we've determined how crossfade should be done, see if we should do it now.
			if( calculated_crossfade && (current_song->Beat() >= crossfade_at_beat) )
			{
				crossfade_now = true;
				
				if( crossfade_for_beats > 0. )
				{
					crossfade = (current_song->Beat() - crossfade_at_beat) / crossfade_for_beats;
					if( crossfade > 1. )
						crossfade = 1.;
				}
				else
					crossfade = 1.;
			}
		}
		
		
		if( ! crossfade_now )
		{
			// Not crossfading yet.
			
			if( ((size_t) ud->Spec.freq == current_song->Audio.SampleRate) && (fabs(current_song->BPM - bpm) < 0.01) )
			{
				for( int channel = 0; channel < ud->Spec.channels; channel ++ )
					streamF[ i + channel ] = current_song->NearestFrame( channel ) * volume;
			}
			else
			{
				for( int channel = 0; channel < ud->Spec.channels; channel ++ )
					streamF[ i + channel ] = current_song->CubicFrame( channel ) * volume;
			}
			
			current_song->Advance( bpm, ud->Spec.freq );
		}
		else
		{
			// Crossfading now.
			
			if( ud->SourceBPM )
				bpm = LinearCrossfade( current_song->BPM * ud->SourcePitchScale, next_song->BPM * ud->SourcePitchScale, crossfade );
			
			bool a_nearest = ( ((size_t) ud->Spec.freq == current_song->Audio.SampleRate) && (fabs(current_song->BPM - bpm) < 0.01) );
			bool b_nearest = ( ((size_t) ud->Spec.freq == next_song->Audio.SampleRate   ) && (fabs(next_song->BPM    - bpm) < 0.01) );
			
			for( int channel = 0; channel < ud->Spec.channels; channel ++ )
			{
				float a = (a_nearest ? current_song->NearestFrame( channel ) : current_song->CubicFrame( channel )) * volume;
				float b = (b_nearest ? next_song->NearestFrame( channel ) : next_song->CubicFrame( channel )) * volume_next;
				streamF[ i + channel ] = EqualPowerCrossfade( a, b, crossfade );
			}
			
			current_song->Advance( bpm, ud->Spec.freq );
			next_song->Advance( bpm, ud->Spec.freq );
			
			// If we completed a crossfade, treat the next song as current.
			// Leave a note for the main thread that there's a song to remove.
			if(unlikely( crossfade >= 1. ))
			{
				ud->MostRecentToRemove = current_song;
				
				current_song = next_song;
				next_song = (songs_size >= 3) ? ud->Songs.at( 2 ) : NULL;
				
				volume = volume_next;
				volume_next = (ud->VolumeMatching && next_song) ? ud->Volume * next_song->VolumeAdjustment() : ud->Volume;
				if( ud->VolumeLimit && next_song )
					volume_next = std::min<float>( next_song->VolumeAdjustmentToMax() * ud->VolumeLimit, volume_next );
				
				calculated_crossfade = false;
				crossfade_now = false;
				crossfade = 0.;
			}
		}
		
		ud->Crossfading = crossfade_now;
		
		
		// Add metronome if enabled.
		
		if( ud->Metronome )
		{
			double beat_ipart = 0.;
			double beat_fpart = modf( current_song->Beat(), &beat_ipart );
			
			if( beat_fpart < 0. )
			{
				beat_fpart += 1.;
				beat_ipart -= 2.; // -1, and then another to fix the rounding direction.
			}
			
			if( beat_fpart < 0.25 )
			{
				int ipart = (((int)( beat_ipart + 0.5 )) % 4 + 4) % 4;
				float hz = ipart ? 960. : 1920.;
				float wave_value = (sin( beat_fpart * hz ) * 0.25 * cos( beat_fpart * M_PI * 2. ));
				
				for( int channel = 0; channel < ud->Spec.channels; channel ++ )
				{
					if( (ud->Metronome == 2) || (ipart % ud->Spec.channels == channel) )
						streamF[ i + channel ] = streamF[ i + channel ] + wave_value;
				}
			}
		}
	}
	
	// Determine max volume for post effects (EQ/reverb).
	float max = 0.f;
	if( ud->VolumeLimit )
	{
		float vol = volume * (1.f - crossfade) + volume_next * crossfade;
		max = vol ? (ud->VolumeLimit / vol) : 0.f;
	}
		
	// Apply equalizer if enabled.
	if( ud->EQ )
	{
		GlobalEQ.Process( streamF, ud->Spec.channels, ud->Spec.freq, samples / ud->Spec.channels, ud->EQ, max );
		if( max )
		{
			float highest = ud->EQ->MaxScale();
			if( highest > max )
				max = 1.f;  // The EQ already hit our volume cap, so don't let reverb go louder.
			else
				max /= highest;
		}
	}
		
	// Apply reverb if enabled, or just remember old audio buffer.
	GlobalReverb.Process( streamF, ud->Spec.channels, ud->Spec.freq, samples / ud->Spec.channels, ud->Reverb, max );
	
	// If not outputting float, convert to 16-bit output format.
	if( (void*) streamF != stream )
	{
		Sint16 *stream16 = (Sint16*) stream;
		for( size_t i = 0; i < samples; i ++ )
			stream16[ i ] = FloatTo16( streamF[ i ] );
		
#ifndef ALLOW_VLA
		free( streamF );
		streamF = NULL;
#endif
	}
	
	// If we're supposed to write our output to a file, do it after each buffer we fill.
	if( ud->WriteTo )
	{
		fwrite( stream, 1, len, ud->WriteTo );
		ud->WroteBytes += len;
	}
	
	// Remove tracks that have finished playing.
	while( ud->Songs.size() && (ud->Songs.front()->Finished()) )
	{
		Song *front = ud->Songs.front();
		ud->Songs.pop_front();
		delete front;
	}
}


void BufferedAudioCallback( void *userdata, Uint8 *stream, int len )
{
	UserData *ud = (UserData*) userdata;
	if( ud->Playing )
	{
		int buffered = ud->Buffer.Buffered;
		if( buffered < len )
		{
			// Prevent PlaybackBuffer from calling AudioCallback directly.
			// This was needed because AudioCallback wasn't fully thread-safe when songs ended.
			// It's probably not necessary to do this check anymore.
			ud->Buffer.FillStream( userdata, stream, buffered );
			memset( stream + buffered, 0, len - buffered );
		}
		else
			ud->Buffer.FillStream( userdata, stream, len );
	}
	else
		memset( stream, 0, len );
	
	UpdateVisualizerBuffer( ud, stream, len );
}


void UnbufferedAudioCallback( void *userdata, Uint8 *stream, int len )
{
	AudioCallback( userdata, stream, len );
	
	UpdateVisualizerBuffer( (UserData*) userdata, stream, len );
}


void UpdateVisualizerBuffer( UserData *ud, Uint8 *stream, int len )
{
	memcpy( ud->VisualizerBuffer, stream, std::min<int>( len, ud->VisualizerBufferSize ) );
	ud->VisualizerClock = clock();
}


// --------------------------------------------------------------------------------------


std::deque<std::string> DirSongs( std::string path )
{
	std::deque<std::string> songs;
	
	struct stat stat_result;
	memset( &stat_result, 0, sizeof(stat_result) );
	int error = stat( path.c_str(), &stat_result );
	
	if( S_ISDIR( stat_result.st_mode ) )
	{
		DIR *dir = opendir( path.c_str() );
		if( dir )
		{
			struct dirent *entry = NULL;
			
			while( (entry = readdir(dir)) )
			{
				if( entry->d_name[ 0 ] == '.' )
				{
					// Ignore current directory (.) parent directory (..) and hidden files (.*).
				}
				else if( ! (entry->d_type & DT_DIR) )
				{
					// Not a directory.
					
					songs.push_back( std::string(path) + std::string(PATH_SEPARATOR) + std::string(entry->d_name) );
				}
				else
				{
					// Subdirectory.
					
					std::deque<std::string> subdir_songs = DirSongs( path + std::string(PATH_SEPARATOR) + std::string(entry->d_name) );
					for( std::deque<std::string>::const_iterator song_iter = subdir_songs.begin(); song_iter != subdir_songs.end(); song_iter ++ )
						songs.push_back( *song_iter );
				}
			}
			
			closedir( dir );
		}
	}
	else if( ! error )
	{
		const char *path_str = path.c_str();
		size_t path_len = strlen(path_str);
		if( (path_len > 4) && (strncasecmp( path_str + strlen(path_str) - 4, ".m3u", 4 ) == 0) )
		{
			// Playlist file.
			
			FILE *file = fopen( path_str, "rt" );
			if( file )
			{
				std::string dir_path;
				std::string::size_type last_slash = path.rfind(PATH_SEPARATOR[0]);
				if( last_slash != std::string::npos )
					dir_path = path.substr( 0, last_slash ) + std::string(PATH_SEPARATOR);
				
				char buffer[ 32768 ] = "";
				while( ! feof(file) )
				{
					buffer[ 0 ] = '\0';
					fgets( buffer, 32768, file );
					char *comment = strchr( buffer, '#' );
					if( comment )
						comment[ 0 ] = '\0';
					size_t len = strlen(buffer);
					if( ! len )
						continue;
					if( buffer[ len - 1 ] == '\n' )
					{
						buffer[ len - 1 ] = '\0';
						len --;
						if( ! len )
							continue;
					}
					if( buffer[ len - 1 ] == '\r' )
					{
						buffer[ len - 1 ] = '\0';
						len --;
						if( ! len )
							continue;
					}
					
					// Determine if this playlist item is a relative or absolute path.
					std::string item_path;
#ifdef WIN32
					// Windows paths usually look like "C:\Music" instead of "\Music".
					const char *colon = strchr( buffer, ':' );
					const char *slash = strstr( buffer, PATH_SEPARATOR );
					if( colon && ((colon < slash) || ! slash) )
						item_path = buffer;
					else
#endif
					if( strncmp( buffer, PATH_SEPARATOR, strlen(PATH_SEPARATOR) ) == 0 )
						item_path = buffer;
					else
						item_path = dir_path + std::string(buffer);
					
					std::deque<std::string> subdir_songs = DirSongs( item_path );
					for( std::deque<std::string>::const_iterator song_iter = subdir_songs.begin(); song_iter != subdir_songs.end(); song_iter ++ )
						songs.push_back( *song_iter );
				}
				
				fclose( file );
			}
		}
		else
			// Individual song.
			songs.push_back( path );
	}
	
	return songs;
}


// --------------------------------------------------------------------------------------


#ifdef WAVEOUT

HWAVEOUT WaveOutHandle = NULL;
WAVEHDR WaveOutHeaders[ 2 ];
Uint8 *WaveOutBuffers[ 2 ] = { NULL, NULL };
int WaveOutBufferSize = 0;
int WaveOutBufferNum = 0;
int WaveOutBuffersNeeded = 0;

void CALLBACK WaveOutCallback( HWAVEOUT hwo, UINT uMsg, DWORD_PTR dwInstance, DWORD_PTR dwParam1, DWORD_PTR dwParam2 )
{
	if( uMsg == WOM_DONE )
		WaveOutBuffersNeeded = std::min<int>( 2, WaveOutBuffersNeeded + 1 );
}

void WaveOutCheck( UserData *ud )
{
	while( WaveOutBuffersNeeded > 0 )
	{
		WaveOutBuffersNeeded --;
		
		// Fill the next buffer.
		ud->Spec.callback( ud, WaveOutBuffers[ WaveOutBufferNum ], WaveOutBufferSize );
		
		// Start playing the buffer.
		waveOutWrite( WaveOutHandle, &(WaveOutHeaders[ WaveOutBufferNum ]), sizeof(WAVEHDR) );
		
		// Swap buffers.
		WaveOutBufferNum ++;
		WaveOutBufferNum %= 2;
		
		// Point the visualizer to the other buffer, because that one is playing sooner.
		ud->VisualizerBuffer = WaveOutBuffers[ WaveOutBufferNum ];
		ud->VisualizerBufferSize = WaveOutBufferSize;
		ud->VisualizerClock = clock();
	}
}

#endif


// --------------------------------------------------------------------------------------


void UpdatePlayback( UserData *ud )
{
	if( ud->Buffer.Unfilled() )
		ud->Buffer.AddToBuffer( ud, ud->Buffer.BufferSize );
	
#ifdef WAVEOUT
	WaveOutCheck( ud );
#endif
}


int PlaybackThread( void *userdata )
{
	UserData *ud = (UserData*) userdata;
	
	while( ud->Running )
	{
		UpdatePlayback( ud );
		SDL_Delay( 1 );
	}
	
	return 0;
}


// --------------------------------------------------------------------------------------


Uint32 cycle_color( double cycle, const SDL_PixelFormat *format )
{
	double unused = 0.;
	cycle = modf( cycle, &unused );
	double leg = cycle * 8.;
	double along = modf( leg, &unused );
	Uint8 r = 0xFF, g = 0x00, b = 0x00;
	
	switch( (int) leg )
	{
		case 0: // Red to Yellow
			r = 255;
			g = along * 255. + 0.5;
			b = 0;
			break;
		case 1: // Yellow
			r = 255;
			g = 255;
			b = 0;
			break;
		case 2: // Yellow to Green
			r = (1. - along) * 255. + 0.5;
			g = 255;
			b = 0;
			break;
		case 3: // Green
			r = 0;
			g = 255;
			b = 0;
			break;
		case 4: // Green to Cyan
			r = 0;
			g = 255;
			b = along * 255. + 0.5;
			break;
		case 5: // Cyan to White
			r = along * 255. + 0.5;
			g = 255;
			b = 255;
			break;
		case 6: // White to Magenta
			r = 255;
			g = (1. - along) * 255. + 0.5;
			b = 255;
			break;
		case 7: // Magenta to Red
			r = 255;
			g = 0;
			b = (1. - along) * 255. + 0.5;
			break;
	}
	
	return SDL_MapRGB( format, r, g, b );
}


// --------------------------------------------------------------------------------------


int main( int argc, char **argv )
{
	std::vector<const char*> paths;
	UserData userdata;
	bool window = true;
	bool fullscreen = false;
	bool resize = false;
#ifdef WIN32
	resize = true;
#endif
	int zoom = 1;
	int visualizer = 2;
	bool playback = true;
	bool sdl_audio = true;
	bool buffer1auto = true, buffer2auto = true;
#ifdef WAVEOUT
	// Default to WaveOut floating-point when possible.
	sdl_audio = false;
	userdata.HighRes = true;
#endif
	int visualizer_color1 = 6, visualizer_color2 = 4, visualizer_text_color = 0;
	const char *write = NULL;
	SDL_AudioSpec want;
	memset( &want, 0, sizeof(want) );
	want.freq = 44100;
	want.format = AUDIO_S16;
	want.channels = 2;
	want.samples = 4096;
	want.callback = BufferedAudioCallback;
	want.userdata = &userdata;
	FontDraw font1(FONT1), font2(FONT2);
	
	// Process command-line arguments.
	for( int i = 1; i < argc; i ++ )
	{
		if( strncmp( argv[ i ], "--", 2 ) == 0 )
		{
			if( strcasecmp( argv[ i ], "--no-window" ) == 0 )
				window = false;
			else if( strcasecmp( argv[ i ], "--fullscreen" ) == 0 )
				fullscreen = true;
			else if( strcasecmp( argv[ i ], "--resize" ) == 0 )
				resize = true;
			else if( strcasecmp( argv[ i ], "--no-resize" ) == 0 )
				resize = false;
			else if( strncasecmp( argv[ i ], "--zoom=", strlen("--zoom=") ) == 0 )
				zoom = std::max<int>( 1, atoi( argv[ i ] + strlen("--zoom=") ) );
			else if( strncasecmp( argv[ i ], "--visualizer=", strlen("--visualizer=") ) == 0 )
				visualizer = atoi( argv[ i ] + strlen("--visualizer=") );
			else if( strcasecmp( argv[ i ], "--rainbow" ) == 0 )
			{
				visualizer_color1 = 8;
				visualizer_color2 = 9;
				visualizer_text_color = 8;
			}
			else if( strcasecmp( argv[ i ], "--no-shuffle" ) == 0 )
				userdata.Shuffle = false;
			else if( strcasecmp( argv[ i ], "--no-repeat" ) == 0 )
				userdata.Repeat = false;
			else if( strcasecmp( argv[ i ], "--metronome" ) == 0 )
				userdata.Metronome = 1;
			else if( strncasecmp( argv[ i ], "--pitch=", strlen("--pitch=") ) == 0 )
			{
				const char *pitch = argv[ i ] + strlen("--pitch=");
				if( pitch[ 0 ] == '+' )
					pitch ++;
				userdata.SourcePitchScale = pow( 2., atof(pitch) / 12. );
				userdata.SourceBPM = true;
			}
			else if( strncasecmp( argv[ i ], "--bpm=", strlen("--bpm=") ) == 0 )
			{
				userdata.BPM = atof( argv[ i ] + strlen("--bpm=") );
				userdata.SourceBPM = false;
			}
			else if( strncasecmp( argv[ i ], "--crossfade-in=", strlen("--crossfade-in=") ) == 0 )
				userdata.CrossfadeIn = atoi( argv[ i ] + strlen("--crossfade-in=") );
			else if( strncasecmp( argv[ i ], "--crossfade-out=", strlen("--crossfade-out=") ) == 0 )
				userdata.CrossfadeOut = atoi( argv[ i ] + strlen("--crossfade-out=") );
			else if( strncasecmp( argv[ i ], "--crossfade=", strlen("--crossfade=") ) == 0 )
			{
				userdata.CrossfadeIn = atoi( argv[ i ] + strlen("--crossfade=") );
				userdata.CrossfadeOut = userdata.CrossfadeIn;
			}
			else if( strncasecmp( argv[ i ], "--volume=", strlen("--volume=") ) == 0 )
			{
				const char *volume = argv[ i ] + strlen("--volume=");
				switch( volume[ 0 ] )
				{
					case '+':
						volume ++;
						// fallthrough
					case '-':
						userdata.Volume = pow( 2., atof(volume) / 6. );
						break;
					default:
						userdata.Volume = atof( volume );
				}
			}
			else if( strcasecmp( argv[ i ], "--no-volume-matching" ) == 0 )
				userdata.VolumeMatching = false;
			else if( strcasecmp( argv[ i ], "--prevent-clipping" ) == 0 )
				userdata.VolumeLimit = 1.;
			else if( strncasecmp( argv[ i ], "--eq=", strlen("--eq=") ) == 0 )
			{
				userdata.EQ = new EqualizerParam();
				const char *eq = argv[ i ] + strlen("--eq=");
				while( eq && eq[ 0 ] )
				{
					float freq = atof( eq );
					const char *colon = strchr( eq, ':' );
					if( colon )
					{
						const char *db = colon + 1;
						if( db[ 0 ] == '+' )
							db ++;
						userdata.EQ->FreqScale[ freq ] = pow( 2., atof(db) / 6. );
					}
					eq = strchr( eq, ',' );
					if( eq )
						eq ++;
				}
			}
			else if( strcasecmp( argv[ i ], "--reverb" ) == 0 )
			{
				userdata.Reverb = new ReverbParam();
				userdata.Reverb->Setup( 44100 );
			}
			else if( strncasecmp( argv[ i ], "--reverb=", strlen("--reverb=") ) == 0 )
			{
				userdata.Reverb = new ReverbParam();
				userdata.Reverb->BounceEnergy = atof( argv[ i ] + strlen("--reverb=") );
				userdata.Reverb->Setup( 44100 );
			}
			else if( strcasecmp( argv[ i ], "--sdl-audio" ) == 0 )
				sdl_audio = true;
			else if( strcasecmp( argv[ i ], "--16-bit" ) == 0 )
				userdata.HighRes = false;
			else if( strcasecmp( argv[ i ], "--float" ) == 0 )
				userdata.HighRes = true;
			else if( strcasecmp( argv[ i ], "--load-float" ) == 0 )
				AlwaysLoadFloat = true;
			else if( strncasecmp( argv[ i ], "--rate=", strlen("--rate=") ) == 0 )
				want.freq = atoi( argv[ i ] + strlen("--rate=") );
			else if( strncasecmp( argv[ i ], "--channels=", strlen("--channels=") ) == 0 )
				want.channels = atoi( argv[ i ] + strlen("--channels=") );
			else if( strncasecmp( argv[ i ], "--buffer1=", strlen("--buffer1=") ) == 0 )
			{
				want.samples = atoi( argv[ i ] + strlen("--buffer1=") );
				buffer1auto = false;
			}
			else if( strncasecmp( argv[ i ], "--buffer2=", strlen("--buffer2=") ) == 0 )
			{
				userdata.Buffer.SetSize( atoi( argv[ i ] + strlen("--buffer2=") ) * 2 * want.channels );
				buffer2auto = false;
			}
			else if( strcasecmp( argv[ i ], "--no-playback" ) == 0 )
				playback = false;
			else if( strcasecmp( argv[ i ], "--compile" ) == 0 )
			{
				playback = false;
				userdata.Buffer.SetSize( 0 );
				buffer2auto = false;
				window = false;
				userdata.Repeat = false;
				userdata.Shuffle = false;
			}
			else if( strncasecmp( argv[ i ], "--write=", strlen("--write=") ) == 0 )
				write = argv[ i ] + strlen("--write=");
			else
				fprintf( stderr, "Unknown option: %s\n", argv[ i ] );
		}
		else if( strlen( argv[ i ] ) )
		{
			char *path = argv[ i ];
			
			// Remove trailing slashes from directories.
			size_t len = strlen(path);
			while( (path[ len - 1 ] == '/') || (path[ len - 1 ] == '\\') )
			{
				path[ len - 1 ] = '\0';
				len --;
			}
			
			paths.push_back( path );
		}
	}
	
	if( ! paths.size() )
	{
		// If no music was dropped onto the icon or specified on the command-line, use default paths.
		#ifdef WIN32
		paths.push_back( "M:\\iTunes\\iTunes Music\\Trance" );
		#else
		const char *home = getenv("HOME");
		if( home )
		{
			int path_size = strlen(home) + strlen("/Music/iTunes/iTunes Music/Trance") + 1;
			char *path = (char*) alloca( path_size );
			snprintf( path, path_size, "%s/Music/iTunes/iTunes Music/Trance", home );
			paths.push_back( path );
		}
		paths.push_back( "/Volumes/Media/Music/iTunes/iTunes Music/Trance" );
		#endif
	}
	
	for( size_t i = 0; i < paths.size(); i ++ )
	{
		// Get all playable files from within the selected directory, or add selected songs directly.
		std::deque<std::string> songs = DirSongs( paths[ i ] );
		for( std::deque<std::string>::const_iterator song_iter = songs.begin(); song_iter != songs.end(); song_iter ++ )
			userdata.QueueSong( (*song_iter).c_str() );
	}
	
	// Seed the random number generator (for shuffle).
	srand( time(NULL) );
	
	// If we want to capture input events, we'll need to initialize SDL video.
	Uint32 sdl_flags = SDL_INIT_AUDIO;
	if( window )
		sdl_flags |= SDL_INIT_VIDEO;
	
	// Initialize SDL.
	SDL_Init( sdl_flags );
	SDL_Surface *screen = NULL;
	if( window )
	{
		SDL_WM_SetCaption( "Raptor007's AutoDJ", "AutoDJ" );
		screen = fullscreen ? SDL_SetVideoMode( 0, 0, 0, SDL_SWSURFACE | SDL_FULLSCREEN ) : SDL_SetVideoMode( 256*zoom, 64*zoom, 0, SDL_SWSURFACE | (resize ? SDL_RESIZABLE : 0) );
	}
	
	// Prepare audio file input.
	av_log_set_level( AV_LOG_FATAL );
	av_register_all();
	
	// Automatic buffer1 is number of samples, just based on sample rate.
	if( buffer1auto )
		want.samples = 1 << (int) log2(want.freq * 0.095);
	
	// Prepare audio output.
	userdata.Buffer.Callback = AudioCallback;
	if( userdata.Buffer.BufferSize > 0 )
		want.callback = BufferedAudioCallback;
	else
		want.callback = UnbufferedAudioCallback;
	memcpy( &(userdata.Spec), &want, sizeof(want) );
	
	if( playback )
	{
#ifdef WAVEOUT
		if( ! sdl_audio )
		{
			// We will fall back to SDL_audio unless WaveOut succeeds.
			sdl_audio = true;
			
			int wave_out_devices = waveOutGetNumDevs();
			if( wave_out_devices )
			{
				// Start WaveOut.
				WAVEFORMATEX wfx;
				memset( &wfx, 0, sizeof(WAVEFORMATEX) );
				wfx.nChannels = want.channels;
				wfx.nSamplesPerSec = want.freq;
				wfx.cbSize = 0;
				MMRESULT wave_out_result = ~MMSYSERR_NOERROR;
				if( userdata.HighRes )
				{
					wfx.wFormatTag = WAVE_FORMAT_IEEE_FLOAT;
					wfx.wBitsPerSample = 32;
					wfx.nBlockAlign = wfx.nChannels * wfx.wBitsPerSample / 8;
					wfx.nAvgBytesPerSec = wfx.nBlockAlign * wfx.nSamplesPerSec;
					wave_out_result = waveOutOpen( &WaveOutHandle, WAVE_MAPPER, &wfx, (DWORD_PTR) &WaveOutCallback, 0, CALLBACK_FUNCTION );
				}
				if( wave_out_result != MMSYSERR_NOERROR )
				{
					wfx.wFormatTag = WAVE_FORMAT_PCM;
					wfx.wBitsPerSample = 16;
					wfx.nBlockAlign = wfx.nChannels * wfx.wBitsPerSample / 8;
					wfx.nAvgBytesPerSec = wfx.nBlockAlign * wfx.nSamplesPerSec;
					wave_out_result = waveOutOpen( &WaveOutHandle, WAVE_MAPPER, &wfx, (DWORD_PTR) &WaveOutCallback, 0, CALLBACK_FUNCTION );
					if( wave_out_result == MMSYSERR_NOERROR )
						userdata.HighRes = false;
				}
				if( wave_out_result == MMSYSERR_NOERROR )
				{
					// Allocate 2 output buffers.
					WaveOutBufferSize = want.samples * wfx.nChannels * wfx.wBitsPerSample / 8;
					WaveOutBuffers[ 0 ] = (Uint8*) malloc( WaveOutBufferSize );
					WaveOutBuffers[ 1 ] = (Uint8*) malloc( WaveOutBufferSize );
					memset( WaveOutBuffers[ 0 ], 0, WaveOutBufferSize );
					memset( WaveOutBuffers[ 1 ], 0, WaveOutBufferSize );
					memset( &(WaveOutHeaders[ 0 ]), 0, sizeof(WAVEHDR) );
					memset( &(WaveOutHeaders[ 1 ]), 0, sizeof(WAVEHDR) );
					WaveOutHeaders[ 0 ].lpData = (HPSTR) WaveOutBuffers[ 0 ];
					WaveOutHeaders[ 1 ].lpData = (HPSTR) WaveOutBuffers[ 1 ];
					WaveOutHeaders[ 0 ].dwBufferLength = WaveOutBufferSize;
					WaveOutHeaders[ 1 ].dwBufferLength = WaveOutBufferSize;
					int error = waveOutPrepareHeader( WaveOutHandle, &(WaveOutHeaders[ 0 ]), sizeof(WAVEHDR) );
					if( ! error )
						error = waveOutPrepareHeader( WaveOutHandle, &(WaveOutHeaders[ 1 ]), sizeof(WAVEHDR) );
					
					if( ! error )
					{
						// WaveOut is ready to go with 2 output buffers, so don't use SDL_audio.
						WaveOutBuffersNeeded = 2;
						sdl_audio = false;
						
						if( screen )
						{
							userdata.VisualizerBufferSize = WaveOutBufferSize;
							userdata.VisualizerBuffer = WaveOutBuffers[ 0 ];
						}
					}
				}
			}
		}
#endif
		if( sdl_audio )
		{
			userdata.HighRes = false;
			SDL_OpenAudio( &want, &(userdata.Spec) );
			
			if( screen )
			{
				userdata.VisualizerBufferSize = userdata.Spec.samples * userdata.Spec.channels * userdata.BytesPerSample();
				userdata.VisualizerBuffer = (Uint8*) malloc( userdata.VisualizerBufferSize );
			}
		}
	}
	
	size_t bytes_per_sample = userdata.BytesPerSample();
	
	// Automatic buffer2 is number of bytes, based on data rate.
	if( buffer2auto )
		userdata.Buffer.SetSize( 4 * userdata.Spec.samples * userdata.Spec.channels * bytes_per_sample );
	
	if( write )
	{
		// Write WAV header.
		unsigned char wave_header[ 44 ] = { 'R','I','F','F', 36,0xFF,0xFF,0x7F, 'W','A','V','E', 'f','m','t',' ', 16,0,0,0, (userdata.HighRes?3:1),0, userdata.Spec.channels,0, userdata.Spec.freq%256,userdata.Spec.freq/256,0,0, (userdata.Spec.freq*userdata.Spec.channels*bytes_per_sample)%256,(userdata.Spec.freq*userdata.Spec.channels*bytes_per_sample)/256,0,0, userdata.Spec.channels*bytes_per_sample,0, 8*bytes_per_sample,0, 'd','a','t','a', 0,0xFF,0xFF,0x7F };
		userdata.WriteTo = fopen( write, "wb" );
		fwrite( wave_header, 1, 44, userdata.WriteTo );
		fflush( userdata.WriteTo );
	}
	
	// Keep track of equalizer parameters for toggling.
	EqualizerParam *disabled_eq = userdata.EQ;
	if( ! disabled_eq )
	{
		// Start with a flat 10-band EQ if user didn't specify.
		disabled_eq = new EqualizerParam();
		disabled_eq->FreqScale[    32. ] = 1.;
		disabled_eq->FreqScale[    64. ] = 1.;
		disabled_eq->FreqScale[   125. ] = 1.;
		disabled_eq->FreqScale[   250. ] = 1.;
		disabled_eq->FreqScale[   500. ] = 1.;
		disabled_eq->FreqScale[  1000. ] = 1.;
		disabled_eq->FreqScale[  2000. ] = 1.;
		disabled_eq->FreqScale[  4000. ] = 1.;
		disabled_eq->FreqScale[  8000. ] = 1.;
		disabled_eq->FreqScale[ 16000. ] = 1.;
	}
	
	// Keep track of disabled reverb parameters for toggling.
	ReverbParam *disabled_reverb = userdata.Reverb;
	if( ! disabled_reverb )
		disabled_reverb = new ReverbParam();
	
	// Visualizer variables.
	SDL_Surface *drawto = screen;
	if( zoom > 1 )
		drawto = SDL_CreateRGBSurface( SDL_SWSURFACE, (screen->w+zoom-1)/zoom, (screen->h+zoom-1)/zoom, screen->format->BitsPerPixel,
			screen->format->Rmask, screen->format->Gmask, screen->format->Bmask, screen->format->Amask );
	clock_t visualizer_updated_clock = 0, visualizer_prev_clock = 0;
	size_t visualizer_start_sample = 0;
	size_t visualizer_loading_frame = 0;
	int visualizer_color_cycle = 128;
	int visualizer_backgrounds[ VISUALIZERS ] = { 0, 2, 3, 3, 1 };
	int visualizer_frames = userdata.VisualizerBufferSize / (userdata.Spec.channels * bytes_per_sample);
	int visualizer_fft_frames = std::min<int>( visualizer_frames, 1 << (int) log2(userdata.Spec.freq * 0.095 / 4.) );
	FFTContext *visualizer_fft_context = av_fft_init( log2(visualizer_fft_frames), false );
	FFTComplex *visualizer_fft_complex_l = (FFTComplex*) av_mallocz( visualizer_fft_frames * sizeof(FFTComplex) );
	FFTComplex *visualizer_fft_complex_r = (FFTComplex*) av_mallocz( visualizer_fft_frames * sizeof(FFTComplex) );
	int visualizer_fft_width_offset = 0;
	int visualizer_fft_height_offset = 0;
	char visualizer_message[ 128 ] = "";
	int visualizer_scoot = 0;
	
	// Keep running until playback is complete.
	userdata.Running = (userdata.Queue.size() || userdata.Songs.size());
	SDL_Thread *playback_thread = SDL_CreateThread( &PlaybackThread, &userdata );
	while( userdata.Running )
	{
		size_t prev_song_count = userdata.Songs.size();
		
		userdata.CheckSongLoading();
		
		if( sdl_audio && userdata.Playing && userdata.Songs.size() && ! prev_song_count )
			SDL_PauseAudio( 0 );
		
		if( window )
		{
			SDL_Event event;
			while( SDL_PollEvent( &event ) )
			{
				if( event.type == SDL_QUIT )
				{
					userdata.Running = false;
					if( sdl_audio )
						SDL_PauseAudio( 1 );
				}
				else if( event.type == SDL_VIDEORESIZE )
				{
					if( (event.resize.w != screen->w) || (event.resize.h != screen->h) )
					{
						int scale = std::max<int>( 1, screen->w / drawto->w );
						SDL_FreeSurface( screen );
						screen = SDL_SetVideoMode( (event.resize.w/scale)*scale, (event.resize.h/scale)*scale, 0, SDL_SWSURFACE | SDL_RESIZABLE );
						if( scale <= 1 )
							drawto = screen;
						else
						{
							SDL_FreeSurface( drawto );
							drawto = SDL_CreateRGBSurface( SDL_SWSURFACE, screen->w/scale, screen->h/scale, screen->format->BitsPerPixel,
								screen->format->Rmask, screen->format->Gmask, screen->format->Bmask, screen->format->Amask );
						}
					}
				}
				else if( event.type == SDL_KEYDOWN )
				{
					SDLKey key = event.key.keysym.sym;
					bool shift = event.key.keysym.mod & (KMOD_LSHIFT | KMOD_RSHIFT);
					if( key == SDLK_q )
					{
						userdata.Running = false;
						if( sdl_audio )
							SDL_PauseAudio( 1 );
					}
					else if( key == SDLK_SPACE )
					{
						userdata.Playing = ! userdata.Playing;
						if( sdl_audio )
							SDL_PauseAudio( userdata.Playing ? 0 : 1 );
					}
					else if( key == SDLK_PAGEUP )
					{
						int scale = (screen->w / drawto->w) + 1;
						if( ! fullscreen )
						{
							if( screen == drawto )
							{
								drawto = SDL_CreateRGBSurface( SDL_SWSURFACE, screen->w, screen->h, screen->format->BitsPerPixel,
									screen->format->Rmask, screen->format->Gmask, screen->format->Bmask, screen->format->Amask );
							}
							SDL_FreeSurface( screen );
							screen = SDL_SetVideoMode( drawto->w*scale, drawto->h*scale, 0, SDL_SWSURFACE | (resize ? SDL_RESIZABLE : 0) );
						}
						else
						{
							SDL_FreeSurface( drawto );
							drawto = SDL_CreateRGBSurface( SDL_SWSURFACE, (screen->w+scale-1)/scale, (screen->h+scale-1)/scale, screen->format->BitsPerPixel,
								screen->format->Rmask, screen->format->Gmask, screen->format->Bmask, screen->format->Amask );
						}
					}
					else if( key == SDLK_PAGEDOWN )
					{
						int scale = (screen->w / drawto->w) - 1;
						if( scale >= 1 )
						{
							if( ! fullscreen )
							{
								SDL_FreeSurface( screen );
								screen = SDL_SetVideoMode( drawto->w*scale, drawto->h*scale, 0, SDL_SWSURFACE | (resize ? SDL_RESIZABLE : 0) );
							}
							if( scale == 1 )
							{
								SDL_FreeSurface( drawto );
								drawto = screen;
							}
							else if( fullscreen )
							{
								SDL_FreeSurface( drawto );
								drawto = SDL_CreateRGBSurface( SDL_SWSURFACE, (screen->w+scale-1)/scale, (screen->h+scale-1)/scale, screen->format->BitsPerPixel,
									screen->format->Rmask, screen->format->Gmask, screen->format->Bmask, screen->format->Amask );
							}
						}
					}
					else if( key == SDLK_MINUS )
					{
						userdata.Volume *= pow( 2., shift ? -3./6. : -1./6. );
						float db = 6. * log2(userdata.Volume);
						snprintf( visualizer_message, 128, "Volume: %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_EQUALS )
					{
						userdata.Volume *= pow( 2., shift ? 3./6. : 1./6. );
						float db = 6. * log2(userdata.Volume);
						snprintf( visualizer_message, 128, "Volume: %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_BACKSPACE )
					{
						userdata.VolumeMatching = ! userdata.VolumeMatching;
						snprintf( visualizer_message, 128, "Volume Matching: %s\n", userdata.VolumeMatching ? "On" : "Off" );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_p )
					{
						if( ! userdata.VolumeLimit )
							userdata.VolumeLimit = shift ? 8. : 1.;
						else if( userdata.VolumeLimit < 1.1 )
							userdata.VolumeLimit = shift ? 0. : sqrt(2.);
						else if( userdata.VolumeLimit > 7.9 )
							userdata.VolumeLimit *= shift ? sqrt(0.5) : 0.;
						else
							userdata.VolumeLimit *= shift ? sqrt(0.5) : sqrt(2.);
						
						if( userdata.VolumeLimit )
						{
							float db = 6. * log2(userdata.VolumeLimit);
							snprintf( visualizer_message, 128, "Clipping Limit: +%.0fdB\n", db );
						}
						else
							snprintf( visualizer_message, 128, "Clipping Limit: Off\n" );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_LEFTBRACKET || key == SDLK_LEFT )
					{
						if( userdata.Songs.size() )
						{
							Song *current_song = userdata.Songs.front();
							current_song->CurrentFrame = 0.;
							current_song->FirstOutroFrame = 0.;
							current_song->SetIntroOutroBeats( userdata.CrossfadeIn, userdata.CrossfadeOut );
						}
					}
					else if( key == SDLK_RIGHTBRACKET || key == SDLK_RIGHT )
					{
						if( userdata.Songs.size() )
						{
							Song *current_song = userdata.Songs.front();
							int beat = current_song->Beat();
							if( userdata.Crossfading )
							{
								current_song->FirstOutroFrame = 1.;
								current_song->OutroBeats = 1;
							}
							else
							{
								int rewind_to_major_beat = beat % 16;
								double expected_beat = beat - rewind_to_major_beat;
								if( rewind_to_major_beat >= 8 )
									expected_beat += 16.;
								double actual_beat = current_song->NearestBeatAtBeat( expected_beat );
								current_song->FirstOutroFrame = current_song->FrameAtBeat( (fabs( actual_beat - expected_beat ) < 32.) ? actual_beat : expected_beat );
								current_song->OutroBeats = std::min<int>( 64, current_song->TotalBeats() - beat );
							}
							if( userdata.Songs.size() < 2 )
								userdata.SetMessage( "Loading...", 4 );
						}
					}
					else if( key == SDLK_BACKSLASH )
					{
						if( userdata.Songs.size() )
						{
							userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->TotalFrames - 2;
							if( userdata.Songs.size() < 2 )
								userdata.SetMessage( "Loading...", 4 );
						}
					}
					else if( key == SDLK_k )
					{
						if( userdata.SourceBPM && userdata.Songs.size() )
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
						
						userdata.BPM -= 1.;
						userdata.SourceBPM = false;
						snprintf( visualizer_message, 128, "Playback BPM: %.0f\n", userdata.BPM );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_l )
					{
						if( userdata.SourceBPM && userdata.Songs.size() )
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
						
						userdata.BPM += 1.;
						userdata.SourceBPM = false;
						snprintf( visualizer_message, 128, "Playback BPM: %.0f\n", userdata.BPM );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_SEMICOLON || key == SDLK_DOWN )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale /= pow( 2., 1./12. );
						double st = log2(userdata.SourcePitchScale) / log2(pow( 2., 1./12. ));
						snprintf( visualizer_message, 128, "Pitch: %s%.0f semitone%s\n", (st >= 0.) ? "+" : "", st, ((int)(fabs(st)+0.5) == 1) ? "" : "s" );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_QUOTE || key == SDLK_UP )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale *= pow( 2., 1./12. );
						double st = log2(userdata.SourcePitchScale) / log2(pow( 2., 1./12. ));
						snprintf( visualizer_message, 128, "Pitch: %s%.0f semitone%s\n", (st >= 0.) ? "+" : "", st, ((int)(fabs(st)+0.5) == 1) ? "" : "s" );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
					}
					else if( key == SDLK_i )
					{
						float source = userdata.Songs.size() ? userdata.Songs.front()->BPM : 0.;
						float playback = userdata.SourceBPM ? source * userdata.SourcePitchScale : userdata.BPM;
						snprintf( visualizer_message, 128, "Song %.1f BPM, Playback %.1f BPM\n", source, playback );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_x )
					{
						if( userdata.CrossfadeIn < 96 )
						{
							userdata.CrossfadeIn = 96;
							userdata.CrossfadeOut = 128;
						}
						else
						{
							userdata.CrossfadeIn = 1;
							userdata.CrossfadeOut = 1;
						}
						
						for( std::deque<Song*>::iterator song_iter = userdata.Songs.begin(); song_iter != userdata.Songs.end(); song_iter ++ )
							(*song_iter)->SetIntroOutroBeats( userdata.CrossfadeIn, userdata.CrossfadeOut );
						
						int beats = std::min<int>( userdata.CrossfadeIn, userdata.CrossfadeOut );
						snprintf( visualizer_message, 128, "Crossfade: %i Beat%s\n", beats, (beats == 1) ? "" : "s" );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_m )
					{
						if( userdata.Metronome )
							userdata.Metronome = 0;
						else if( shift )
							userdata.Metronome = 2;
						else
							userdata.Metronome = 1;
					}
					else if( key == SDLK_COMMA )
					{
						if( userdata.Songs.size() )
						{
							Song *current_song = userdata.Songs.front();
							current_song->CurrentFrame = current_song->FrameAtBeat( current_song->Beat() - 32. );
							if( current_song->CurrentFrame < 0. )
								current_song->CurrentFrame = 0.;
						}
					}
					else if( key == SDLK_PERIOD )
					{
						if( userdata.Songs.size() )
						{
							Song *current_song = userdata.Songs.front();
							current_song->CurrentFrame = current_song->FrameAtBeat( current_song->Beat() + 32. );
						}
					}
					else if( key == SDLK_SLASH )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale = 1.;
						if( userdata.Songs.size() )
						{
							userdata.BPM = userdata.Songs.front()->BPM;
							snprintf( visualizer_message, 128, "Pitch/BPM: Match Source (%.1f BPM)\n", userdata.BPM );
						}
						else
							snprintf( visualizer_message, 128, "Pitch/BPM: Match Source\n" );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_v )
					{
						if( screen )
						{
							if( ! shift )
								visualizer ++;
							else
								visualizer += VISUALIZERS - 1;
							visualizer %= VISUALIZERS;
						}
						else
							visualizer = 0;
						
						if( screen && ! visualizer )
						{
							Uint32 bg = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0xFF );
							SDL_LockSurface( screen );
							Uint32* pixels = (Uint32*) screen->pixels;
							for( int x = 0; x < screen->w; x ++ )
							{
								for( int y = 0; y < screen->h; y ++ )
									pixels[ y * screen->w + x ] = bg;
							}
							SDL_UnlockSurface( screen );
							SDL_UpdateRect( screen, 0, 0, screen->w, screen->h );
						}
					}
					else if( key == SDLK_f )
					{
						if( shift )
							visualizer_color1 += VISUALIZER_COLORS - 1;
						else
							visualizer_color1 ++;
						visualizer_color1 %= VISUALIZER_COLORS;
					}
					else if( key == SDLK_g )
					{
						if( shift )
							visualizer_color2 += VISUALIZER_COLORS - 1;
						else
							visualizer_color2 ++;
						visualizer_color2 %= VISUALIZER_COLORS;
					}
					else if( key == SDLK_c )
					{
						if( ! shift )
							visualizer_color_cycle /= 2;
						else
							visualizer_color_cycle *= 2;
						if( visualizer_color_cycle < 32 )
							visualizer_color_cycle = 512;
						else if( visualizer_color_cycle > 512 )
							visualizer_color_cycle = 32;
					}
					else if( key == SDLK_d )
					{
						if( ! shift )
							visualizer_text_color ++;
						else
							visualizer_text_color += VISUALIZER_COLORS - 1;
						visualizer_text_color %= VISUALIZER_COLORS;
					}
					else if( key == SDLK_b )
					{
						if( visualizer && userdata.Playing )
						{
							if( ! shift )
								visualizer_backgrounds[ visualizer ] ++;
							else
								visualizer_backgrounds[ visualizer ] += VISUALIZER_BACKGROUNDS - 1;
							visualizer_backgrounds[ visualizer ] %= VISUALIZER_BACKGROUNDS;
						}
					}
					else if( key == SDLK_e )
					{
						if( shift || ! userdata.EQ )
						{
							if( ! userdata.EQ )
								userdata.EQ = disabled_eq;
							
							bool was_flat = (userdata.EQ->MaxScale() == 1.) && (userdata.EQ->MinScale() == 1.);
							
							if( shift || was_flat )
							{
								userdata.EQ->FreqScale.clear();
								userdata.EQ->FreqScale[    32. ] = 1.;
								userdata.EQ->FreqScale[    64. ] = was_flat ? pow( 2.,  1./6. ) : 1.;
								userdata.EQ->FreqScale[   125. ] = 1.;
								userdata.EQ->FreqScale[   250. ] = 1.;
								userdata.EQ->FreqScale[   500. ] = 1.;
								userdata.EQ->FreqScale[  1000. ] = was_flat ? pow( 2., -1./6. ) : 1.;
								userdata.EQ->FreqScale[  2000. ] = was_flat ? pow( 2., -1./6. ) : 1.;
								userdata.EQ->FreqScale[  4000. ] = was_flat ? pow( 2., -1./6. ) : 1.;
								userdata.EQ->FreqScale[  8000. ] = was_flat ? pow( 2., -2./6. ) : 1.;
								userdata.EQ->FreqScale[ 16000. ] = was_flat ? pow( 2., -1./6. ) : 1.;
								
								snprintf( visualizer_message, 128, "Equalizer: %s\n", was_flat ? "Earbuds" : "Flat" );
							}
							else
								snprintf( visualizer_message, 128, "Equalizer: On\n" );
						}
						else
						{
							disabled_eq = userdata.EQ;
							userdata.EQ = NULL;
							
							snprintf( visualizer_message, 128, "Equalizer: Off\n" );
						}
						
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_1 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 32. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 32. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 32Hz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_2 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 64. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 64. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 64Hz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_3 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 125. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 125. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 125Hz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_4 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 250. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 250. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 250Hz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_5 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 500. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 500. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 500Hz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_6 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 1000. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 1000. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 1KHz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_7 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 2000. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 2000. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 2KHz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_8 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 4000. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 4000. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 4KHz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_9 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 8000. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 8000. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 8KHz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_0 )
					{
						if( ! userdata.EQ )
							userdata.EQ = disabled_eq;
						
						float scale = userdata.EQ->GetScale( 16000. ) * pow( 2., shift ? -1./6. : 1./6. );
						userdata.EQ->FreqScale[ 16000. ] = scale;
						float db = 6. * log2(scale);
						
						snprintf( visualizer_message, 128, "Equalizer: 16KHz %s%.0fdB\n", (db >= 0.) ? "+" : "", db );
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
					else if( key == SDLK_r )
					{
						if( ! userdata.Reverb )
							userdata.Reverb = disabled_reverb;
						else if( ! shift )
						{
							disabled_reverb = userdata.Reverb;
							userdata.Reverb = NULL;
						}
						else
						{
							userdata.Reverb->BounceEnergy += 0.0625;
							if( userdata.Reverb->BounceEnergy > 1. )
								userdata.Reverb->BounceEnergy = 0.125;
							userdata.Reverb->Setup( userdata.Spec.freq );
						}
						
						if( userdata.Reverb )
							snprintf( visualizer_message, 128, "Reverb Bounce: %0.1f%%\n", userdata.Reverb->BounceEnergy * 100. );
						else
							snprintf( visualizer_message, 128, "Reverb: Off\n" );
						
						userdata.SetMessage( visualizer_message, 4 );
						printf( "%s", visualizer_message );
						fflush( stdout );
					}
				}
			}
		}
		
		// If the main thread is handling playback, do that now.
		if( ! playback_thread )
			UpdatePlayback( &userdata );
		
		// Wait 10ms.
		SDL_Delay( 10 );
		
		// Prepare to draw.
		Uint32* pixels = NULL;
		if( screen && drawto )
		{
			if( screen == drawto )
				SDL_LockSurface( screen );
			pixels = (Uint32*) drawto->pixels;
		}
		
		// Visualize the waveform in the window.
		if( pixels && visualizer && (userdata.Playing || userdata.Songs.empty()) )
		{
			clock_t now = clock();
			
			// Keep the visualizer synchronized with the audio playback.
			if( userdata.Playing )
			{
				if( userdata.VisualizerClock != visualizer_updated_clock )
				{
					visualizer_updated_clock = userdata.VisualizerClock;
					visualizer_start_sample = 0;
					visualizer_prev_clock = now;
				}
				else
				{
					size_t frames = userdata.Spec.freq * (now - visualizer_prev_clock) / (double) CLOCKS_PER_SEC;
					visualizer_prev_clock = now;
					visualizer_start_sample += userdata.Spec.channels * frames;
					visualizer_start_sample %= userdata.VisualizerBufferSize / (userdata.Spec.channels * userdata.BytesPerSample());
				}
			}
			else
			{
				visualizer_loading_frame ++;
				visualizer_loading_frame %= drawto->h;
				visualizer_prev_clock = now;
			}
			
			Uint32 colors[ VISUALIZER_COLORS ];
			colors[ 0 ] = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0xFF ); // White to Blue
			colors[ 1 ] = SDL_MapRGB( screen->format, 0x00, 0xFF, 0xFF ); // Cyan to Blue
			colors[ 2 ] = SDL_MapRGB( screen->format, 0x00, 0x00, 0xFF ); // Blue
			colors[ 3 ] = SDL_MapRGB( screen->format, 0xFF, 0x00, 0xFF ); // Magenta to Blue
			colors[ 4 ] = SDL_MapRGB( screen->format, 0xFF, 0x00, 0x00 ); // Red
			colors[ 5 ] = SDL_MapRGB( screen->format, 0xFF, 0x80, 0x00 ); // Orange to Red
			colors[ 6 ] = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0x00 ); // Yellow to Red
			colors[ 7 ] = SDL_MapRGB( screen->format, 0x00, 0xFF, 0x00 ); // Green
			double cycle = (now / (double) CLOCKS_PER_SEC) * (userdata.BPM / 60.) / (double) visualizer_color_cycle;
			colors[ 8 ] = cycle_color( cycle, screen->format );
			colors[ 9 ] = cycle_color( cycle - 0.03125, screen->format );
			
			font1.SetTarget( pixels, drawto->w, drawto->h );
			font2.SetTarget( pixels, drawto->w, drawto->h );
			
			#define FADE_R(r,g,b) (std::max<int>( 0, std::min<int>( 0xFF, (r - b) * 0.875f ) ))
			#define FADE_G(r,g,b) (std::min<int>( 0xFF, g * 0.5f + std::max<float>( 0.f, (g - r - b) * 0.25f ) ))
			#define FADE_B(r,g,b) (std::min<int>( 0xFF, b * 0.9375f ))
			
			if( (visualizer_backgrounds[ visualizer ] == 0) || ! userdata.Playing )
			{
				// Fade in-place.
				for( int x = 0; x < drawto->w; x ++ )
				{
					for( int y = 0; y < drawto->h; y ++ )
					{
						Uint8 r = 0x00, g = 0x00, b = 0x00;
						SDL_GetRGB( pixels[ y * drawto->w + x ], screen->format, &r, &g, &b );
						pixels[ y * drawto->w + x ] = SDL_MapRGB( screen->format, FADE_R(r,g,b), FADE_G(r,g,b), FADE_B(r,g,b) );
					}
				}
			}
			else if( visualizer_backgrounds[ visualizer ] == 1 )
			{
				// Fade downward.
				for( int x = 0; x < drawto->w; x ++ )
				{
					Uint8 r = 0x00, g = 0x00, b = 0x00;
					for( int y = drawto->h - 1; y >= 1; y -- )
					{
						Uint8 r1 = 0x00, g1 = 0x00, b1 = 0x00, r2 = 0x00, g2 = 0x00, b2 = 0x00;
						SDL_GetRGB( pixels[  y    * drawto->w + x ], screen->format, &r1, &g1, &b1 );
						SDL_GetRGB( pixels[ (y-1) * drawto->w + x ], screen->format, &r2, &g2, &b2 );
						r = std::max<Uint8>( r1, r2 );
						g = std::max<Uint8>( g1, g2 );
						b = std::max<Uint8>( b1, b2 );
						pixels[ y * drawto->w + x ] = SDL_MapRGB( screen->format, FADE_R(r,g,b), FADE_G(r,g,b), FADE_B(r,g,b) );
					}
					SDL_GetRGB( pixels[ x ], screen->format, &r, &g, &b );
					pixels[ x ] = SDL_MapRGB( screen->format, FADE_R(r,g,b), FADE_G(r,g,b), FADE_B(r,g,b) );
				}
			}
			else if( visualizer_backgrounds[ visualizer ] == 2 )
			{
				// Fade outward.
				for( int x = 0; x < (drawto->w + 1) / 2; x ++ )
				{
					Uint8 r = 0x00, g = 0x00, b = 0x00;
					for( int y = 0; y < (drawto->h + 1) / 2; y ++ )
					{
						for( int q = 0; q < 4; q ++ )
						{
							int x1 = x, y1 = y, dx = 1, dy = 1;
							if( q >= 2 )
							{
								x1 = drawto->w - x - 1;
								dx = -1;
							}
							if( q % 2 )
							{
								y1 = drawto->h - y - 1;
								dy = -1;
							}
							Uint8 r1 = 0x00, g1 = 0x00, b1 = 0x00, r2 = 0x00, g2 = 0x00, b2 = 0x00;
							SDL_GetRGB( pixels[  y1     * drawto->w + x1    ], screen->format, &r1, &g1, &b1 );
							SDL_GetRGB( pixels[ (y1+dy) * drawto->w + x1+dx ], screen->format, &r2, &g2, &b2 );
							r = std::max<Uint8>( r1, r2 );
							g = std::max<Uint8>( g1, g2 );
							b = std::max<Uint8>( b1, b2 );
							pixels[ y1 * drawto->w + x1 ] = SDL_MapRGB( screen->format, FADE_R(r,g,b), FADE_G(r,g,b), FADE_B(r,g,b) );
						}
					}
				}
			}
			else if( visualizer_backgrounds[ visualizer ] == 3 )
			{
				// Fade in all directions.
				Uint8 r = 0x00, g = 0x00, b = 0x00;
				Uint8 r1 = 0x00, g1 = 0x00, b1 = 0x00, r2 = 0x00, g2 = 0x00, b2 = 0x00, r3 = 0x00, g3 = 0x00, b3 = 0x00;
				for( int y = 0; y < drawto->h; y ++ )
				{
					r2 = 0x00; g2 = 0x00; b2 = 0x00;
					SDL_GetRGB( pixels[ y * drawto->w ], screen->format, &r3, &g3, &b3 );
					for( int x = 0; x < drawto->w; x ++ )
					{
						r1 = r2; g1 = g2; b1 = b2;
						r2 = r3; g2 = g3; b2 = b3;
						if( x + 1 < drawto->w )
							SDL_GetRGB( pixels[ y * drawto->w + (x+1) ], screen->format, &r3, &g3, &b3 );
						else
							{ r3 = 0x00; g3 = 0x00; b3 = 0x00; }
						r = std::max<Uint8>( std::max<Uint8>( FADE_R(r1,g1,b1), r2 ), FADE_R(r3,g3,b3) );
						g = std::max<Uint8>( std::max<Uint8>( FADE_G(r1,g1,b1), g2 ), FADE_G(r3,g3,b3) );
						b = std::max<Uint8>( std::max<Uint8>( FADE_B(r1,g1,b1), b2 ), FADE_B(r3,g3,b3) );
						pixels[ y * drawto->w + x ] = SDL_MapRGB( screen->format, r, g, b );
					}
				}
				for( int x = 0; x < drawto->w; x ++ )
				{
					r2 = 0x00; g2 = 0x00; b2 = 0x00;
					SDL_GetRGB( pixels[ x ], screen->format, &r3, &g3, &b3 );
					for( int y = 0; y < drawto->h; y ++ )
					{
						r1 = r2; g1 = g2; b1 = b2;
						r2 = r3; g2 = g3; b2 = b3;
						if( y + 1 < drawto->h )
							SDL_GetRGB( pixels[ (y+1) * drawto->w + x ], screen->format, &r3, &g3, &b3 );
						else
							{ r3 = 0x00; g3 = 0x00; b3 = 0x00; }
						r = std::max<Uint8>( std::max<Uint8>( r1, r2 ), r3 );
						g = std::max<Uint8>( std::max<Uint8>( g1, g2 ), g3 );
						b = std::max<Uint8>( std::max<Uint8>( b1, b2 ), b3 );
						pixels[ y * drawto->w + x ] = SDL_MapRGB( screen->format, FADE_R(r,g,b), FADE_G(r,g,b), FADE_B(r,g,b) );
					}
				}
			}
			
			if( ! userdata.Playing )
			{
				// Loading animation.
				for( int x = 0; x < drawto->w; x ++ )
				{
					int y1 = ((drawto->h - x) % drawto->h + drawto->h + visualizer_loading_frame) % drawto->h;
					int y2 = (y1 + drawto->h / 2) % drawto->h;
					pixels[ y1 * drawto->w + x ] = colors[ visualizer_color1 ];
					pixels[ y2 * drawto->w + x ] = colors[ visualizer_color2 ];
				}
			}
			
			else if( visualizer == 1 )
			{
				// Fill in the waveform.
				for( int x = 0; x < drawto->w; x ++ )
				{
					for( int ch = userdata.Spec.channels - 1; ch >= 0; ch -- )
					{
						float amplitude = userdata.VisualizerSample( visualizer_start_sample + (x * userdata.Spec.channels) + ch );
						int y = std::max<int>( 0, std::min<int>( drawto->h - 1, drawto->h * (0.5f - amplitude * 0.5f) + 0.5f ) );
						pixels[ y * drawto->w + x ] = (ch % 2) ? colors[ visualizer_color2 ] : colors[ visualizer_color1 ];
					}
				}
			}
			
			else if( visualizer == 2 )
			{
				// Show spectrum of frequencies.
				int frames = std::min<int>( visualizer_fft_frames, visualizer_frames - visualizer_start_sample / userdata.Spec.channels );
				float base = pow( 2., log2( visualizer_fft_frames/2 ) / (drawto->w + visualizer_fft_width_offset) );
				for( int ch = userdata.Spec.channels - 1; ch >= 0; ch -- )
				{
					FFTComplex *visualizer_fft_complex = (ch % 2) ? visualizer_fft_complex_r : visualizer_fft_complex_l;
					memset( visualizer_fft_complex, 0, visualizer_fft_frames * sizeof(FFTComplex) );
					for( int i = 0; i < frames; i ++ )
						visualizer_fft_complex[ i ].re = userdata.VisualizerSample( visualizer_start_sample + (i * userdata.Spec.channels) + ch );
					av_fft_permute( visualizer_fft_context, visualizer_fft_complex );
					av_fft_calc( visualizer_fft_context, visualizer_fft_complex );
					int offset = 0;
					for( int x = 0; x < drawto->w; x ++ )
					{
						int band_min = std::min<int>( visualizer_fft_frames - 1, offset + pow( base, x ) );
						int band_max = std::min<int>( visualizer_fft_frames, offset + pow( base, x + 1 ) );
						if( band_min == band_max )
						{
							offset ++;
							band_max ++;
						}
						float harmonics = 1.f + log2( (visualizer_fft_frames/2) / band_max );
						float amplitude = 0.f;
						for( int band = band_min; band < band_max; band ++ )
							amplitude += sqrt( visualizer_fft_complex[ band ].re * visualizer_fft_complex[ band ].re + visualizer_fft_complex[ band ].im * visualizer_fft_complex[ band ].im );
						amplitude /= std::max<float>( 7.f, pow(2.f,harmonics) );
						int h = std::max<int>( 0, std::min<int>( drawto->h - 1, drawto->h * amplitude ) );
						for( int y = drawto->h - 1 - h; y < drawto->h; y ++ )
							pixels[ y * drawto->w + x ] = (ch % 2) ? colors[ visualizer_color2 ] : colors[ visualizer_color1 ];
					}
					visualizer_fft_width_offset = offset;
				}
			}
			
			else if( visualizer == 3 )
			{
				// Show left/right spectrum of frequencies.
				int frames = std::min<int>( visualizer_fft_frames, visualizer_frames - visualizer_start_sample / userdata.Spec.channels );
				float base = pow( 2., log2( visualizer_fft_frames*0.4375 ) / (drawto->h + visualizer_fft_height_offset) );
				for( int ch = (userdata.Spec.channels > 1) ? 1 : 0; ch >= 0; ch -- )
				{
					FFTComplex *visualizer_fft_complex = (ch % 2) ? visualizer_fft_complex_r : visualizer_fft_complex_l;
					memset( visualizer_fft_complex, 0, visualizer_fft_frames * sizeof(FFTComplex) );
					for( int i = 0; i < frames; i ++ )
						for( int in_ch = ch; in_ch < userdata.Spec.channels; in_ch += 2 )
							visualizer_fft_complex[ i ].re += userdata.VisualizerSample( visualizer_start_sample + (i * userdata.Spec.channels) + in_ch );
					av_fft_permute( visualizer_fft_context, visualizer_fft_complex );
					av_fft_calc( visualizer_fft_context, visualizer_fft_complex );
				}
				const FFTComplex *spectrum_l = visualizer_fft_complex_l;
				const FFTComplex *spectrum_r = (userdata.Spec.channels == 1) ? visualizer_fft_complex_l : visualizer_fft_complex_r;
				int offset = 0;
				for( int y = 0; y < drawto->h; y ++ )
				{
					int band_min = std::min<int>( visualizer_fft_frames - 1, offset + pow( base, y ) );
					int band_max = std::min<int>( visualizer_fft_frames, offset + pow( base, y + 1 ) );
					if( band_min == band_max )
					{
						offset ++;
						band_max ++;
					}
					float harmonics = 1.f + log2( (visualizer_fft_frames/2) / band_max );
					float amplitude_l = 0.f, amplitude_r = 0.f;
					for( int band = band_min; band < band_max; band ++ )
					{
						amplitude_l += sqrt( spectrum_l[ band ].re * spectrum_l[ band ].re + spectrum_l[ band ].im * spectrum_l[ band ].im );
						amplitude_r += sqrt( spectrum_r[ band ].re * spectrum_r[ band ].re + spectrum_r[ band ].im * spectrum_r[ band ].im );
					}
					float amplitude = amplitude_l + amplitude_r;
					if( ! amplitude )
						continue;
					float right_ratio = amplitude_r / amplitude;
					float separation = fabs( 0.5f - right_ratio ) * 2.f;
					amplitude /= std::max<float>( 7.f, pow(2.f,harmonics) );
					for( int ch = userdata.Spec.channels - 1; ch >= 0; ch -- )
					{
						float this_amplitude = (ch % 2) ? amplitude_r : amplitude_l;
						Uint32 color = (separation * separation * this_amplitude * 2.f > 0.015625f) ? colors[ visualizer_color1 ] : colors[ visualizer_color2 ];
						this_amplitude /= std::max<float>( 7.f, pow(2.f,harmonics) );
						int h = std::max<int>( 0, std::min<int>( drawto->w / 2 - 1, drawto->w / 2 * this_amplitude ) );
						int dx = (ch % 2) ? 1 : -1;
						for( int x = drawto->w / 2 - (ch+1)%2; (ch % 2) ? (x < drawto->w / 2 + h) : (x > drawto->w / 2 - h); x += dx )
							pixels[ (drawto->h - y - 1) * drawto->w + x ] = color;
					}
				}
				visualizer_fft_height_offset = offset;
			}
			
			else if( visualizer == 4 )
			{
				// Show left/right separation of frequencies.
				int frames = std::min<int>( visualizer_fft_frames, visualizer_frames - visualizer_start_sample / userdata.Spec.channels );
				float base = pow( 2., log2( visualizer_fft_frames*0.4375 ) / (drawto->h + visualizer_fft_height_offset) );
				for( int ch = (userdata.Spec.channels > 1) ? 1 : 0; ch >= 0; ch -- )
				{
					FFTComplex *visualizer_fft_complex = (ch % 2) ? visualizer_fft_complex_r : visualizer_fft_complex_l;
					memset( visualizer_fft_complex, 0, visualizer_fft_frames * sizeof(FFTComplex) );
					for( int i = 0; i < frames; i ++ )
						for( int in_ch = ch; in_ch < userdata.Spec.channels; in_ch += 2 )
							visualizer_fft_complex[ i ].re += userdata.VisualizerSample( visualizer_start_sample + (i * userdata.Spec.channels) + in_ch );
					av_fft_permute( visualizer_fft_context, visualizer_fft_complex );
					av_fft_calc( visualizer_fft_context, visualizer_fft_complex );
				}
				const FFTComplex *spectrum_l = visualizer_fft_complex_l;
				const FFTComplex *spectrum_r = (userdata.Spec.channels == 1) ? visualizer_fft_complex_l : visualizer_fft_complex_r;
				int offset = 0;
				for( int y = 0; y < drawto->h; y ++ )
				{
					int band_min = std::min<int>( visualizer_fft_frames - 1, offset + pow( base, y ) );
					int band_max = std::min<int>( visualizer_fft_frames, offset + pow( base, y + 1 ) );
					if( band_min == band_max )
					{
						offset ++;
						band_max ++;
					}
					float harmonics = 1.f + log2( (visualizer_fft_frames/2) / band_max );
					float amplitude_l = 0.f, amplitude_r = 0.f;
					for( int band = band_min; band < band_max; band ++ )
					{
						amplitude_l += sqrt( spectrum_l[ band ].re * spectrum_l[ band ].re + spectrum_l[ band ].im * spectrum_l[ band ].im );
						amplitude_r += sqrt( spectrum_r[ band ].re * spectrum_r[ band ].re + spectrum_r[ band ].im * spectrum_r[ band ].im );
					}
					float amplitude = amplitude_l + amplitude_r;
					if( ! amplitude )
						continue;
					float right_ratio = amplitude_r / amplitude;
					float separation = fabs( 0.5f - right_ratio ) * 2.f;
					amplitude /= 4.f * std::max<float>( 7.f, pow(2.f,harmonics) );
					int w = std::max<int>( 1, std::min<int>( drawto->w - 1, drawto->w * std::min<float>(0.5f,amplitude) ) );
					int c = std::max<int>( 0, std::min<int>( drawto->w - 1, drawto->w * right_ratio ) );
					int x_min = std::max<int>( 0, c - w/2 );
					int x_max = std::min<int>( drawto->w - 1, c + w/2 );
					Uint32 color = (separation * separation * amplitude > 0.015625f) ? colors[ visualizer_color1 ] : colors[ visualizer_color2 ];
					for( int x = x_min; x <= x_max; x ++ )
						pixels[ (drawto->h - y - 1) * drawto->w + x ] = color;
				}
				visualizer_fft_height_offset = offset;
			}
			
			font1.Color = colors[ visualizer_text_color ];
			font2.Color = font1.Color;
			uint32_t shadow = SDL_MapRGB( screen->format, 0x00, 0x00, 0x00 );
			
			if( userdata.Crossfading )
				visualizer_scoot ++;
			else
				visualizer_scoot = 0;
			int title_scoot = std::max<int>( 0, (visualizer_scoot - 50) / 2 );
			int artist_scoot = visualizer_scoot / 3;
			
			font2.Draw( 2 - title_scoot, 3, userdata.Title, shadow );
			font2.Draw( 3 - title_scoot, 3, userdata.Title, shadow );
			font2.Draw( 2 - title_scoot, 2, userdata.Title );
			font1.Draw( 2 - artist_scoot, 3 + font2.CharH, userdata.Artist, shadow );
			font1.Draw( 3 - artist_scoot, 3 + font2.CharH, userdata.Artist, shadow );
			font1.Draw( 2 - artist_scoot, 2 + font2.CharH, userdata.Artist );
			int x_offset = font1.CharW * strlen(userdata.Artist);
			if( x_offset && userdata.Album[ 0 ] )
			{
				font1.Draw( 2 - artist_scoot + x_offset, 3 + font2.CharH, " - ", shadow );
				font1.Draw( 3 - artist_scoot + x_offset, 3 + font2.CharH, " - ", shadow );
				font1.Draw( 2 - artist_scoot + x_offset, 2 + font2.CharH, " - " );
				x_offset += 3 * font1.CharW;
			}
			font1.Draw( 2 - artist_scoot + x_offset, 3 + font2.CharH, userdata.Album, shadow );
			font1.Draw( 3 - artist_scoot + x_offset, 3 + font2.CharH, userdata.Album, shadow );
			font1.Draw( 2 - artist_scoot + x_offset, 2 + font2.CharH, userdata.Album );
			
			if( time(NULL) <= userdata.MessageUntil )
			{
				font2.Draw( 2, drawto->h - font2.CharH + 1, userdata.Message, shadow );
				font2.Draw( 3, drawto->h - font2.CharH + 1, userdata.Message, shadow );
				font2.Draw( 2, drawto->h - font2.CharH, userdata.Message );
			}
		}
		else if( pixels && ! visualizer )
		{
			// No visualizer, so just draw the text black-on-white.
			
			Uint32 bg = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0xFF );
			
			for( int x = 0; x < drawto->w; x ++ )
			{
				for( int y = 0; y < drawto->h; y ++ )
					pixels[ y * drawto->w + x ] = bg;
			}
			
			font1.SetTarget( pixels, drawto->w, drawto->h );
			font2.SetTarget( pixels, drawto->w, drawto->h );			
			font1.Color = SDL_MapRGB( screen->format, 0x00, 0x00, 0x00 );
			font2.Color = font1.Color;
			font2.Draw( 2, 2, userdata.Title );
			font1.Draw( 2, 2 + font2.CharH, userdata.Artist );
			int x_offset = font1.CharW * strlen(userdata.Artist);
			if( x_offset && userdata.Album[ 0 ] )
			{
				font1.Draw( 2 + x_offset, 2 + font2.CharH, " - " );
				x_offset += 3 * font1.CharW;
			}
			font1.Draw( 2 + x_offset, 2 + font2.CharH, userdata.Album );
			
			if( time(NULL) <= userdata.MessageUntil )
				font2.Draw( 2, drawto->h - font2.CharH, userdata.Message );
		}
		
		// Finish drawing.
		if( pixels )
		{
			if( screen != drawto )
			{
				int scale = screen->w / drawto->w;
				SDL_LockSurface( screen );
				Uint32 *dest = (Uint32*) screen->pixels;
				
				for( int x = 0; x < screen->w; x ++ )
				{
					for( int y = 0; y < screen->h; y ++ )
						dest[ y * screen->w + x ] = pixels[ (y/scale) * drawto->w + (x/scale) ];
				}
			}
			
			SDL_UnlockSurface( screen );
			SDL_UpdateRect( screen, 0, 0, screen->w, screen->h );
		}
		
		// If we're not playing back audio and the songs are ready, advance by 1MB.
		if( unlikely( ! playback ) && ( (userdata.Songs.size() >= 2) || (userdata.Loader.Finished && userdata.Queue.empty()) ) )
		{
			Uint8 buffer[ 1024*1024 ];
			want.callback( (void*) &userdata, buffer, 1024*1024 );
		}
		
		userdata.Running = userdata.Running && (userdata.Queue.size() || userdata.Songs.size() || userdata.Buffer.Buffered);
	}
	
	if( userdata.WriteTo )
	{
		// Fix the size we wrote to the WAV.
		fflush( userdata.WriteTo );
		fseek( userdata.WriteTo, 4, SEEK_SET );
		size_t bytes = userdata.WroteBytes + 36;
		putc( bytes & 0xFF, userdata.WriteTo );
		putc( bytes/256 & 0xFF, userdata.WriteTo );
		putc( bytes/(256*256) & 0xFF, userdata.WriteTo );
		putc( bytes/(256*256*256) & 0xFF, userdata.WriteTo );
		fseek( userdata.WriteTo, 40, SEEK_SET );
		bytes -= 36;
		putc( bytes & 0xFF, userdata.WriteTo );
		putc( bytes/256 & 0xFF, userdata.WriteTo );
		putc( bytes/(256*256) & 0xFF, userdata.WriteTo );
		putc( bytes/(256*256*256) & 0xFF, userdata.WriteTo );
		fclose( userdata.WriteTo );
	}
	
	// Cleanup before quitting.
	userdata.Loader.Running = userdata.Running;
	av_fft_end( visualizer_fft_context );
	av_free( visualizer_fft_complex_l );
	av_free( visualizer_fft_complex_r );
	if( ! userdata.Loader.Finished )
	{
		printf( "Waiting for loader thread...\n" );
		size_t wait_count = 0;
		while( ! userdata.Loader.Finished )
		{
			if( wait_count >= 10 )
			{
				printf( "Loader thread took too long!  Giving up on it...\n" );
				break;
			}
			SDL_Delay( 1000 );
			wait_count ++;
		}
	}
	userdata.Queue.clear();
	while( userdata.Songs.size() )
	{
		Song *front = userdata.Songs.front();
		userdata.Songs.pop_front();
		delete front;
	}
	
	// Quit SDL.
	if( sdl_audio )
	{
		SDL_LockAudio();
		SDL_CloseAudio();
		SDL_UnlockAudio();
	}
	SDL_Quit();
	
	printf( "Done quitting.\n" );
	
	return 0;
}
