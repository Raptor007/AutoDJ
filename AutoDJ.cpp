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
#ifdef WAVEOUT
#include <windows.h>
#endif
extern "C" {
#include <libavcodec/avfft.h>
}

class Song;
class SongLoader;
class UserData;

namespace ResampleMethod
{
	enum
	{
		Auto = 0,
		Nearest,
		Cubic
	};
}

int SongLoaderThread( void *data_ptr );
double LinearCrossfade( double a, double b, double b_percent );
double EqualPowerCrossfade( double a, double b, double b_percent );
Sint16 Finalize( double value );
bool CalculateCrossfade( const Song *current_song, const Song *next_song, double *crossfade_for_beats, double *crossfade_at_beat );
void AudioCallback( void *userdata, Uint8* stream, int len );
void BufferedAudioCallback( void *userdata, Uint8* stream, int len );
void UnbufferedAudioCallback( void *userdata, Uint8* stream, int len );
void UpdateVisualizerBuffer( UserData *ud, Uint8 *stream, int len );
std::deque<std::string> DirSongs( std::string path );

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


class Song
{
public:
	AudioFile Audio;
	size_t TotalFrames;
	double BPM;
	double CurrentFrame, FirstBeatFrame, FirstOutroFrame;
	int IntroBeats, OutroBeats;
	std::map<size_t,double> Beats;
	double VolumeMax, VolumeAverage;
	
	Song( std::string filename )
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
		
		if( Audio.Load( filename.c_str() ) )
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
	
	double NearestFrame( Uint8 channel ) const
	{
		channel %= Audio.Channels;
		Sint16 *buffer16 = (Sint16*) Audio.Data;  // FIXME: Use BytesPerSample.
		size_t a_index = Audio.Channels * (size_t) CurrentFrame + channel;
		size_t b_index = a_index + Audio.Channels;
		Sint16 a = likely(a_index < Audio.Size / Audio.BytesPerSample) ? buffer16[ a_index ] : 0;
		Sint16 b = likely(b_index < Audio.Size / Audio.BytesPerSample) ? buffer16[ b_index ] : 0;
		double unused = 0.;
		double b_part = modf( CurrentFrame, &unused );
		return (b_part >= 0.5) ? b : a;
	}
	
	double CubicFrame( Uint8 channel ) const
	{
		channel %= Audio.Channels;
		Sint16 *buffer16 = (Sint16*) Audio.Data;  // FIXME: Use BytesPerSample.
		size_t a_index = Audio.Channels * (size_t) CurrentFrame + channel;
		size_t b_index = a_index + Audio.Channels;
		long prev_index = a_index - Audio.Channels;
		size_t next_index = b_index + Audio.Channels;
		double a = 0., b = 0., prev = 0., next = 0.;
		if(likely( (prev_index >= 0) && (next_index < Audio.Size / Audio.BytesPerSample) ))
		{
			prev = buffer16[ prev_index ];
			a = buffer16[ a_index ];
			b = buffer16[ b_index ];
			next = buffer16[ next_index ];
		}
		else
		{
			prev = ((prev_index >= 0) && (prev_index < (long)( Audio.Size / Audio.BytesPerSample ))) ? buffer16[ prev_index ] : 0.;
			a = (a_index < Audio.Size / Audio.BytesPerSample) ? buffer16[ a_index ] : 0.;
			b = (b_index < Audio.Size / Audio.BytesPerSample) ? buffer16[ b_index ] : 0.;
			next = (next_index < Audio.Size / Audio.BytesPerSample) ? buffer16[ next_index ] : 0.;
		}
		double unused = 0.;
		double b_part = modf( CurrentFrame, &unused );
		/*
		double along_a_tangent = a + b_part * (b - prev) / 2.;
		double along_b_tangent = b - (1. - b_part) * (next - a) / 2.;
		double remaining_portion = 1. - b_part * b_part - (1. - b_part) * (1. - b_part);
		double linear = a + (b - a) * b_part;
		return along_a_tangent * (1. - b_part) * (1. - b_part) + along_b_tangent * b_part * b_part + linear * remaining_portion;
		*/
		// Paul Breeuwsma came up with a simpler cubic interpolation than mine, so I'm using it.
		// http://www.paulinternet.nl/?page=bicubic
		return a + 0.5 * b_part * (b - prev + b_part * (2. * prev - 5. * a + 4. * b - next + b_part * (3. * (a - b) + next - prev)));
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
	
	double VolumeAdjustment( void ) const
	{
		return 0.25 / VolumeAverage;
	}
	
	double VolumeAdjustmentToMax( void ) const
	{
		return 1. / VolumeMax;
	}
	
	void Analyze( size_t first_frame = 0 )
	{
		if( !( Audio.Data && Audio.Size ) )
			return;
		
		Sint16 *buffer16 = (Sint16*) Audio.Data;  // FIXME: Use BytesPerSample.
		
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
				double this_amp = fabs( buffer16[ frame * Audio.Channels + channel ] ) / 32767.;
				if( this_amp > VolumeMax )
					VolumeMax = this_amp;
			}
		}
		VolumeAverage = 0.;
		size_t samples_in_average = 0;
		for( size_t frame = first_frame; frame < max_analysis_frame; frame ++ )
		{
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
			{
				double this_amp = fabs( buffer16[ frame * Audio.Channels + channel ] ) / 32767.;
				if( this_amp > VolumeMax / 8. )
				{
					VolumeAverage += this_amp;
					samples_in_average ++;
				}
			}
		}
		VolumeAverage /= samples_in_average;
		
		SDL_Delay( 1 );
		
		
		// Apply low-pass filter and keep track of average peak height over a period of samples.
		
		size_t lpf_samples = std::min<size_t>( max_analysis_frame - first_frame, Audio.SampleRate * 2 / 49 ); // 44100 -> 1800
		double prev = 0., prev_abs = 0., prev_prev_abs = 0.;
		std::map<size_t,double> avg_peak;
		int num_nearby_peaks = 1;
		double highest = 0.;
		
		// Get starting block sum for the low-pass filter.
		for( size_t frame = first_frame; frame < lpf_samples + first_frame; frame ++ )
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
				prev += buffer16[ frame * Audio.Channels + channel ] / (double) Audio.Channels;
		
		SDL_Delay( 1 );
		
		// Process the rest of the audio.
		for( size_t frame = lpf_samples + first_frame; frame < max_analysis_frame; frame ++ )
		{
			// Subtract the oldest sample of the block and add the new one.
			double point = prev;
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
			{
				point -= buffer16[ (frame - lpf_samples) * Audio.Channels + channel ] / (double) Audio.Channels;
				point += buffer16[ frame * Audio.Channels + channel ] / (double) Audio.Channels;
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
		}
		
		// If it's close, round to the nearest whole BPM.
		double bpm_rounding = 0.106;
		double bpm_fpart = modf( Audio.SampleRate * 60. / (double) best_frame_skip, &BPM );
		if( bpm_fpart >= (1. - bpm_rounding) )
			BPM += 1.;
		else if( bpm_fpart > bpm_rounding )
			BPM += bpm_fpart;
		
		FirstBeatFrame = best_first_beat;
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
		Filename = filename;
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
	}
};


// --------------------------------------------------------------------------------------


class UserData
{
public:
	bool Playing;
	bool Running;
	SDL_AudioSpec Spec;
	double BPM;
	bool SourceBPM;
	double SourcePitchScale;
	uint8_t Resample;
	double Volume;
	bool VolumeMatching;
	bool PreventClipping;
	bool Repeat;
	bool Metronome;
	PlaybackBuffer Buffer;
	FILE *WriteTo;
	size_t WroteBytes;
	std::deque<std::string> Queue;
	std::deque<Song*> Songs;
	SongLoader Loader;
	Uint8 *VisualizerBuffer;
	int VisualizerBufferSize;
	clock_t VisualizerClock;
	
	UserData( void )
	{
		Playing = false;
		Running = true;
		memset( &Spec, 0, sizeof(Spec) );
		BPM = 140.;
		SourceBPM = true;
		SourcePitchScale = 1.;
		Resample = ResampleMethod::Auto;
		Volume = 1.;
		VolumeMatching = true;
		PreventClipping = false;
		Repeat = true;
		Metronome = false;
		Buffer.SetSize( 131072 );
		WriteTo = NULL;
		WroteBytes = 0;
		VisualizerBuffer = NULL;
		VisualizerBufferSize = 0;
		VisualizerClock = 0;
	}
	
	void QueueSong( const char *filename )
	{
		Queue.push_back( filename );
	}
	
	void CheckSongLoading( void )
	{
		if( Running )
		{
			if( Loader.Finished )
			{
				// Get the song pointer and take it away from the loader.
				Song *song = Loader.LoadingSong;
				Loader.LoadingSong = NULL;
				
				// Make sure the song loaded okay.
				if( song )
				{
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
						Playing = true;
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
		}
		else
			Loader.Running = false;
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
			loader->LoadingSong = new Song( loader->Filename );
			
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
			
			// Notify main thread that it should now clean up memory.
			loader->Finished = true;
		}
		
		// The loader can wait 5sec each time.
		SDL_Delay( 5000 );
	}
	
	return 0;
}


// --------------------------------------------------------------------------------------


double LinearCrossfade( double a, double b, double b_percent )
{
	return a * (1. - b_percent) + b * b_percent;
}


double EqualPowerCrossfade( double a, double b, double b_percent )
{
	return a * cos( b_percent * M_PI * 0.5 ) + b * cos( (1. - b_percent) * M_PI * 0.5 );
}


Sint16 Finalize( double value )
{
	if( value >= 0. )
	{
		if( value > 32767. )
			return 32767;
		else
			return value + 0.5;
	}
	else
	{
		if( value < -32768. )
			return -32768;
		else
			return value - 0.5;
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
	else
	{
		// If FirstOutroFrame is not specified, guess a good crossfade point.
		double expected_beat = current_total_beats - fmod( current_total_beats, 16. ) - current_song->OutroBeats;
		double actual_beat = current_song->NearestBeatAtBeat( expected_beat );
		*crossfade_at_beat = (fabs( actual_beat - expected_beat ) < 32.) ? actual_beat : expected_beat;
	}
	
	// Make sure the current song doesn't end before the cross-fade is complete.
	if( *crossfade_for_beats > current_total_beats - *crossfade_at_beat )
		*crossfade_for_beats = current_total_beats - *crossfade_at_beat;
	
	return true;
}


// --------------------------------------------------------------------------------------


void AudioCallback( void *userdata, Uint8* stream, int len )
{
	memset( stream, 0, len );
	
	UserData *ud = (UserData*) userdata;
	
	// If we're not filling an intermediate buffer and playback is paused, leave it silent.
	if( !( ud->Playing || ud->Buffer.BufferSize ) )
		return;
	
	// For simplicity's sake, assume it's always 16-bit audio output.
	Sint16 *stream16 = (Sint16*) stream;
	
	bool calculated_crossfade = false;
	double crossfade_for_beats = 96.;
	double crossfade_at_beat = 0.;
	bool crossfade_now = false;
	double crossfade = 0.;
	
	if( ud->Songs.size() )
	{
		Song *current_song = ud->Songs.front();
		Song *next_song = (ud->Songs.size() >= 2) ? ud->Songs[ 1 ] : NULL;
		
		double bpm = ud->SourceBPM ? (current_song->BPM * ud->SourcePitchScale) : ud->BPM;
		
		double volume = ud->VolumeMatching ? ud->Volume * current_song->VolumeAdjustment() : ud->Volume;
		double volume_next = (ud->VolumeMatching && next_song) ? ud->Volume * next_song->VolumeAdjustment() : ud->Volume;
		if( ud->PreventClipping )
		{
			volume = std::min<double>( current_song->VolumeAdjustmentToMax(), volume );
			if( next_song )
				volume_next = std::min<double>( next_song->VolumeAdjustmentToMax(), volume_next );
		}
		
		for( int i = 0; i < len / 2; i += ud->Spec.channels )
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
					crossfade = (current_song->Beat() - crossfade_at_beat) / crossfade_for_beats;
					if( crossfade > 1. )
						crossfade = 1.;
				}
			}
			
			
			if( ! crossfade_now )
			{
				// Not crossfading yet.
				
				if( (ud->Resample == ResampleMethod::Nearest) || ((ud->Resample == ResampleMethod::Auto) && (fabs(current_song->BPM - bpm) < 0.05)) )
				{
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( current_song->NearestFrame( channel ) * volume );
				}
				else
				{
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( current_song->CubicFrame( channel ) * volume );
				}
				
				current_song->Advance( bpm, ud->Spec.freq );
			}
			else
			{
				// Crossfading now.
				
				if( ud->SourceBPM )
					bpm = LinearCrossfade( current_song->BPM * ud->SourcePitchScale, next_song->BPM * ud->SourcePitchScale, crossfade );
				
				bool a_nearest = ( (ud->Resample == ResampleMethod::Nearest) || ((ud->Resample == ResampleMethod::Auto) && (fabs(current_song->BPM - bpm) < 0.05)) );
				bool b_nearest = ( (ud->Resample == ResampleMethod::Nearest) || ((ud->Resample == ResampleMethod::Auto) && (fabs(next_song->BPM - bpm) < 0.05)) );
				
				for( int channel = 0; channel < ud->Spec.channels; channel ++ )
				{
					double a = (a_nearest ? current_song->NearestFrame( channel ) : current_song->CubicFrame( channel )) * volume;
					double b = (b_nearest ? next_song->NearestFrame( channel ) : next_song->CubicFrame( channel )) * volume_next;
					stream16[ i + channel ] = Finalize( EqualPowerCrossfade( a, b, crossfade ) );
				}
				
				current_song->Advance( bpm, ud->Spec.freq );
				next_song->Advance( bpm, ud->Spec.freq );
				
				// If we completed a crossfade, remove the song we faded from.
				if(unlikely( crossfade >= 1. ))
				{
					ud->Songs.pop_front();
					delete current_song;
					current_song = next_song;
					next_song = (ud->Songs.size() >= 2) ? ud->Songs[ 1 ] : NULL;
					
					calculated_crossfade = false;
					crossfade_now = false;
					crossfade = 0.;
				}
			}
			
			
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
					double hz = (((int)( beat_ipart + 0.5 )) % 4) ? 960. : 1920.;
					double wave_value = (sin( beat_fpart * hz ) * 8000. * cos( beat_fpart * M_PI * 2. ));
					
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( stream16[ i + channel ] + wave_value );
				}
			}
		}
		
		// If we're supposed to write our output to a file, do it after each buffer we fill.
		if( ud->WriteTo )
		{
			fwrite( stream, 1, len, ud->WriteTo );
			ud->WroteBytes += len;
		}
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
			// This is because AudioCallback isn't fully thread-safe when songs end!
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
			
			#ifdef WIN32
				#define PATH_SEPARATOR "\\"
			#else
				#define PATH_SEPARATOR "/"
			#endif
			
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
		songs.push_back( path );
	
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


int main( int argc, char **argv )
{
	// If no arguments were given, search the default music directory.
	if( argc == 1 )
	{
		char cmd[ 102400 ] = "";
		#ifdef WIN32
			const char *DEFAULT_ARGS = "\"M:\\iTunes\\iTunes\\ Music\\Trance\"";
		#else
			const char *DEFAULT_ARGS = "~/Music/iTunes/iTunes\\ Music/Trance /Volumes/Media/Music/iTunes/iTunes\\ Music/Trance &";
		#endif
		snprintf( cmd, 102400, "\"%s\" %s", argv[ 0 ], DEFAULT_ARGS );
		system( cmd );
		return 0;
	}
	
	UserData userdata;
	bool shuffle = true;
	bool window = true;
	bool fullscreen = false;
	int visualizer = 1;
	bool playback = true;
	bool sdl_audio = true;
#ifdef WAVEOUT
	// Default to WaveOut.
	sdl_audio = false;
#endif
	const char *write = NULL;
	SDL_AudioSpec want;
	memset( &want, 0, sizeof(want) );
	want.freq = 44100;
	want.format = AUDIO_S16;
	want.channels = 2;
	want.samples = 4096;
	want.callback = BufferedAudioCallback;
	want.userdata = &userdata;
	
	// Process command-line arguments.
	for( int i = 1; i < argc; i ++ )
	{
		if( strncmp( argv[ i ], "--", 2 ) == 0 )
		{
			if( strcasecmp( argv[ i ], "--no-window" ) == 0 )
				window = false;
#ifdef WAVEOUT
			else if( strcasecmp( argv[ i ], "--sdl-audio" ) == 0 )
				sdl_audio = true;
#endif
			else if( strcasecmp( argv[ i ], "--fullscreen" ) == 0 )
			{
				window = true;
				fullscreen = true;
			}
			else if( strncasecmp( argv[ i ], "--visualizer=", strlen("--visualizer=") ) == 0 )
				visualizer = atoi( argv[ i ] + strlen("--visualizer=") );
			else if( strcasecmp( argv[ i ], "--no-shuffle" ) == 0 )
				shuffle = false;
			else if( strcasecmp( argv[ i ], "--no-repeat" ) == 0 )
				userdata.Repeat = false;
			else if( strcasecmp( argv[ i ], "--metronome" ) == 0 )
				userdata.Metronome = true;
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
			else if( strncasecmp( argv[ i ], "--volume=", strlen("--volume=") ) == 0 )
				userdata.Volume = atof( argv[ i ] + strlen("--volume=") );
			else if( strcasecmp( argv[ i ], "--no-volume-matching" ) == 0 )
				userdata.VolumeMatching = false;
			else if( strcasecmp( argv[ i ], "--prevent-clipping" ) == 0 )
				userdata.PreventClipping = true;
			else if( strncasecmp( argv[ i ], "--rate=", strlen("--rate=") ) == 0 )
				want.freq = atoi( argv[ i ] + strlen("--rate=") );
			else if( strncasecmp( argv[ i ], "--channels=", strlen("--channels=") ) == 0 )
				want.channels = atoi( argv[ i ] + strlen("--channels=") );
			else if( strncasecmp( argv[ i ], "--buffer1=", strlen("--buffer1=") ) == 0 )
				want.samples = atoi( argv[ i ] + strlen("--buffer1=") );
			else if( strncasecmp( argv[ i ], "--buffer2=", strlen("--buffer2=") ) == 0 )
				userdata.Buffer.SetSize( atoi( argv[ i ] + strlen("--buffer2=") ) * 2 * want.channels );
			else if( strcasecmp( argv[ i ], "--no-playback" ) == 0 )
				playback = false;
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
			
			// Get all playable files from within the selected directory, or add selected songs directly.
			std::deque<std::string> songs = DirSongs( argv[ i ] );
			for( std::deque<std::string>::const_iterator song_iter = songs.begin(); song_iter != songs.end(); song_iter ++ )
				userdata.QueueSong( (*song_iter).c_str() );
		}
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
		screen = fullscreen ? SDL_SetVideoMode( 0, 0, 0, SDL_SWSURFACE | SDL_FULLSCREEN ) : SDL_SetVideoMode( 256, 64, 0, SDL_SWSURFACE );
	}
	
	// Prepare audio file input.
	av_log_set_level( AV_LOG_FATAL );
	av_register_all();
	
	// Shuffle song list if requested.
	if( shuffle )
		std::random_shuffle( userdata.Queue.begin(), userdata.Queue.end() );
	
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
				wfx.wFormatTag = WAVE_FORMAT_PCM;
				wfx.wBitsPerSample = 16;
				wfx.nChannels = want.channels;
				wfx.nSamplesPerSec = want.freq;
				wfx.nBlockAlign = wfx.nChannels * wfx.wBitsPerSample / 8;
				wfx.nAvgBytesPerSec = wfx.nBlockAlign * wfx.nSamplesPerSec;
				wfx.cbSize = 0;
				MMRESULT wave_out_result = waveOutOpen( &WaveOutHandle, WAVE_MAPPER, &wfx, (DWORD_PTR) &WaveOutCallback, 0, CALLBACK_FUNCTION );
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
			SDL_OpenAudio( &want, &(userdata.Spec) );
			
			if( screen )
			{
				userdata.VisualizerBufferSize = userdata.Spec.samples * userdata.Spec.channels * 2;
				userdata.VisualizerBuffer = (Uint8*) malloc( userdata.VisualizerBufferSize );
			}
		}
	}
	
	if( write )
	{
		// Write WAV header.
		unsigned char wave_header[ 44 ] = { 'R','I','F','F', 36,0xFF,0xFF,0x7F, 'W','A','V','E', 'f','m','t',0x20, 16,0,0,0, 1,0, userdata.Spec.channels,0, userdata.Spec.freq%256,userdata.Spec.freq/256,0,0, (userdata.Spec.freq*userdata.Spec.channels*2)%256,(userdata.Spec.freq*userdata.Spec.channels*2)/256,0,0, 4,0, 16,0, 'd','a','t','a', 0,0xFF,0xFF,0x7F };
		userdata.WriteTo = fopen( write, "wb" );
		fwrite( wave_header, 1, 44, userdata.WriteTo );
		fflush( userdata.WriteTo );
	}
	
	// Visualizer variables.
	clock_t visualizer_updated_clock = 0, visualizer_prev_clock = 0;
	size_t visualizer_start_sample = 0;
	size_t visualizer_loading_frame = 0;
	int visualizer_color1 = 1, visualizer_color2 = 0;
	int visualizer_backgrounds[ 4 ] = { 0, 2, 0, 1 };
	int visualizer_frames = userdata.VisualizerBufferSize / (userdata.Spec.channels * 2);
	int visualizer_fft_frames = visualizer_frames / 4;
	FFTContext *visualizer_fft_context = av_fft_init( log2(visualizer_fft_frames), false );
	FFTComplex *visualizer_fft_complex_l = (FFTComplex*) av_mallocz( visualizer_fft_frames * sizeof(FFTComplex) );
	FFTComplex *visualizer_fft_complex_r = (FFTComplex*) av_mallocz( visualizer_fft_frames * sizeof(FFTComplex) );
	int visualizer_fft_width_offset = 0;
	int visualizer_fft_height_offset = 0;
	
	// Keep running until playback is complete.
	userdata.Running = (userdata.Queue.size() || userdata.Songs.size());
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
					else if( key == SDLK_MINUS )
					{
						userdata.Volume -= 0.0625;
						printf( "Volume: %.3f\n", userdata.Volume );
						fflush( stdout );
					}
					else if( key == SDLK_EQUALS )
					{
						userdata.Volume += 0.0625;
						printf( "Volume: %.3f\n", userdata.Volume );
						fflush( stdout );
					}
					else if( key == SDLK_BACKSPACE )
					{
						userdata.VolumeMatching = ! userdata.VolumeMatching;
						printf( "Volume Matching: %s\n", userdata.VolumeMatching ? "On" : "Off" );
						fflush( stdout );
					}
					else if( key == SDLK_p )
					{
						userdata.PreventClipping = ! userdata.PreventClipping;
						printf( "Prevent Clipping: %s\n", userdata.PreventClipping ? "On" : "Off" );
						fflush( stdout );
					}
					else if( key == SDLK_LEFTBRACKET || key == SDLK_LEFT )
					{
						Song *current_song = userdata.Songs.front();
						current_song->CurrentFrame = 0.;
						current_song->FirstOutroFrame = 0.;
						current_song->OutroBeats = 128;
					}
					else if( key == SDLK_RIGHTBRACKET || key == SDLK_RIGHT )
					{
						Song *current_song = userdata.Songs.front();
						int beat = current_song->Beat();
						if( current_song->FirstOutroFrame && (current_song->FrameAtBeat(beat) >= current_song->FirstOutroFrame) )
							current_song->OutroBeats /= 2;
						else
						{
							double expected_beat = beat - (beat % 16);
							double actual_beat = current_song->NearestBeatAtBeat( expected_beat );
							current_song->FirstOutroFrame = current_song->FrameAtBeat( (fabs( actual_beat - expected_beat ) < 32.) ? actual_beat : expected_beat );
							current_song->OutroBeats = std::min<int>( 64, current_song->TotalBeats() - beat );
						}
					}
					else if( key == SDLK_BACKSLASH )
						userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->TotalFrames;
					else if( key == SDLK_k )
					{
						if( userdata.SourceBPM && userdata.Songs.size() )
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
						
						userdata.BPM -= 1.;
						userdata.SourceBPM = false;
						printf( "Playback BPM: %.2f\n", userdata.BPM );
						fflush( stdout );
					}
					else if( key == SDLK_l )
					{
						if( userdata.SourceBPM && userdata.Songs.size() )
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
						
						userdata.BPM += 1.;
						userdata.SourceBPM = false;
						printf( "Playback BPM: %.2f\n", userdata.BPM );
						fflush( stdout );
					}
					else if( key == SDLK_SEMICOLON || key == SDLK_DOWN )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale /= pow( 2., 1./12. );
						printf( "Pitch: %.1f%%\n", userdata.SourcePitchScale * 100. );
						fflush( stdout );
					}
					else if( key == SDLK_QUOTE || key == SDLK_UP )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale *= pow( 2., 1./12. );
						printf( "Pitch: %.1f%%\n", userdata.SourcePitchScale * 100. );
						fflush( stdout );
					}
					else if( key == SDLK_m )
						userdata.Metronome = ! userdata.Metronome;
					else if( key == SDLK_COMMA )
					{
						int beat = userdata.Songs.front()->Beat();
						userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->FrameAtBeat( beat - (beat % 16) - 32 );
					}
					else if( key == SDLK_PERIOD )
					{
						int beat = userdata.Songs.front()->Beat();
						userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->FrameAtBeat( beat - (beat % 16) + 32 );
					}
					else if( key == SDLK_SLASH )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale = 1.;
						printf( "Pitch/BPM: Match Source\n" );
						fflush( stdout );
					}
					else if( key == SDLK_v )
					{
						if( screen )
						{
							if( ! shift )
								visualizer ++;
							else
								visualizer += 3;
							visualizer %= 4;
						}
						else
							visualizer = 0;
						
						if( screen && ! visualizer )
						{
							Uint32 bg = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0xFF );
							Uint32* pixels = (Uint32*) screen->pixels;
							SDL_LockSurface( screen );
							for( int x = 0; x < screen->w; x ++ )
							{
								for( int y = 0; y < screen->h; y ++ )
									pixels[ y * screen->w + x ] = bg;
							}
							SDL_UnlockSurface( screen );
							SDL_UpdateRect( screen, 0, 0, screen->w, screen->h );
						}
					}
					else if( key == SDLK_c )
					{
						int visualizer_colors = visualizer_color1 * 3 + visualizer_color2;
						if( ! shift )
							visualizer_colors ++;
						else
							visualizer_colors += 8;
						visualizer_colors %= 9;
						visualizer_color1 = visualizer_colors / 3;
						visualizer_color2 = visualizer_colors % 3;
					}
					else if( key == SDLK_b )
					{
						if( visualizer && userdata.Playing )
						{
							if( ! shift )
								visualizer_backgrounds[ visualizer ] ++;
							else
								visualizer_backgrounds[ visualizer ] += 2;
							visualizer_backgrounds[ visualizer ] %= 3;
						}
					}
				}
			}
		}
		
		if( userdata.Buffer.BufferSize )
		{
			// Add some data to the playback buffer.
			userdata.Buffer.AddToBuffer( &userdata, 16384 );
		}
		
		// Wait 10ms.
		SDL_Delay( 10 );
		
		// Visualize the waveform in the window.
		if( visualizer && screen && (userdata.Playing || userdata.Songs.empty()) )
		{
			// Keep the visualizer synchronized with the audio playback.
			if( userdata.Playing )
			{
				if( userdata.VisualizerClock != visualizer_updated_clock )
				{
					visualizer_updated_clock = userdata.VisualizerClock;
					visualizer_start_sample = 0;
					visualizer_prev_clock = clock();
				}
				else
				{
					clock_t now = clock();
					size_t frames = userdata.Spec.freq * (now - visualizer_prev_clock) / (double) CLOCKS_PER_SEC;
					visualizer_prev_clock = now;
					visualizer_start_sample += userdata.Spec.channels * frames;
					visualizer_start_sample %= userdata.VisualizerBufferSize / (userdata.Spec.channels * 2);
				}
			}
			else
			{
				visualizer_loading_frame ++;
				visualizer_loading_frame %= screen->h;
				visualizer_prev_clock = clock();
			}
			
			Uint32 colors[ 3 ];
			colors[ 0 ] = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0xFF );
			colors[ 1 ] = SDL_MapRGB( screen->format, 0xFF, 0xFF, 0x00 );
			colors[ 2 ] = SDL_MapRGB( screen->format, 0x00, 0xFF, 0x00 );
			Uint32* pixels = (Uint32*) screen->pixels;
			SDL_LockSurface( screen );
			
			if( (visualizer_backgrounds[ visualizer ] == 0) || ! userdata.Playing )
			{
				// Fade in-place.
				for( int x = 0; x < screen->w; x ++ )
				{
					for( int y = 0; y < screen->h; y ++ )
					{
						Uint8 r = 0x00, g = 0x00, b = 0x00;
						SDL_GetRGB( pixels[ y * screen->w + x ], screen->format, &r, &g, &b );
						pixels[ y * screen->w + x ] = SDL_MapRGB( screen->format, std::max<int>( 0, std::min<int>( 0xFF, (r - b) * 0.875f ) ), std::min<int>( 0xFF, g * 0.5f + std::max<float>( 0.f, (g - r - b) * 0.25f ) ), std::min<int>( 0xFF, b * 0.9375f ) );
					}
				}
			}
			else if( visualizer_backgrounds[ visualizer ] == 1 )
			{
				// Fade downward.
				for( int x = 0; x < screen->w; x ++ )
				{
					Uint8 r = 0x00, g = 0x00, b = 0x00;
					for( int y = screen->h - 1; y >= 1; y -- )
					{
						Uint8 r1 = 0x00, g1 = 0x00, b1 = 0x00, r2 = 0x00, g2 = 0x00, b2 = 0x00;
						SDL_GetRGB( pixels[  y    * screen->w + x ], screen->format, &r1, &g1, &b1 );
						SDL_GetRGB( pixels[ (y-1) * screen->w + x ], screen->format, &r2, &g2, &b2 );
						r = std::max<Uint8>( r1, r2 );
						g = std::max<Uint8>( g1, g2 );
						b = std::max<Uint8>( b1, b2 );
						pixels[ y * screen->w + x ] = SDL_MapRGB( screen->format, std::max<int>( 0, std::min<int>( 0xFF, (r - b) * 0.875f ) ), std::min<int>( 0xFF, g * 0.5f + std::max<float>( 0.f, (g - r - b) * 0.25f ) ), std::min<int>( 0xFF, b * 0.9375f ) );
					}
					SDL_GetRGB( pixels[ x ], screen->format, &r, &g, &b );
					pixels[ x ] = SDL_MapRGB( screen->format, std::max<int>( 0, std::min<int>( 0xFF, (r - b) * 0.875f ) ), std::min<int>( 0xFF, g * 0.5f + std::max<float>( 0.f, (g - r - b) * 0.25f ) ), std::min<int>( 0xFF, b * 0.9375f ) );
				}
			}
			else if( visualizer_backgrounds[ visualizer ] == 2 )
			{
				// Fade outward.
				for( int x = 0; x < screen->w / 2; x ++ )
				{
					Uint8 r = 0x00, g = 0x00, b = 0x00;
					for( int y = 0; y < screen->h / 2; y ++ )
					{
						for( int q = 0; q < 4; q ++ )
						{
							int x1 = x, y1 = y, dx = 1, dy = 1;
							if( q >= 2 )
							{
								x1 = screen->w - x - 1;
								dx = -1;
							}
							if( q % 2 )
							{
								y1 = screen->h - y - 1;
								dy = -1;
							}
							Uint8 r1 = 0x00, g1 = 0x00, b1 = 0x00, r2 = 0x00, g2 = 0x00, b2 = 0x00;
							SDL_GetRGB( pixels[  y1     * screen->w + x1    ], screen->format, &r1, &g1, &b1 );
							SDL_GetRGB( pixels[ (y1+dy) * screen->w + x1+dx ], screen->format, &r2, &g2, &b2 );
							r = std::max<Uint8>( r1, r2 );
							g = std::max<Uint8>( g1, g2 );
							b = std::max<Uint8>( b1, b2 );
							pixels[ y1 * screen->w + x1 ] = SDL_MapRGB( screen->format, std::max<int>( 0, std::min<int>( 0xFF, (r - b) * 0.875f ) ), std::min<int>( 0xFF, g * 0.5f + std::max<float>( 0.f, (g - r - b) * 0.25f ) ), std::min<int>( 0xFF, b * 0.9375f ) );
						}
					}
				}
			}
			
			if( ! userdata.Playing )
			{
				// Loading animation.
				for( int x = 0; x < screen->w; x ++ )
				{
					int y1 = ((screen->h - x) % screen->h + screen->h + visualizer_loading_frame) % screen->h;
					int y2 = (y1 + screen->h / 2) % screen->h;
					pixels[ y1 * screen->w + x ] = colors[ visualizer_color1 ];
					pixels[ y2 * screen->w + x ] = colors[ visualizer_color2 ];
				}
			}
			
			else if( visualizer == 1 )
			{
				// Fill in the waveform.
				for( int x = 0; x < screen->w; x ++ )
				{
					for( int ch = userdata.Spec.channels - 1; ch >= 0; ch -- )
					{
						Sint16 raw = ((Sint16*)( userdata.VisualizerBuffer ))[ (visualizer_start_sample + (x * userdata.Spec.channels) + ch) % (userdata.VisualizerBufferSize / 2) ];
						float amplitude = raw / 32767.f;
						float volume = userdata.Volume;
						if( (volume > 0.f) && (volume < 1.f) )
							amplitude /= volume;
						int y = std::max<int>( 0, std::min<int>( screen->h - 1, screen->h * (0.5f - amplitude * 0.5f) + 0.5f ) );
						pixels[ y * screen->w + x ] = (ch % 2) ? colors[ visualizer_color2 ] : colors[ visualizer_color1 ];
					}
				}
			}
			
			else if( visualizer == 2 )
			{
				// Show spectrum of frequencies.
				int frames = std::min<int>( visualizer_fft_frames, visualizer_frames - visualizer_start_sample / userdata.Spec.channels );
				float base = pow( 2., log2( visualizer_fft_frames/2 ) / (screen->w + visualizer_fft_width_offset) );
				for( int ch = userdata.Spec.channels - 1; ch >= 0; ch -- )
				{
					FFTComplex *visualizer_fft_complex = (ch % 2) ? visualizer_fft_complex_r : visualizer_fft_complex_l;
					memset( visualizer_fft_complex, 0, visualizer_fft_frames * sizeof(FFTComplex) );
					float scale = 1.f / 32768.f;
					float volume = userdata.Volume;
					if( (volume > 0.f) && (volume < 1.f) )
						scale /= volume;
					for( int i = 0; i < frames; i ++ )
						visualizer_fft_complex[ i ].re = ((Sint16*)( userdata.VisualizerBuffer ))[ visualizer_start_sample + (i * userdata.Spec.channels) + ch ] * scale;
					av_fft_permute( visualizer_fft_context, visualizer_fft_complex );
					av_fft_calc( visualizer_fft_context, visualizer_fft_complex );
					int offset = 0;
					for( int x = 0; x < screen->w; x ++ )
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
						amplitude /= std::max<float>( 8.f, pow(2.f,harmonics) );
						int h = std::max<int>( 0, std::min<int>( screen->h - 1, screen->h * amplitude ) );
						for( int y = screen->h - 1 - h; y < screen->h; y ++ )
							pixels[ y * screen->w + x ] = (ch % 2) ? colors[ visualizer_color1 ] : colors[ visualizer_color2 ];
					}
					visualizer_fft_width_offset = offset;
				}
			}
			
			else if( visualizer == 3 )
			{
				// Show left/right separation of frequencies.
				int frames = std::min<int>( visualizer_fft_frames, visualizer_frames - visualizer_start_sample / userdata.Spec.channels );
				float base = pow( 2., log2( visualizer_fft_frames/8 ) / (screen->h + visualizer_fft_height_offset) );
				for( int ch = (userdata.Spec.channels > 1) ? 1 : 0; ch >= 0; ch -- )
				{
					FFTComplex *visualizer_fft_complex = (ch % 2) ? visualizer_fft_complex_r : visualizer_fft_complex_l;
					memset( visualizer_fft_complex, 0, visualizer_fft_frames * sizeof(FFTComplex) );
					float scale = 1.f / 32768.f;
					float volume = userdata.Volume;
					if( (volume > 0.f) && (volume < 1.f) )
						scale /= volume;
					for( int i = 0; i < frames; i ++ )
						for( int in_ch = ch; in_ch < userdata.Spec.channels; in_ch += 2 )
							visualizer_fft_complex[ i ].re += ((Sint16*)( userdata.VisualizerBuffer ))[ visualizer_start_sample + (i * userdata.Spec.channels) + in_ch ] * scale;
					av_fft_permute( visualizer_fft_context, visualizer_fft_complex );
					av_fft_calc( visualizer_fft_context, visualizer_fft_complex );
				}
				const FFTComplex *spectrum_l = visualizer_fft_complex_l;
				const FFTComplex *spectrum_r = (userdata.Spec.channels == 1) ? visualizer_fft_complex_l : visualizer_fft_complex_r;
				int offset = 0;
				for( int y = 0; y < screen->h; y ++ )
				{
					float amplitude_l = 0.f, amplitude_r = 0.f;
					int band_min = std::min<int>( visualizer_fft_frames - 1, offset + pow( base, y ) );
					int band_max = std::min<int>( visualizer_fft_frames, offset + pow( base, y + 1 ) );
					if( band_min == band_max )
					{
						offset ++;
						band_max ++;
					}
					float harmonics = 1.f + log2( (visualizer_fft_frames/8) / band_max );
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
					amplitude /= std::max<float>( 128.f, pow(2.f,harmonics+5.f) );
					int w = std::max<int>( 1, std::min<int>( screen->w - 1, screen->w * amplitude ) );
					int c = std::max<int>( 0, std::min<int>( screen->w - 1, screen->w * right_ratio ) );
					int x_min = std::max<int>( 0, c - w/2 );
					int x_max = std::min<int>( screen->w - 1, c + w/2 );
					Uint32 color = (separation * amplitude > 0.01f) ? colors[ visualizer_color2 ] : colors[ visualizer_color1 ];
					for( int x = x_min; x <= x_max; x ++ )
						pixels[ (screen->h - y - 1) * screen->w + x ] = color;
				}
				visualizer_fft_height_offset = offset;
			}
			
			SDL_UnlockSurface( screen );
			SDL_UpdateRect( screen, 0, 0, screen->w, screen->h );
		}
		
		// If we're not playing back audio, advance by 1024 bytes.
		// This is mostly for valgrind debugging to avoid alsa buffer underrun warnings.
		if(unlikely( ! playback ))
		{
			Uint8 buffer[ 1024 ];
			want.callback( (void*) &userdata, buffer, 1024 );
		}
#ifdef WAVEOUT
		else if( ! sdl_audio )
			WaveOutCheck( &userdata );
#endif
		
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
		while( ! userdata.Loader.Finished )
			SDL_Delay( 100 );
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
