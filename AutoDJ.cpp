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
#include <algorithm>
#include <string>
#include <SDL/SDL.h>
#include <SDL/SDL_audio.h>

#ifdef USE_SDL_MIXER
#include <SDL/SDL_mixer.h>
#endif

#include "AudioFile.h"
#include "PlaybackBuffer.h"

class Song;
class SongLoadData;
class UserData;

int SongLoad( void *data_ptr );
double LinearCrossfade( double a, double b, double b_percent );
double EqualPowerCrossfade( double a, double b, double b_percent );
Sint16 Finalize( double value );
bool CalculateCrossfade( Song *current_song, const Song *next_song, double *crossfade_for_beats, double *crossfade_at_beat );
void AudioCallback( void *userdata, Uint8* stream, int len );
void BufferedAudioCallback( void *userdata, Uint8* stream, int len );
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


namespace ResampleMethod
{
	enum
	{
		Auto = 0,
		Nearest,
		Cubic
	};
}


// --------------------------------------------------------------------------------------


class Song
{
public:
	AudioFile Audio;
	double BPM;
	double CurrentFrame, MaxFrame, FirstBeatFrame, FirstOutroFrame;
	int IntroBeats, OutroBeats;
	std::map<size_t,double> Beats;
	
	Song( std::string filename, const SDL_AudioSpec *spec )
	{
		BPM = 140.;
		CurrentFrame = 0.;
		MaxFrame = 0.;
		FirstBeatFrame = 0.;
		FirstOutroFrame = 0.;
		IntroBeats = 0;
		OutroBeats = 0;
		
		bool loaded = false;
#ifdef USE_SDL_MIXER
		if( ! loaded )
			loaded = Audio.LoadWithSDLMixer( filename.c_str(), spec, false );
#endif
#ifdef USE_LIBAV
		if( ! loaded )
			loaded = Audio.LoadWithLibAV( filename.c_str() );
#endif
#if defined(USE_SDL_MIXER) && defined(USE_EXTERNAL_FFMPEG)
		if( ! loaded )
			loaded = Audio.LoadWithSDLMixer( filename.c_str(), spec, true );
#endif
		if( loaded )
		{
			MaxFrame = Audio.Size / (Audio.BytesPerSample * Audio.Channels);
			
			if( MaxFrame < 0.5 )
			{
				printf( "%s: zero-length\n", filename.c_str() );
				Audio.Clear();
			}
			else if( MaxFrame / (double) Audio.SampleRate < 30. )
			{
				printf( "%s: only %.1f sec long\n", filename.c_str(), MaxFrame / (double) Audio.SampleRate );
				Audio.Clear();
			}
		}
		else
		{
			MaxFrame = 0.;
			printf( "%s: %s\n", filename.c_str(), Audio.Error );
		}
	}
	
	~Song()
	{
	}
	
	void SetFirstBeat( int minutes, double seconds )
	{
		FirstBeatFrame = Audio.SampleRate * ((60. * minutes) + seconds);
	}
	
	double NearestFrame( Uint8 channel ) const
	{
		channel %= Audio.Channels;
		Sint16 *buffer16 = (Sint16*) Audio.Data;  // FIXME: Use BytesPerSample.
		size_t a_index = Audio.Channels * (size_t) CurrentFrame;
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
		size_t a_index = Audio.Channels * (size_t) CurrentFrame;
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
	
	void Advance( double playback_bpm )
	{
		CurrentFrame += playback_bpm / BPM;
	}
	
	double Beat( void ) const
	{
		return BeatAtFrame( CurrentFrame );
	}
	
	double TotalBeats( void ) const
	{
		return BeatAtFrame( MaxFrame );
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
	
	double NearestBeatFrameAtFrame( double frame )
	{
		std::map<size_t,double>::iterator exact_beat_found = Beats.find( (size_t)( frame + 0.5 ) );
		if( exact_beat_found != Beats.end() )
		{
			// We found an exact match.
			
			return exact_beat_found->first;
		}
		else
		{
			// No exact match, check for beats near where we guessed.
			
			Beats[ frame ] = 0.;
			std::map<size_t,double>::iterator temp = Beats.find( frame );
			std::map<size_t,double>::iterator ahead = temp;
			ahead ++;
			std::map<size_t,double>::reverse_iterator behind(temp);
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
			
			Beats.erase( temp );
			
			return nearest_beat;
		}
	}
	
	double NearestBeatFrameAtBeat( double beat )
	{
		return NearestBeatFrameAtFrame( FrameAtBeat( beat ) );
	}
	
	double NearestBeatAtFrame( double frame )
	{
		return BeatAtFrame( NearestBeatFrameAtFrame( frame ) );
	}
	
	double NearestBeatAtBeat( double beat )
	{
		return BeatAtFrame( NearestBeatFrameAtBeat( beat ) );
	}
	
	bool Finished( void ) const
	{
		return (CurrentFrame >= MaxFrame);
	}
	
	void Analyze( void )
	{
		if( !( Audio.Data && Audio.Size ) )
			return;
		
		// Apply low-pass filter and keep track of average peak height over a period of samples.
		
		Sint16 *buffer16 = (Sint16*) Audio.Data;  // FIXME: Use BytesPerSample.
		
		// Don't analyze more than 10 minutes of audio per track.
		size_t max_analysis_frame = Audio.SampleRate * 60 * 10;
		if( max_analysis_frame > MaxFrame )
			max_analysis_frame = MaxFrame;
		
		double prev = 0., prev_abs = 0.;
		std::map<size_t,double> avg_peak;
		int num_peaks = 1;
		double highest = 0.;
		
		// Parameters for low-pass filter.
		size_t samples = 1800;
		if( samples > max_analysis_frame )
			samples = max_analysis_frame;
		double prev_sample_power = 1.0;
		
		// Get starting block sum for the low-pass filter.
		for( size_t frame = 0; frame < samples; frame ++ )
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
				prev += buffer16[ frame * Audio.Channels + channel ] * pow( prev_sample_power, samples - frame - 1. ) / (double) Audio.Channels;
		
		// Process the rest of the audio.
		for( size_t frame = samples; frame < max_analysis_frame; frame ++ )
		{
			// Subtract the oldest sample of the block and add the new one.
			double point = prev * prev_sample_power;
			for( size_t channel = 0; channel < Audio.Channels; channel ++ )
			{
				point -= buffer16[ (frame - samples) * Audio.Channels + channel ] * pow( prev_sample_power, samples ) / (double) Audio.Channels;
				point += buffer16[ frame * Audio.Channels + channel ] / (double) Audio.Channels;
			}
			
			if( fabs(point) < prev_abs )
			{
				// Found a peak on the previous frame.
				size_t peak_frame = frame - (samples / 2) - 1;
				
				if( avg_peak.size() )
				{
					if( avg_peak.rbegin()->first + 900 < peak_frame )
					{
						avg_peak.rbegin()->second += prev_abs;
						num_peaks ++;
					}
					else
					{
						avg_peak.rbegin()->second /= num_peaks;
						if( avg_peak.rbegin()->second > highest )
							highest = avg_peak.rbegin()->second;
						
						avg_peak[ peak_frame ] = prev_abs;
						num_peaks = 1;
					}
				}
				else
				{
					avg_peak[ peak_frame ] = prev_abs;
					num_peaks = 1;
				}
				
				prev_abs = 0.;
			}
			else
				prev_abs = fabs(point);
			
			prev = point;
		}
		
		if( avg_peak.size() )
		{
			avg_peak.rbegin()->second /= num_peaks;
			if( avg_peak.rbegin()->second > highest )
				highest = avg_peak.rbegin()->second;
		}
		
		
		// Look for significant increases in peak averages (bass hits).
		
		double min_value = highest / 8., min_inc_factor = 128.;
		prev = 0.;
		Beats.clear();
		bool prev_adjacent = false;
		
		for( std::map<size_t,double>::iterator avg_iter = avg_peak.begin(); avg_iter != avg_peak.end(); avg_iter ++ )
		{
			if( prev_adjacent && (avg_iter->second > prev) )
			{
				std::map<size_t,double>::iterator last = Beats.end();
				last --;
				Beats.erase( last );
				Beats[ avg_iter->first ] = avg_iter->second;
			}
			else if( (avg_iter->second >= min_value) && (avg_iter->second > prev * min_inc_factor) )
			{
				Beats[ avg_iter->first ] = avg_iter->second;
				prev_adjacent = true;
			}
			else
				prev_adjacent = false;
		
			prev = avg_iter->second;
		}
		
		size_t first_beat = avg_peak.begin()->first;
		
		
		// Search a range of likely BPMs and see how well the beats match.
		
		size_t best_frame_skip = Audio.SampleRate * 60 / 140;
		size_t best_first_beat = first_beat;
		double best_error = FLT_MAX;
		int best_doff = INT_MAX;
		
		for( int bpm = 120; (bpm <= 150) && (best_error > 0.); bpm ++ )
		{
			// Search at a BPM and see how closely it matches.
			
			size_t frame_skip = (Audio.SampleRate * 60. / (double) bpm) + 0.5;
			size_t first_beat_matched = 0;
			
			// Try a few different starting points.
			for( std::map<size_t,double>::iterator beat_iter = Beats.begin(); (beat_iter != Beats.end()) && (beat_iter->first < max_analysis_frame / 2); beat_iter ++ )
			{
				size_t start = beat_iter->first;
				std::vector<int> doff_behind, doff_ahead;
				int off_behind = 0, off_ahead = 0;
				int misses_behind = 0, misses_ahead = 0;
				
				// Look at each place we expect the beat to hit.
				for( size_t frame = start; frame < max_analysis_frame; frame += frame_skip )
				{
					std::map<size_t,double>::iterator exact_beat_found = Beats.find( frame );
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
						
						Beats[ frame ] = 0.;
						std::map<size_t,double>::iterator temp = Beats.find( frame );
						std::map<size_t,double>::iterator ahead = temp;
						ahead ++;
						std::map<size_t,double>::reverse_iterator behind(temp);
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
						
						Beats.erase( temp );
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
					for( std::vector<int>::iterator doff_iter = doff_behind.begin(); doff_iter != doff_behind.end(); doff_iter ++ )
						mean_behind += *doff_iter / (double) doff_behind.size();
					
					error_behind = 0.;
					
					for( std::vector<int>::iterator doff_iter = doff_behind.begin(); doff_iter != doff_behind.end(); doff_iter ++ )
					{
						double ddoff = (*doff_iter - mean_behind);
						error_behind += (ddoff * ddoff) / (double) doff_behind.size();
					}
					
					// Attempt to level the playing field.
					error_behind *= bpm / (double) doff_behind.size();
				}
				
				if( doff_ahead.size() >= 2 )
				{
					for( std::vector<int>::iterator doff_iter = doff_ahead.begin(); doff_iter != doff_ahead.end(); doff_iter ++ )
						mean_ahead += *doff_iter / (double) doff_ahead.size();
					
					error_ahead = 0.;
					
					for( std::vector<int>::iterator doff_iter = doff_ahead.begin(); doff_iter != doff_ahead.end(); doff_iter ++ )
					{
						double ddoff = (*doff_iter - mean_ahead);
						error_ahead += (ddoff * ddoff) / (double) doff_ahead.size();
					}
					
					// Attempt to level the playing field.
					error_behind *= bpm / (double) doff_ahead.size();
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


class SongLoadData
{
public:
	bool Finished;
	time_t StartSecs;
	UserData *UD;
	std::string Filename;
	Song *LoadingSong;
	
	SongLoadData( UserData *ud, std::string filename )
	{
		Finished = false;
		StartSecs = time(NULL);
		UD = ud;
		Filename = filename;
		LoadingSong = NULL;
	}
	
	~SongLoadData(){}
};


// --------------------------------------------------------------------------------------


class UserData
{
public:
	bool Playing;
	SDL_AudioSpec Spec;
	double BPM;
	bool SourceBPM;
	double SourcePitchScale;
	uint8_t Resample;
	double Volume;
	bool Repeat;
	bool Metronome;
	PlaybackBuffer Buffer;
	FILE *WriteTo;
	size_t WroteBytes;
	std::deque<std::string> Queue;
	std::deque<Song*> Songs;
	std::map<SDL_Thread*,SongLoadData*> LoadThreads;
	
	UserData( void )
	{
		Playing = false;
		memset( &Spec, 0, sizeof(Spec) );
		//memset( &Info, 0, sizeof(Info) );
		BPM = 140.;
		SourceBPM = true;
		SourcePitchScale = 1.;
		Resample = ResampleMethod::Auto;
		Volume = 1.;
		Repeat = true;
		Metronome = false;
		Buffer.SetSize( 131072 );
		WriteTo = NULL;
		WroteBytes = 0;
	}
	
	void LoadSong( std::string filename )
	{
		SongLoadData *sld = new SongLoadData( this, filename );
		SDL_Thread *thread = SDL_CreateThread( &SongLoad, sld );
		if( thread )
			LoadThreads[ thread ] = sld;
		else
			delete sld;
	}
	
	void CheckThreads( void )
	{
		for( std::map<SDL_Thread*,SongLoadData*>::iterator thread_iter = LoadThreads.begin(); thread_iter != LoadThreads.end(); thread_iter ++ )
		{
			if( thread_iter->second->Finished )
			{
				// Finished loading.
				
				int unused = 0;
				SDL_WaitThread( thread_iter->first, &unused );
				
				// Make sure the song loaded okay.
				if( thread_iter->second->LoadingSong )
				{
					printf( "%s: %.2f BPM, Start %.4f sec, Beats %.1f\n", thread_iter->second->Filename.c_str(), thread_iter->second->LoadingSong->BPM, thread_iter->second->LoadingSong->FirstBeatFrame / thread_iter->second->LoadingSong->Audio.SampleRate, thread_iter->second->LoadingSong->TotalBeats() );
					fflush( stdout );
					
					if( thread_iter->second->UD->Repeat )
						thread_iter->second->UD->Queue.push_back( thread_iter->second->Filename );
					
					bool first_song = thread_iter->second->UD->Songs.empty();
					
					// Only the first song should play pre-beat lead-in.
					if( ! first_song )
						thread_iter->second->LoadingSong->CurrentFrame = thread_iter->second->LoadingSong->FirstBeatFrame;
					
					// Done loading and analyzing, so put it in the playback queue.
					thread_iter->second->UD->Songs.push_back( thread_iter->second->LoadingSong );
					
					// Start playback after we load the first one successfully.
					if( first_song )
					{
						thread_iter->second->UD->Playing = true;
						SDL_PauseAudio( 0 );
					}
				}
				
				delete thread_iter->second;
				LoadThreads.erase( thread_iter );
				break;
			}
		}
	}
	
	~UserData(){}
};


// --------------------------------------------------------------------------------------


int SongLoad( void *data_ptr )
{
	SongLoadData *data = (SongLoadData*) data_ptr;
	
	// Wait 10ms before beginning to load song.
	SDL_Delay( 10 );
	
	data->LoadingSong = new Song( data->Filename, &(data->UD->Spec) );
	
	// Make sure it loaded okay.
	if( data->LoadingSong->Audio.Data && data->LoadingSong->Audio.Size )
	{
		// Wait 10ms before beginning to analyze song.
		SDL_Delay( 10 );
		
		data->LoadingSong->Analyze();
	}
	else
	{
		// It didn't load correctly, so we won't queue it.  Memory cleanup.
		delete data->LoadingSong;
		data->LoadingSong = NULL;
	}
	
	// Notify main thread that it should now clean up memory.
	data->Finished = true;
	
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


bool CalculateCrossfade( Song *current_song, const Song *next_song, double *crossfade_for_beats, double *crossfade_at_beat )
{
	if( !( current_song && next_song ) )
		return false;
	
	double current_total_beats = current_song->TotalBeats();
	
	// Default to 96 beats for crossfade, but allow override from the songs.
	if( (current_song->OutroBeats > 0) && (next_song->IntroBeats > 0) )
		*crossfade_for_beats = std::min<int>( current_song->OutroBeats, next_song->IntroBeats );
	else if( current_song->OutroBeats > 0 )
		*crossfade_for_beats = current_song->OutroBeats;
	else if( next_song->IntroBeats > 0 )
		*crossfade_for_beats = next_song->IntroBeats;
	else
		*crossfade_for_beats = 96.;
	
	// If FirstOutroFrame is not specified, guess a good crossfade point.
	if( current_song->FirstOutroFrame )
		*crossfade_at_beat = current_song->BeatAtFrame( current_song->FirstOutroFrame );
	else
	{
		double expected_beat = current_total_beats - fmod( current_total_beats, 16. ) - *crossfade_for_beats;
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
	if( ! ud->Playing )
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
						stream16[ i + channel ] = Finalize( current_song->NearestFrame( channel ) * ud->Volume );
				}
				else
				{
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( current_song->CubicFrame( channel ) * ud->Volume );
				}
				
				current_song->Advance( bpm );
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
					double a = a_nearest ? current_song->NearestFrame( channel ) : current_song->CubicFrame( channel );
					double b = b_nearest ? next_song->NearestFrame( channel ) : next_song->CubicFrame( channel );
					stream16[ i + channel ] = Finalize( EqualPowerCrossfade( a, b, crossfade ) * ud->Volume );
				}
				
				current_song->Advance( bpm );
				next_song->Advance( bpm );
				
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


void BufferedAudioCallback( void *userdata, Uint8* stream, int len )
{
	UserData *ud = (UserData*) userdata;
	ud->Buffer.FillStream( userdata, stream, len );
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


int main( int argc, char **argv )
{
	// If no arguments were given, search the default music directory.
	if( argc == 1 )
	{
		char cmd[ 102400 ] = "";
		#ifdef WIN32
			const char *DEFAULT_ARGS = "\"M:\\iTunes\\iTunes\\ Music\\Trance\" \"M:\\BeatPort AIFF\"";
		#else
			const char *DEFAULT_ARGS = "~/Music/iTunes/iTunes\\ Music/*Trance* /Volumes/Media/Music/iTunes/iTunes\\ Music/*Trance* /Volumes/Media/Music/BeatPort\\ AIFF &";
		#endif
		snprintf( cmd, 102400, "\"%s\" %s", argv[ 0 ], DEFAULT_ARGS );
		system( cmd );
		return 0;
	}
	
	UserData userdata;
	bool shuffle = true;
	bool window = true;
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
			if( strcasecmp( argv[ i ], "--metronome" ) == 0 )
				userdata.Metronome = true;
			else if( strcasecmp( argv[ i ], "--no-shuffle" ) == 0 )
				shuffle = false;
			else if( strcasecmp( argv[ i ], "--no-repeat" ) == 0 )
				userdata.Repeat = false;
			else if( strcasecmp( argv[ i ], "--no-window" ) == 0 )
				window = false;
			else if( strcasecmp( argv[ i ], "--bpm=source" ) == 0 )
				userdata.SourceBPM = true;
			else if( strncasecmp( argv[ i ], "--bpm=", strlen("--bpm=") ) == 0 )
			{
				userdata.BPM = atof( argv[ i ] + strlen("--bpm=") );
				userdata.SourceBPM = false;
			}
			else if( strncasecmp( argv[ i ], "--pitch=", strlen("--pitch=") ) == 0 )
			{
				const char *pitch = argv[ i ] + strlen("--pitch=");
				if( pitch[ 0 ] == '+' )
					pitch ++;
				userdata.SourcePitchScale = pow( 2., atof(pitch) / 12. );
				userdata.SourceBPM = true;
			}
			else if( strncasecmp( argv[ i ], "--volume=", strlen("--volume=") ) == 0 )
				userdata.Volume = atof( argv[ i ] + strlen("--volume=") );
			else if( strncasecmp( argv[ i ], "--rate=", strlen("--rate=") ) == 0 )
				want.freq = atoi( argv[ i ] + strlen("--rate=") );
			else if( strncasecmp( argv[ i ], "--buffer1=", strlen("--buffer1=") ) == 0 )
				want.samples = atoi( argv[ i ] + strlen("--buffer1=") );
			else if( strncasecmp( argv[ i ], "--buffer2=", strlen("--buffer2=") ) == 0 )
			{
				userdata.Buffer.SetSize( atoi( argv[ i ] + strlen("--buffer2=") ) * 2 * want.channels );
				if( userdata.Buffer.BufferSize > 0 )
					want.callback = BufferedAudioCallback;
				else
					want.callback = AudioCallback;
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
			
			// Get all playable files from within the selected directory, or add selected songs directly.
			std::deque<std::string> songs = DirSongs( argv[ i ] );
			for( std::deque<std::string>::const_iterator song_iter = songs.begin(); song_iter != songs.end(); song_iter ++ )
				userdata.Queue.push_back( *song_iter );
		}
	}
	
	// Seed the random number generator (for shuffle).
	srand( time(NULL) );
	
	// If we want to capture input events, we'll need to initialize SDL video.
	Uint32 sdl_flags = SDL_INIT_AUDIO;
	if( window )
		sdl_flags |= SDL_INIT_VIDEO;
	
	// Initialize SDL and open the audio output.
	SDL_Init( sdl_flags );
	if( window )
	{
		SDL_WM_SetCaption( "Raptor007's AutoDJ", "AutoDJ" );
		SDL_SetVideoMode( 200, 50, 0, 0 );
	}
	userdata.Buffer.Callback = AudioCallback;
#ifdef USE_SDL_MIXER
	Mix_OpenAudio( want.freq, want.format, want.channels, want.samples );
	SDL_CloseAudio();
#endif
	SDL_OpenAudio( &want, &(userdata.Spec) );
	
	// Shuffle song list if requested.
	if( shuffle )
		std::random_shuffle( userdata.Queue.begin(), userdata.Queue.end() );
	
	if( write )
	{
		// Write WAV header.
		unsigned char wave_header[ 44 ] = { 'R','I','F','F', 36,0xFF,0xFF,0x7F, 'W','A','V','E', 'f','m','t',0x20, 16,0,0,0, 1,0, userdata.Spec.channels,0, userdata.Spec.freq%256,userdata.Spec.freq/256,0,0, (userdata.Spec.freq*userdata.Spec.channels*2)%256,(userdata.Spec.freq*userdata.Spec.channels*2)/256,0,0, 4,0, 16,0, 'd','a','t','a', 0,0xFF,0xFF,0x7F };
		userdata.WriteTo = fopen( write, "wb" );
		fwrite( wave_header, 1, 44, userdata.WriteTo );
		fflush( userdata.WriteTo );
	}
	
	// Keep running until playback is complete.
	bool running = (userdata.Queue.size() || userdata.Songs.size());
	while( running )
	{
		if( userdata.Queue.size() && ((userdata.Songs.size() + userdata.LoadThreads.size()) < 3) )
		{
			std::string song_name = userdata.Queue.front();
			userdata.Queue.pop_front();
			
			userdata.LoadSong( song_name );
		}
		
		userdata.CheckThreads();
		
		if( window )
		{
			SDL_Event event;
			while( SDL_PollEvent( &event ) )
			{
				if( event.type == SDL_QUIT )
				{
					running = false;
					SDL_PauseAudio( 1 );
				}
				else if( event.type == SDL_KEYDOWN )
				{
					Uint8 key = event.key.keysym.sym;
					if( key == SDLK_SPACE )
					{
						userdata.Playing = ! userdata.Playing;
						SDL_PauseAudio( userdata.Playing ? 0 : 1 );
					}
					else if( key == SDLK_m )
						userdata.Metronome = ! userdata.Metronome;
					else if( key == SDLK_LEFTBRACKET )
						userdata.Songs.front()->CurrentFrame = 0;
					else if( key == SDLK_RIGHTBRACKET )
					{
						Song *current_song = userdata.Songs.front();
						int beat = current_song->Beat();
						double expected_beat = beat - (beat % 16);
						double actual_beat = current_song->NearestBeatAtBeat( expected_beat );
						current_song->FirstOutroFrame = current_song->FrameAtBeat( (fabs( actual_beat - expected_beat ) < 32.) ? actual_beat : expected_beat );
						current_song->OutroBeats = std::min<int>( 64, current_song->TotalBeats() - beat );
					}
					else if( key == SDLK_BACKSLASH )
						userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->MaxFrame;
					else if( key == SDLK_SEMICOLON )
					{
						int beat = userdata.Songs.front()->Beat();
						userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->FrameAtBeat( beat - (beat % 16) - 32 );
					}
					else if( key == SDLK_QUOTE )
					{
						int beat = userdata.Songs.front()->Beat();
						userdata.Songs.front()->CurrentFrame = userdata.Songs.front()->FrameAtBeat( beat - (beat % 16) + 32 );
					}
					else if( key == SDLK_COMMA )
					{
						if( userdata.SourceBPM && userdata.Songs.size() )
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
						
						userdata.BPM -= 1.;
						userdata.SourceBPM = false;
						printf( "Playback BPM: %.2f\n", userdata.BPM );
					}
					else if( key == SDLK_PERIOD )
					{
						if( userdata.SourceBPM && userdata.Songs.size() )
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
						
						userdata.BPM += 1.;
						userdata.SourceBPM = false;
						printf( "Playback BPM: %.2f\n", userdata.BPM );
					}
					else if( key == SDLK_SLASH )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale = 1.;
						printf( "Playback BPM: Match Source\n" );
					}
					else if( key == SDLK_k )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale /= pow( 2., 1./12. );
						printf( "Pitch: %.1f%%\n", userdata.SourcePitchScale * 100. );
					}
					else if( key == SDLK_l )
					{
						userdata.SourceBPM = true;
						userdata.SourcePitchScale *= pow( 2., 1./12. );
						printf( "Pitch: %.1f%%\n", userdata.SourcePitchScale * 100. );
					}
					else if( key == SDLK_MINUS )
					{
						userdata.Volume -= 0.125;
						printf( "Playback Volume: %.3f\n", userdata.Volume );
					}
					else if( key == SDLK_EQUALS )
					{
						userdata.Volume += 0.125;
						printf( "Playback Volume: %.3f\n", userdata.Volume );
					}
					else if( key == SDLK_a )
					{
						userdata.Resample = ResampleMethod::Auto;
						printf( "Resample: Auto\n" );
					}
					else if( key == SDLK_s )
					{
						if( userdata.Resample == ResampleMethod::Nearest )
						{
							userdata.Resample = ResampleMethod::Cubic;
							printf( "Resample: Cubic\n" );
						}
						else
						{
							userdata.Resample = ResampleMethod::Nearest;
							printf( "Resample: Nearest\n" );
						}
					}
					else if( key == SDLK_q )
					{
						running = false;
						SDL_PauseAudio( 1 );
					}
				}
			}
		}
		
		if( userdata.Buffer.BufferSize )
		{
			// Add some data to the playback buffer and wait 10ms.
			userdata.Buffer.AddToBuffer( &userdata, 16384 );
			SDL_Delay( 10 );
		}
		else
			// Not buffering playback, so delay 100ms.
			SDL_Delay( 100 );
		
		running = running && (userdata.Queue.size() || userdata.Songs.size() || userdata.Buffer.Buffered);
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
	if( userdata.LoadThreads.size() )
	{
		printf( "Waiting for threads...\n" );
		while( userdata.LoadThreads.size() )
		{
			SDL_Delay( 100 );
			userdata.CheckThreads();
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
	printf( "Closing SDL...\n" );
	SDL_LockAudio();
	SDL_CloseAudio();
	SDL_UnlockAudio();
	SDL_Quit();
	
	printf( "Done quitting.\n" );
	
	return 0;
}
