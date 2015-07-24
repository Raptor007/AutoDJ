#include <SDL/SDL.h>
#include <SDL/SDL_audio.h>
#include <SDL/SDL_sound.h>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <stdint.h>
#include <dirent.h>
#include <sys/stat.h>
#include <deque>
#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include "PlaybackBuffer.h"


class Song;
class SongLoadData;
class UserData;

void AudioCallback( void *userdata, Uint8* stream, int len );
void BufferedAudioCallback( void *userdata, Uint8* stream, int len );

int SongLoad( void *data_ptr );


namespace ResampleMethod
{
	enum
	{
		Auto = 0,
		Nearest,
		Linear,
		Triangular
	};
}


class Song
{
public:
	double BPM;
	Sound_Sample *Sample;
	double CurrentFrame, MaxFrame, FirstBeatFrame, FirstOutroFrame;
	int IntroBeats, OutroBeats;
	std::map<size_t,double> Beats;
	
	Song( std::string filename, Sound_AudioInfo *info )
	{
		Sample = Sound_NewSampleFromFile( filename.c_str(), info, 1024*1024*32 );
		if( Sample )
		{
			Sound_DecodeAll( Sample );
			
			// For simplicity's sake, assume it's always 16-bit audio.
			MaxFrame = Sample->buffer_size / (2 * Sample->actual.channels);
			
			if( MaxFrame < 0.5 )
			{
				Sound_FreeSample( Sample );
				Sample = NULL;
				printf( "%s: zero-length\n", filename.c_str() );
			}
		}
		else
		{
			MaxFrame = 0.;
			printf( "%s: %s\n", filename.c_str(), Sound_GetError() );
		}
		
		CurrentFrame = 0.;
		BPM = 140.;
		FirstBeatFrame = 0.;
		IntroBeats = 0;
		OutroBeats = 0;
		FirstOutroFrame = 0.;
	}
	
	virtual ~Song()
	{
		if( Sample )
			Sound_FreeSample( Sample );
		Sample = NULL;
	}
	
	void SetFirstBeat( int minutes, double seconds )
	{
		if( Sample )
			FirstBeatFrame = Sample->actual.rate * ((60. * minutes) + seconds);
		else
			FirstBeatFrame = 0.;
	}
	
	double NearestFrame( Uint8 channel ) const
	{
#ifdef __GNUC__
		if(__builtin_expect( !!( Sample && (channel < Sample->actual.channels) ), 1 ))
#else
		if( Sample && (channel < Sample->actual.channels) )
#endif
		{
			Sint16 *buffer16 = (Sint16*) Sample->buffer;
			size_t a_index = Sample->actual.channels * (size_t) CurrentFrame;
			size_t b_index = a_index + Sample->actual.channels;
#ifdef __GNUC__
			Sint16 a = __builtin_expect( !!( a_index < Sample->buffer_size / 2 ), 1 ) ? buffer16[ a_index ] : 0;
			Sint16 b = __builtin_expect( !!( b_index < Sample->buffer_size / 2 ), 1 ) ? buffer16[ b_index ] : 0;
#else
			Sint16 a = (a_index < Sample->buffer_size / 2) ? buffer16[ a_index ] : 0;
			Sint16 b = (b_index < Sample->buffer_size / 2) ? buffer16[ b_index ] : 0;
#endif
			double unused = 0.;
			double b_part = modf( CurrentFrame, &unused );
			return (b_part >= 0.5) ? b : a;
		}
		else
			return 0.;
	}
	
	double InterpolatedFrame( Uint8 channel ) const
	{
#ifdef __GNUC__
		if(__builtin_expect( !!( Sample && (channel < Sample->actual.channels) ), 1 ))
#else
		if( Sample && (channel < Sample->actual.channels) )
#endif
		{
			Sint16 *buffer16 = (Sint16*) Sample->buffer;
			size_t a_index = Sample->actual.channels * (size_t) CurrentFrame;
			size_t b_index = a_index + Sample->actual.channels;
#ifdef __GNUC__
			Sint16 a = __builtin_expect( !!( a_index < Sample->buffer_size / 2 ), 1 ) ? buffer16[ a_index ] : 0;
			Sint16 b = __builtin_expect( !!( b_index < Sample->buffer_size / 2 ), 1 ) ? buffer16[ b_index ] : 0;
#else
			Sint16 a = (a_index < Sample->buffer_size / 2) ? buffer16[ a_index ] : 0;
			Sint16 b = (b_index < Sample->buffer_size / 2) ? buffer16[ b_index ] : 0;
#endif
			double unused = 0.;
			double b_part = modf( CurrentFrame, &unused );
			return a * (1 - b_part) + b * b_part;
		}
		else
			return 0.;
	}
	
	double TriangularFrame( Uint8 channel ) const
	{
#ifdef __GNUC__
		if(__builtin_expect( !!( (CurrentFrame <= 1.) || (CurrentFrame + 2. >= MaxFrame) ), 0 ))
			return InterpolatedFrame( channel );
		else if(__builtin_expect( !!( Sample && (channel < Sample->actual.channels) ), 1 ))
#else
		if( (CurrentFrame <= 1.) || (CurrentFrame + 2. >= MaxFrame) )
			return InterpolatedFrame( channel );
		else if( Sample && (channel < Sample->actual.channels) )
#endif
		{
			Sint16 *buffer16 = (Sint16*) Sample->buffer;
			size_t a_index = Sample->actual.channels * (size_t) CurrentFrame;
			size_t b_index = a_index + Sample->actual.channels;
			size_t prev_index = a_index - Sample->actual.channels;
			size_t next_index = b_index + Sample->actual.channels;
			double a = buffer16[ a_index ];
			double b = buffer16[ b_index ];
			double prev = buffer16[ prev_index ];
			double next = buffer16[ next_index ];
			double mid = a + b - (prev / 2.) - (next / 2.);
			double unused = 0.;
			double b_part = modf( CurrentFrame, &unused );
			if( b_part <= 0.5 )
				return a * (1 - b_part * 2.) + mid * b_part * 2.;
			else
				return mid * (1 - (b_part - 0.5) * 2.) + b * (b_part - 0.5) * 2.;
		}
		else
			return 0.;
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
		if( Sample )
			return (frame - FirstBeatFrame) * BPM / (60. * Sample->actual.rate);
		else
			return 0.;
	}
	
	double FrameAtBeat( double beat ) const
	{
		if( Sample )
			return beat * (60. * Sample->actual.rate) / BPM + FirstBeatFrame;
		else
			return 0.;
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
	
	void Analyze( void )
	{
		if( Sample )
		{
			// Apply low-pass filter and keep track of average peak height over a period of samples.
			
			Sint16 *buffer16 = (Sint16*) Sample->buffer;
			
			// Don't analyze more than 10 minutes of audio per track.
			size_t max_analysis_frame = Sample->actual.rate * 60 * 10;
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
				for( size_t channel = 0; channel < Sample->actual.channels; channel ++ )
					prev += buffer16[ frame * Sample->actual.channels + channel ] * pow( prev_sample_power, samples - frame - 1. ) / (double) Sample->actual.channels;
			
			// Process the rest of the audio.
			for( size_t frame = samples; frame < max_analysis_frame; frame ++ )
			{
				// Subtract the oldest sample of the block and add the new one.
				double point = prev * prev_sample_power;
				for( size_t channel = 0; channel < Sample->actual.channels; channel ++ )
				{
					point -= buffer16[ (frame - samples) * Sample->actual.channels + channel ] * pow( prev_sample_power, samples ) / (double) Sample->actual.channels;
					point += buffer16[ frame * Sample->actual.channels + channel ] / (double) Sample->actual.channels;
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
			
			size_t best_frame_skip = Sample->actual.rate * 60 / 140;
			int best_bpm = 0;
			size_t best_first_beat = first_beat;
			double best_error = FLT_MAX;
			int best_doff = INT_MAX;
			
			for( int bpm = 120; (bpm <= 150) && (best_error > 0.); bpm ++ )
			{
				// Search at a BPM and see how closely it matches.
				
				size_t frame_skip = (44100. * 60. / (double) bpm) + 0.5;
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
						best_bpm = bpm;
						best_doff = doff_behind[ doff_behind.size() / 2 ];
						best_frame_skip = frame_skip + best_doff;
						best_first_beat = first_beat_matched ? first_beat_matched : first_beat;
						best_error = error_behind;
					}
					if( doff_ahead.size() && (error_ahead < best_error) )
					{
						best_bpm = bpm;
						best_doff = doff_ahead[ doff_ahead.size() / 2 ];
						best_frame_skip = frame_skip + best_doff;
						best_first_beat = first_beat_matched ? first_beat_matched : first_beat;
						best_error = error_ahead;
					}
				}
			}
			
			
			//printf( "Detected: %.2f BPM, Start: %.4f, Error from %i BPM: %.2f\n", 44100. * 60. / (double) best_frame_skip, first_beat / 44100., best_bpm, best_error );
			
			
			// If it's close, round to the nearest whole BPM.
			double bpm_rounding = 0.106;
			double bpm_fpart = modf( 44100. * 60. / (double) best_frame_skip, &BPM );
			if( bpm_fpart >= (1. - bpm_rounding) )
				BPM += 1.;
			else if( bpm_fpart > bpm_rounding )
				BPM += bpm_fpart;
			FirstBeatFrame = best_first_beat;
		}
	}
};


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
	
	virtual ~SongLoadData(){}
};


class UserData
{
public:
	bool Playing;
	SDL_AudioSpec Spec;
	Sound_AudioInfo Info;
	double BPM;
	bool SourceBPM;
	uint8_t Resample;
	double Volume;
	bool Repeat;
	bool Metronome;
	PlaybackBuffer Buffer;
	time_t LoadTimeoutSecs;
	FILE *WriteTo;
	size_t WroteBytes;
	std::deque<std::string> Queue;
	std::deque<Song*> Songs;
	std::map<SDL_Thread*,SongLoadData*> LoadThreads;
	
	UserData( void )
	{
		Playing = false;
		memset( &Spec, 0, sizeof(Spec) );
		memset( &Info, 0, sizeof(Info) );
		BPM = 140.;
		SourceBPM = true;
		Resample = ResampleMethod::Auto;
		Volume = 1.;
		Repeat = true;
		Metronome = false;
		Buffer.SetSize( 131072 );
		LoadTimeoutSecs = 120;
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
		time_t now_secs = time(NULL);
		
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
					printf( "%s: %.2f BPM, Start: %.4f sec\n", thread_iter->second->Filename.c_str(), thread_iter->second->LoadingSong->BPM, thread_iter->second->LoadingSong->FirstBeatFrame / thread_iter->second->LoadingSong->Sample->actual.rate );
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
			else if( now_secs > thread_iter->second->StartSecs + LoadTimeoutSecs )
			{
				// We gave it some time to try loading, so if it's not done yet, kill it.
				// This is because SDL_sound can get stuck when loading certain songs.
				
				printf( "%s: Loading Timeout\n", thread_iter->second->Filename.c_str() );
				fflush( stdout );
				
				SDL_KillThread( thread_iter->first );
				
				delete thread_iter->second;
				LoadThreads.erase( thread_iter );
				break;
			}
		}
	}
	
	virtual ~UserData(){}
};


int SongLoad( void *data_ptr )
{
	SongLoadData *data = (SongLoadData*) data_ptr;
	
	// Wait 10ms before beginning to load song.
	SDL_Delay( 10 );
	
	data->LoadingSong = new Song( data->Filename, &(data->UD->Info) );
	
	// Make sure it loaded okay.
	if( data->LoadingSong->Sample && data->LoadingSong->Sample->buffer_size )
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


void AudioCallback( void *userdata, Uint8* stream, int len )
{
	memset( stream, 0, len );
	
	UserData *ud = (UserData*) userdata;
	if( ! ud->Playing )
		return;
	
	// For simplicity's sake, assume it's always 16-bit audio.
	Sint16 *stream16 = (Sint16*) stream;
	
	bool crossfade_now = false;
	double crossfade = 0.;
	
	if( ud->Songs.size() )
	{
		Song *current_song = ud->Songs.front();
		Song *next_song = (ud->Songs.size() >= 2) ? ud->Songs[ 1 ] : NULL;
		
		double bpm = ud->SourceBPM ? current_song->BPM : ud->BPM;
		
		for( int i = 0; i < len / 2; i += ud->Spec.channels )
		{
			// Check for crossfade into second song.
			
			// FIXME: Could save some math per sample by checking if these have been set already.
			crossfade_now = false;
			crossfade = 0.;
			
			if( next_song )
			{
				double current_total_beats = current_song->TotalBeats();
				
				// Default to 64 beats for crossfade, but allow override from the songs.
				double crossfade_for_beats = 96.;
				if( current_song->OutroBeats && next_song->IntroBeats )
					crossfade_for_beats = std::min<int>( current_song->OutroBeats, next_song->IntroBeats );
				else if( current_song->OutroBeats )
					crossfade_for_beats = current_song->OutroBeats;
				else if( next_song->IntroBeats )
					crossfade_for_beats = next_song->IntroBeats;
				
				// If FirstOutroFrame is not specified, guess a good crossfade point.
				double crossfade_at_beat = current_total_beats - crossfade_for_beats;
				if( current_song->FirstOutroFrame )
					crossfade_at_beat = current_song->BeatAtFrame( current_song->FirstOutroFrame );
				else
					crossfade_at_beat = current_song->NearestBeatAtBeat( current_total_beats - fmod( current_total_beats, 16. ) - crossfade_for_beats );
				
				// Make sure the current song doesn't end before the cross-fade is complete.
				if( crossfade_for_beats > current_total_beats - crossfade_at_beat )
					crossfade_for_beats = current_total_beats - crossfade_at_beat;
				
				// Now that we've determined how crossfade should be done, see if we should do it now.
				if( current_song->Beat() >= crossfade_at_beat )
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
				
				if( (ud->Resample == ResampleMethod::Nearest) || ((ud->Resample == ResampleMethod::Auto) && (fabs(current_song->BPM - bpm) < 0.5)) )
				{
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( current_song->NearestFrame( channel ) * ud->Volume );
				}
				else if( ud->Resample == ResampleMethod::Linear )
				{
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( current_song->InterpolatedFrame( channel ) * ud->Volume );
				}
				else
				{
					for( int channel = 0; channel < ud->Spec.channels; channel ++ )
						stream16[ i + channel ] = Finalize( current_song->TriangularFrame( channel ) * ud->Volume );
				}
				
				current_song->Advance( bpm );
			}
			else
			{
				// Crossfading now.
				
				if( ud->SourceBPM )
					bpm = LinearCrossfade( current_song->BPM, next_song->BPM, crossfade );
				
				bool a_nearest = ( (ud->Resample == ResampleMethod::Nearest) || ((ud->Resample == ResampleMethod::Auto) && (fabs(current_song->BPM - bpm) < 0.5)) );
				bool b_nearest = ( (ud->Resample == ResampleMethod::Nearest) || ((ud->Resample == ResampleMethod::Auto) && (fabs(next_song->BPM - bpm) < 0.5)) );
				
				for( int channel = 0; channel < ud->Spec.channels; channel ++ )
				{
					double a, b;
					if( ud->Resample == ResampleMethod::Linear )
					{
						a = current_song->InterpolatedFrame( channel );
						b = next_song->InterpolatedFrame( channel );
					}
					else
					{
						a = a_nearest ? current_song->NearestFrame( channel ) : current_song->TriangularFrame( channel );
						b = b_nearest ? next_song->NearestFrame( channel ) : next_song->TriangularFrame( channel );
					}
					
					stream16[ i + channel ] = Finalize( EqualPowerCrossfade( a, b, crossfade ) * ud->Volume );
				}
				
				current_song->Advance( bpm );
				next_song->Advance( bpm );
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
	
	// If we completed a crossfade, remove the song we faded from.
	if( crossfade_now && (crossfade >= 1.) )
	{
		Song *front = ud->Songs.front();
		ud->Songs.pop_front();
		delete front;
	}
	
	// Remove tracks that have finished playing.
	while( ud->Songs.size() && (ud->Songs.front()->CurrentFrame > ud->Songs.front()->MaxFrame) )
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


std::deque<std::string> DirSongs( std::string path )
{
	std::deque<std::string> songs;
	
	struct stat stat_result;
	memset( &stat_result, 0, sizeof(stat_result) );
	stat( path.c_str(), &stat_result );
	
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
					
					// FIXME: Check for supported file extensions?
					songs.push_back( std::string(path) + std::string("/") + std::string(entry->d_name) );
				}
				else
				{
					// Subdirectory.
					
					std::deque<std::string> subdir_songs = DirSongs( path + std::string("/") + std::string(entry->d_name) );
					for( std::deque<std::string>::const_iterator song_iter = subdir_songs.begin(); song_iter != subdir_songs.end(); song_iter ++ )
						songs.push_back( *song_iter );
				}
			}
			
			closedir( dir );
		}
	}
	else
		songs.push_back( path );
	
	return songs;
}


int main( int argc, char **argv )
{
	// If no arguments were given, search the default music directory.
	if( argc == 1 )
	{
		char cmd[ 102400 ] = "";
		snprintf( cmd, 102400, "\"%s\" ~/Music/iTunes/iTunes\\ Music/*Trance* /Volumes/Media/Music/iTunes/iTunes\\ Music/*Trance* &", argv[ 0 ] );
		system( cmd );
		return 0;
	}
	
	UserData userdata;
	bool shuffle = true;
	bool window = true;
	const char *write = NULL;
	
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
			else if( strncasecmp( argv[ i ], "--volume=", strlen("--volume=") ) == 0 )
				userdata.Volume = atof( argv[ i ] + strlen("--volume=") );
			else if( strncasecmp( argv[ i ], "--buffer=", strlen("--buffer=") ) == 0 )
				userdata.Buffer.SetSize( atoi( argv[ i ] + strlen("--buffer=") ) );
			else if( strncasecmp( argv[ i ], "--timeout=", strlen("--timeout=") ) == 0 )
				userdata.LoadTimeoutSecs = atoi( argv[ i ] + strlen("--timeout=") );
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
			if( path[ len - 1 ] == '/' )
				path[ len - 1 ] = '\0';
			
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
	SDL_AudioSpec want;
	memset( &want, 0, sizeof(want) );
	want.freq = 44100;
	want.format = AUDIO_S16;
	want.channels = 2;
	want.samples = 4096;
	want.callback = BufferedAudioCallback;
	want.userdata = &userdata;
	userdata.Buffer.Callback = AudioCallback;
	int dev = SDL_OpenAudio( &want, &(userdata.Spec) );
	Sound_Init();
	userdata.Info.format = userdata.Spec.format;
	userdata.Info.channels = userdata.Spec.channels;
	userdata.Info.rate = userdata.Spec.freq;
	
	// Shuffle song list if requested.
	if( shuffle )
		std::random_shuffle( userdata.Queue.begin(), userdata.Queue.end() );
	
	if( write )
	{
		// Write WAV header.
		unsigned char wave_header[ 44 ] = { 'R','I','F','F', 36,0xFF,0xFF,0x7F, 'W','A','V','E', 'f','m','t',0x20, 16,0,0,0, 1,0, userdata.Info.channels,0, userdata.Info.rate%256,userdata.Info.rate/256,0,0, (userdata.Info.rate*userdata.Info.channels*2)%256,(userdata.Info.rate*userdata.Info.channels*2)/256,0,0, 4,0, 16,0, 'd','a','t','a', 0,0xFF,0xFF,0x7F };
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
					running = false;
				else if( event.type == SDL_KEYDOWN )
				{
					Uint8 key = event.key.keysym.sym;
					if( key == SDLK_SPACE )
						userdata.Playing = ! userdata.Playing;
					else if( key == SDLK_m )
						userdata.Metronome = ! userdata.Metronome;
					else if( key == SDLK_LEFTBRACKET )
						userdata.Songs.front()->CurrentFrame = 0;
					else if( key == SDLK_RIGHTBRACKET )
					{
						int beat = userdata.Songs.front()->Beat();
						userdata.Songs.front()->OutroBeats = 64;
						userdata.Songs.front()->FirstOutroFrame = userdata.Songs.front()->NearestBeatFrameAtBeat( beat - (beat % 16) );
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
						if( userdata.SourceBPM )
						{
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
							userdata.SourceBPM = false;
						}
						
						userdata.BPM -= 1.;
						printf( "Playback BPM: %.2f\n", userdata.BPM );
					}
					else if( key == SDLK_PERIOD )
					{
						if( userdata.SourceBPM )
						{
							userdata.BPM = (int)( userdata.Songs.front()->BPM + 0.5 );
							userdata.SourceBPM = false;
						}
						
						userdata.BPM += 1.;
						printf( "Playback BPM: %.2f\n", userdata.BPM );
					}
					else if( key == SDLK_SLASH )
					{
						userdata.SourceBPM = true;
						printf( "Playback BPM: Match Source\n" );
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
							userdata.Resample = ResampleMethod::Linear;
							printf( "Resample: Linear\n" );
						}
						else if( userdata.Resample == ResampleMethod::Linear )
						{
							userdata.Resample = ResampleMethod::Triangular;
							printf( "Resample: Triangular\n" );
						}
						else
						{
							userdata.Resample = ResampleMethod::Nearest;
							printf( "Resample: Nearest\n" );
						}
					}
					else if( key == SDLK_q )
						running = false;
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
	userdata.LoadTimeoutSecs = 0;
	userdata.CheckThreads();
	userdata.Queue.clear();
	userdata.Songs.clear();
	
	// Quit SDL.
	Sound_Quit();
	SDL_CloseAudio();
	SDL_Quit();
	return 0;
}
