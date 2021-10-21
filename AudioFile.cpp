#include "AudioFile.h"


AudioFile::AudioFile( const char *filename, volatile bool *running_ptr )
{
	Data = NULL;
	Clear();
	
	if( filename && running_ptr )
		Load( filename, running_ptr );
}


AudioFile::~AudioFile()
{
	Clear();
}


void AudioFile::Clear( void )
{
	Allocated = 0;
	Size = 0;
	if( Data )
		free( Data );
	Data = NULL;
	
	Channels = 1;
	SampleRate = 44100;
	BytesPerSample = 0;
	SampleFormat = AV_SAMPLE_FMT_NONE;
	
	fmt_ctx = NULL;
	audio_dec_ctx = NULL;
	audio_stream = NULL;
	audio_stream_idx = -1;
	frame = NULL;
	memset( &pkt, 0, sizeof(AVPacket) );
	audio_frame_count = 0;
	decoded = 0;
	got_frame = 0;
	avr = NULL;
	
	Tags.clear();
}


bool AudioFile::SetAllocation( size_t new_alloc )
{
	if( Allocated == new_alloc )
		return true;
	
	if( ! new_alloc )
	{
		free( Data );
		Data = NULL;
		Allocated = 0;
		return true;
	}
	
	if( Data )
	{
		uint8_t *new_data = (uint8_t*) realloc( Data, new_alloc );
		if( new_data )
		{
			Data = new_data;
			Allocated = new_alloc;
			return true;
		}
	}
	else
	{
		Data = (uint8_t*) malloc( new_alloc );
		if( Data )
		{
			Allocated = new_alloc;
			return true;
		}
	}
	
	return false;
}


bool AudioFile::AddData( uint8_t *add_data, size_t add_size )
{
	if( ! add_size )
		return true;
	
#ifdef WIN32
	// First attempt 80 minutes to better handle compilations; it's okay if this fails.
	#define FIRST_ALLOC (808*1024*1024)
	if( (! Allocated) && (add_size <= FIRST_ALLOC) && (SampleFormat == AV_SAMPLE_FMT_FLT) )
		SetAllocation( FIRST_ALLOC * 2 );
	if( (! Allocated) && (add_size <= FIRST_ALLOC) )
		SetAllocation( FIRST_ALLOC );
#endif
	
	// Make sure we have enough room for the new data.
	if( Allocated < Size + add_size )
	{
		size_t new_alloc = add_size + Size;
		
		// Additional allocations are in 32MB chunks to reduce the number of times we have to do it.
		#define CHUNK_SIZE (32*1024*1024)
		if( new_alloc % CHUNK_SIZE )
			new_alloc += CHUNK_SIZE - (new_alloc % CHUNK_SIZE);
		
		// Try to reallocate, and handle failure.
		if( ! SetAllocation( new_alloc ) )
		{
			if( Data )
			{
				fprintf( stderr, "Couldn't realloc %iMB buffer!\n", (int)(new_alloc/(1024*1024)) );
				// FIXME: Should this retain the old buffer and wait for more memory to free up?
				free( Data );
				Data = NULL;
				Allocated = 0;
				Size = 0;
			}
			else
				fprintf( stderr, "Couldn't malloc %iMB buffer!\n", (int)(new_alloc/(1024*1024)) );
			
			return false;
		}
	}
	
	if( Data && add_data )
	{
		// Add data to the buffer.
		memcpy( Data + Size, add_data, add_size );
		Size += add_size;
		return true;
	}
	
	return false;
}


bool AudioFile::Load( const char *filename, volatile bool *running_ptr )
{
	AVDictionaryEntry *tag = NULL;
	
	if( ! *running_ptr )
		goto end;
	
	// open input file, and allocate format context
	if( avformat_open_input( &fmt_ctx, filename, NULL, NULL ) < 0 )
		goto end;
	
	// retrieve stream information
	if( avformat_find_stream_info( fmt_ctx, NULL ) < 0 )
		goto end;
	
	if( OpenAudioCodecContext() )
		audio_stream = fmt_ctx->streams[ audio_stream_idx ];
	if( ! audio_stream )
		goto end;
	
	audio_dec_ctx = audio_stream->codec;
	
	frame = av_frame_alloc();
	if( ! frame )
	{
		audio_stream = NULL;
		goto end;
	}
	
	// If we didn't specify a sample format, default to 16-bit or float automatically.
	if( SampleFormat == AV_SAMPLE_FMT_NONE )
	{
		if( (audio_dec_ctx->sample_fmt == AV_SAMPLE_FMT_S16)
		||  (audio_dec_ctx->sample_fmt == AV_SAMPLE_FMT_S16P) )
			SampleFormat = AV_SAMPLE_FMT_S16;
		else
			SampleFormat = AV_SAMPLE_FMT_FLT;
	}
	
	// We want interleaved in the original channel layout and sample rate.
	if( audio_dec_ctx->sample_fmt != SampleFormat )
	{
		avr = avresample_alloc_context();
		av_opt_set_int( avr, "in_channel_layout",  audio_dec_ctx->channel_layout, 0 );
		av_opt_set_int( avr, "out_channel_layout", audio_dec_ctx->channel_layout, 0 );
		av_opt_set_int( avr, "in_sample_rate",     audio_dec_ctx->sample_rate,    0 );
		av_opt_set_int( avr, "out_sample_rate",    audio_dec_ctx->sample_rate,    0 );
		av_opt_set_int( avr, "in_sample_fmt",      audio_dec_ctx->sample_fmt,     0 );
		av_opt_set_int( avr, "out_sample_fmt",     SampleFormat,                  0 );
		if( avresample_open(avr) < 0 )
		{
			avresample_free( &avr );
			avr = NULL;
		}
	}
	
	// Keep track of the decoded audio format.
	BytesPerSample = av_get_bytes_per_sample( avr ? SampleFormat : audio_dec_ctx->sample_fmt );
	SampleRate = audio_dec_ctx->sample_rate;
	if( (! avr) && av_sample_fmt_is_planar(audio_dec_ctx->sample_fmt) )
		Channels = 1;
	else
		Channels = audio_dec_ctx->channels;
	
	// initialize packet, set data to NULL, let the demuxer fill it
	av_init_packet( &pkt );
	pkt.data = NULL;
	pkt.size = 0;
	
	// read frames from the file
	while( av_read_frame( fmt_ctx, &pkt ) >= 0 )
	{
		AVPacket orig_pkt = pkt;
		do
		{
			if( ! DecodePacket() )
				break;
			pkt.data += decoded;
			pkt.size -= decoded;
		}
		while( pkt.size > 0 );
		av_packet_unref( &orig_pkt );
		
		if( ! *running_ptr )
			goto end;
	}
	
	// flush cached frames
	pkt.data = NULL;
	pkt.size = 0;
	do
	{
		DecodePacket();
		
		if( ! *running_ptr )
			goto end;
	}
	while( got_frame );
	
	// Get metadata tags.
	while(( tag = av_dict_get( fmt_ctx->metadata, "", tag, AV_DICT_IGNORE_SUFFIX ) ))
		Tags[ tag->key ] = tag->value;
	
end:
	if( audio_dec_ctx )
		avcodec_close( audio_dec_ctx );
	if( fmt_ctx )
		avformat_close_input( &fmt_ctx );
	if( frame )
		av_frame_free( &frame );
	
	if( avr )
	{
		avresample_close( avr );
		avresample_free( &avr );
	}
	
	// Attempt to shrink the buffer to the used size, but don't throw it away if realloc fails.
	if( Size && ! SetAllocation( Size ) )
		fprintf( stderr, "Couldn't shrink to %iMB buffer!\n", (int)(Size/(1024*1024)) );
	
	return audio_stream;
}


// --------------------------------------------------------------------------------------


bool AudioFile::DecodePacket( void )
{
	decoded = 0;
	got_frame = 0;
	
	if( pkt.stream_index == audio_stream_idx )
	{
		// decode audio frame
		int ret = avcodec_decode_audio4( audio_dec_ctx, frame, &got_frame, &pkt );
		
		if( got_frame )
		{
			size_t unpadded_linesize = frame->nb_samples * audio_dec_ctx->channels * av_get_bytes_per_sample( (AVSampleFormat) frame->format );
			audio_frame_count ++;
			
			if( avr )
			{
				// Convert to interleaved in our desired sample format.
				uint8_t *output = NULL;
				int out_linesize = 0;
				av_samples_alloc( &output, &out_linesize, audio_dec_ctx->channels, frame->nb_samples, SampleFormat, 0 );
				avresample_convert( avr, &output, 0, frame->nb_samples, frame->data, 0, frame->nb_samples );
				AddData( output, out_linesize );
				av_freep( &output );
			}
			else
				AddData( frame->extended_data[ 0 ], unpadded_linesize );
		}
		
		if( ret > 0 )
		{
			// Some audio decoders decode only part of the packet, and have to be
			// called again with the remainder of the packet data.
			// Sample: fate-suite/lossless-audio/luckynight-partial.shn
			// Also, some decoders might over-read the packet.
			decoded = FFMIN( ret, pkt.size );
			return true;
		}
	}
	
	return false;
}


bool AudioFile::OpenAudioCodecContext( void )
{
	int stream_index = 0;
	AVStream *st = NULL;
	AVCodecContext *dec_ctx = NULL;
	AVCodec *dec = NULL;
	AVDictionary *opts = NULL;

	stream_index = av_find_best_stream( fmt_ctx, AVMEDIA_TYPE_AUDIO, -1, -1, NULL, 0 );
	if( stream_index < 0 )
		return false;
	else
	{
		st = fmt_ctx->streams[ stream_index ];
		
		// find decoder for the stream
		dec_ctx = st->codec;
		dec = avcodec_find_decoder( dec_ctx->codec_id );
		if( ! dec )
			return false;
		
		// Init the decoders, without reference counting
		av_dict_set( &opts, "refcounted_frames", "0", 0 );
		
		if( avcodec_open2( dec_ctx, dec, &opts ) < 0 )
			return false;
		
		audio_stream_idx = stream_index;
	}

	return true;
}
