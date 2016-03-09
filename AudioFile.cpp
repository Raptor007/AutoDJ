#include "AudioFile.h"


AudioFile::AudioFile( const char *filename )
{
	Data = NULL;
	Clear();
	
	if( filename )
		Load( filename );
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
	BytesPerSample = 2;
	
	fmt_ctx = NULL;
	audio_dec_ctx = NULL;
	audio_stream = NULL;
	audio_stream_idx = -1;
	frame = NULL;
	audio_frame_count = 0;
	avr = NULL;
}


void AudioFile::AddData( uint8_t *add_data, size_t add_size )
{
	if( !( add_data && add_size ) )
		return;
	
	// Make sure we have enough room for the new data.
	if( Allocated < Size + add_size )
	{
		size_t need_size = add_size + Size - Allocated;
		
		// Allocate in 1MB chunks to reduce the number of times we have to do it.
		#define CHUNK_SIZE (1024*1024)
		if( need_size % CHUNK_SIZE )
			need_size += CHUNK_SIZE - (need_size % CHUNK_SIZE);
		
		Allocated += need_size;
		if( Data )
		{
			uint8_t *new_data = (uint8_t*) realloc( Data, Allocated );
			if( new_data )
				Data = new_data;
			else
			{
				free( Data );
				Data = NULL;
				Allocated = 0;
				Size = 0;
			}
		}
		else
			Data = (uint8_t*) malloc( Allocated );
	}
	
	if( Data )
	{
		// Add data to the buffer.
		memcpy( Data + Size, add_data, add_size );
		Size += add_size;
	}
}


bool AudioFile::Load( const char *filename )
{
	int ret = 0, got_frame = 0;
	
	// open input file, and allocate format context
	if( avformat_open_input( &fmt_ctx, filename, NULL, NULL ) < 0 )
		return false;
	
	// retrieve stream information
	if( avformat_find_stream_info( fmt_ctx, NULL ) < 0 )
		return false;
	
	if( open_codec_context( &audio_stream_idx, fmt_ctx, AVMEDIA_TYPE_AUDIO ) >= 0 )
		audio_stream = fmt_ctx->streams[ audio_stream_idx ];
	if( ! audio_stream )
	{
		ret = 1;
		goto end;
	}
	
	audio_dec_ctx = audio_stream->codec;
	
	// We want interleaved 16-bit in the original channel layout and sample rate.
	if( audio_dec_ctx->sample_fmt != AV_SAMPLE_FMT_S16 )
	{
		avr = avresample_alloc_context();
		av_opt_set_int( avr, "in_channel_layout",  audio_dec_ctx->channel_layout, 0 );
		av_opt_set_int( avr, "out_channel_layout", audio_dec_ctx->channel_layout, 0 );
		av_opt_set_int( avr, "in_sample_rate",     audio_dec_ctx->sample_rate,    0 );
		av_opt_set_int( avr, "out_sample_rate",    audio_dec_ctx->sample_rate,    0 );
		av_opt_set_int( avr, "in_sample_fmt",      audio_dec_ctx->sample_fmt,     0 );
		av_opt_set_int( avr, "out_sample_fmt",     AV_SAMPLE_FMT_S16,             0 );
		if( avresample_open(avr) )
		{
			avresample_free( &avr );
			avr = NULL;
		}
	}

	frame = av_frame_alloc();
	if( ! frame )
	{
		ret = AVERROR(ENOMEM);
		goto end;
	}

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
			ret = decode_packet( &got_frame );
			if( ret < 0 )
				break;
			pkt.data += ret;
			pkt.size -= ret;
		}
		while( pkt.size > 0 );
		av_packet_unref( &orig_pkt );
	}

	// flush cached frames
	pkt.data = NULL;
	pkt.size = 0;
	do
	{
		decode_packet( &got_frame );
	}
	while( got_frame );

	if( audio_stream )
	{
		enum AVSampleFormat sfmt = audio_dec_ctx->sample_fmt;
		
		if( av_sample_fmt_is_planar(sfmt) && ! avr )
			Channels = 1;
		else
			Channels = audio_dec_ctx->channels;
		
		SampleRate = audio_dec_ctx->sample_rate;
		
		if( ! avr )
			BytesPerSample = av_get_bytes_per_sample( sfmt );
		
		ret = 0;
	}

end:
	if( audio_dec_ctx )
		avcodec_close( audio_dec_ctx );
	avformat_close_input( &fmt_ctx );
	if( frame )
		av_frame_free( &frame );
	
	if( avr )
	{
		avresample_close( avr );
		avresample_free( &avr );
	}

	return !(ret < 0);
}


// --------------------------------------------------------------------------------------


int AudioFile::decode_packet( int *got_frame )
{
	int ret = 0;
	int decoded = pkt.size;

	*got_frame = 0;

	if( pkt.stream_index == audio_stream_idx )
	{
		// decode audio frame
		ret = avcodec_decode_audio4( audio_dec_ctx, frame, got_frame, &pkt );
		if( ret < 0 )
			return ret;
		
		// Some audio decoders decode only part of the packet, and have to be
		// called again with the remainder of the packet data.
		// Sample: fate-suite/lossless-audio/luckynight-partial.shn
		// Also, some decoders might over-read the packet.
		decoded = FFMIN( ret, pkt.size );
		
		if( *got_frame )
		{
			size_t unpadded_linesize = frame->nb_samples * av_get_bytes_per_sample( (AVSampleFormat) frame->format );
			audio_frame_count ++;
			
			if( avr )
			{
				// Convert to 16-bit interleaved.
				uint8_t *output = NULL;
				int out_linesize = 0;
				av_samples_alloc( &output, &out_linesize, audio_dec_ctx->channels, frame->nb_samples, AV_SAMPLE_FMT_S16, 0 );
				avresample_convert( avr, &output, 0, frame->nb_samples, frame->data, 0, frame->nb_samples );
				AddData( output, out_linesize );
				av_freep( &output );
			}
			else
				AddData( frame->extended_data[ 0 ], unpadded_linesize );
		}
	}

	return decoded;
}


int AudioFile::open_codec_context( int *stream_idx, AVFormatContext *fmt_ctx, enum AVMediaType type )
{
	int ret = 0, stream_index = 0;
	AVStream *st = NULL;
	AVCodecContext *dec_ctx = NULL;
	AVCodec *dec = NULL;
	AVDictionary *opts = NULL;

	ret = av_find_best_stream( fmt_ctx, type, -1, -1, NULL, 0 );
	if( ret < 0 )
		return ret;
	else
	{
		stream_index = ret;
		st = fmt_ctx->streams[ stream_index ];
		
		// find decoder for the stream
		dec_ctx = st->codec;
		dec = avcodec_find_decoder( dec_ctx->codec_id );
		if( ! dec )
			return AVERROR(EINVAL);
		
		// Init the decoders, without reference counting
		av_dict_set( &opts, "refcounted_frames", "0", 0 );
		
		if( (ret = avcodec_open2( dec_ctx, dec, &opts )) < 0)
			return ret;
		
		*stream_idx = stream_index;
	}

	return 0;
}
