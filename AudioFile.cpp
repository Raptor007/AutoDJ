/*
 * Copyright (c) 2012 Stefano Sabatini
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#include "AudioFile.h"
#ifdef USE_SDL_MIXER
#include <SDL/SDL_mixer.h>
#endif


AudioFile::AudioFile( const char *filename )
{
	Data = NULL;
	Clear();
	
	if( filename )
		LoadWithLibAV( filename );
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
	Error = NULL;
	
#ifdef USE_LIBAV
	fmt_ctx = NULL;
	audio_dec_ctx = NULL;
	width = 0;
	height = 0;
	audio_stream = NULL;
	src_filename = NULL;
	audio_stream_idx = -1;
	frame = NULL;
	audio_frame_count = 0;
#endif
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
		Data = (uint8_t*)( Data ? realloc( Data, Allocated ) : malloc( Allocated ) );
	}
	
	// Add data to the buffer.
	memcpy( Data + Size, add_data, add_size );
	Size += add_size;
}


#ifdef USE_SDL_MIXER
bool AudioFile::LoadWithSDLMixer( const char *filename, const SDL_AudioSpec *spec, bool use_ffmpeg )
{
	Channels = spec->channels;
	SampleRate = spec->freq;
	
	Mix_Chunk *chunk = Mix_LoadWAV( filename );
#ifdef USE_EXTERNAL_FFMPEG
	if( use_ffmpeg && ! chunk )
	{
		#ifdef WIN32
			#define TEMPDIR "C:\\Windows\\Temp\\"
			#define FFMPEG "ffmpeg.exe"
		#else
			#define TEMPDIR "/tmp/"
			#define FFMPEG "ffmpeg"
		#endif
		
		const char *last_slash = strrchr( filename, '/' );
#ifdef WIN32
		if( ! last_slash )
			last_slash = strrchr( filename, '\\' );
#endif
		char new_filename[ 1024 ] = "";
		snprintf( new_filename, 1024, "%s%s.wav", TEMPDIR, last_slash ? (last_slash + 1) : filename );
		
		char cmd[ 1024*128 ] = "";
		
		snprintf( cmd, 1024*128, "%s -loglevel quiet -i \"%s\" \"%s\"", FFMPEG, filename, new_filename );
		system( cmd );
		
		chunk = Mix_LoadWAV( new_filename );
		
		#ifdef WIN32
			#define REMOVE "del"
		#else
			#define REMOVE "rm"
		#endif
		snprintf( cmd, 1024*128, "%s \"%s\"", REMOVE, new_filename );
		system( cmd );
	}
#endif
	if( chunk )
	{
		AddData( chunk->abuf, chunk->alen );
		
		Mix_FreeChunk( chunk );
		chunk = NULL;
		
		return true;
	}
	else
		Error = Mix_GetError();
	
	return false;
}
#endif


#ifdef USE_LIBAV
bool AudioFile::LoadWithLibAV( const char *filename )
{
	int ret = 0, got_frame;

	src_filename = filename;

	/* register all formats and codecs */
	av_register_all();

	/* open input file, and allocate format context */
	if (avformat_open_input(&fmt_ctx, src_filename, NULL, NULL) < 0) {
		//fprintf(stderr, "Could not open source file %s\n", src_filename);
		return false;
	}

	/* retrieve stream information */
	if (avformat_find_stream_info(fmt_ctx, NULL) < 0) {
		//fprintf(stderr, "Could not find stream information\n");
		return false;
	}

	if (open_codec_context(&audio_stream_idx, fmt_ctx, AVMEDIA_TYPE_AUDIO) >= 0) {
		audio_stream = fmt_ctx->streams[audio_stream_idx];
		audio_dec_ctx = audio_stream->codec;
	}

	if (!audio_stream) {
		//fprintf(stderr, "Could not find audio stream in the input, aborting\n");
		ret = 1;
		goto end;
	}

	frame = av_frame_alloc();
	if (!frame) {
		//fprintf(stderr, "Could not allocate frame\n");
		ret = AVERROR(ENOMEM);
		goto end;
	}

	/* initialize packet, set data to NULL, let the demuxer fill it */
	av_init_packet(&pkt);
	pkt.data = NULL;
	pkt.size = 0;

	/* read frames from the file */
	while (av_read_frame(fmt_ctx, &pkt) >= 0) {
		AVPacket orig_pkt = pkt;
		do {
			ret = decode_packet( &got_frame );
			if (ret < 0)
				break;
			pkt.data += ret;
			pkt.size -= ret;
		} while (pkt.size > 0);
		av_packet_unref(&orig_pkt);
	}

	/* flush cached frames */
	pkt.data = NULL;
	pkt.size = 0;
	do {
		decode_packet( &got_frame );
	} while (got_frame);

	if (audio_stream) {
		enum AVSampleFormat sfmt = audio_dec_ctx->sample_fmt;
		int n_channels = audio_dec_ctx->channels;
		const char *fmt;

		if (av_sample_fmt_is_planar(sfmt)) {
			av_get_sample_fmt_name(sfmt);
			sfmt = av_get_packed_sample_fmt(sfmt);
			n_channels = 1;
		}

		if ((ret = get_format_from_sample_fmt(&fmt, sfmt)) < 0)
			goto end;
		
		Channels = n_channels;
		SampleRate = audio_dec_ctx->sample_rate;
		BytesPerSample = av_get_bytes_per_sample(sfmt);
	}

end:
	avcodec_close(audio_dec_ctx);
	avformat_close_input(&fmt_ctx);
	av_frame_free(&frame);

	return !(ret < 0);
}


// --------------------------------------------------------------------------------------


int AudioFile::decode_packet( int *got_frame )
{
	int ret = 0;
	int decoded = pkt.size;

	*got_frame = 0;

	if (pkt.stream_index == audio_stream_idx) {
		/* decode audio frame */
		ret = avcodec_decode_audio4(audio_dec_ctx, frame, got_frame, &pkt);
		if (ret < 0) {
			//fprintf(stderr, "Error decoding audio frame (%s)\n", av_err2str(ret));
			return ret;
		}
		/* Some audio decoders decode only part of the packet, and have to be
		 * called again with the remainder of the packet data.
		 * Sample: fate-suite/lossless-audio/luckynight-partial.shn
		 * Also, some decoders might over-read the packet. */
		decoded = FFMIN(ret, pkt.size);

		if (*got_frame) {
			size_t unpadded_linesize = frame->nb_samples * av_get_bytes_per_sample((AVSampleFormat)frame->format);
			audio_frame_count ++;

			/* Write the raw audio data samples of the first plane. This works
			 * fine for packed formats (e.g. AV_SAMPLE_FMT_S16). However,
			 * most audio decoders output planar audio, which uses a separate
			 * plane of audio samples for each channel (e.g. AV_SAMPLE_FMT_S16P).
			 * In other words, this code will write only the first audio channel
			 * in these cases.
			 * You should use libswresample or libavfilter to convert the frame
			 * to packed data. */
			AddData( frame->extended_data[0], unpadded_linesize );
		}
	}

	return decoded;
}


int AudioFile::open_codec_context( int *stream_idx, AVFormatContext *fmt_ctx, enum AVMediaType type )
{
	int ret, stream_index;
	AVStream *st;
	AVCodecContext *dec_ctx = NULL;
	AVCodec *dec = NULL;
	AVDictionary *opts = NULL;

	ret = av_find_best_stream(fmt_ctx, type, -1, -1, NULL, 0);
	if (ret < 0) {
		//fprintf(stderr, "Could not find %s stream in input file '%s'\n",
		//		av_get_media_type_string(type), src_filename);
		return ret;
	} else {
		stream_index = ret;
		st = fmt_ctx->streams[stream_index];

		/* find decoder for the stream */
		dec_ctx = st->codec;
		dec = avcodec_find_decoder(dec_ctx->codec_id);
		if (!dec) {
			//fprintf(stderr, "Failed to find %s codec\n",
			//		av_get_media_type_string(type));
			return AVERROR(EINVAL);
		}

		/* Init the decoders, with or without reference counting */
		av_dict_set(&opts, "refcounted_frames", "0", 0);
		if ((ret = avcodec_open2(dec_ctx, dec, &opts)) < 0) {
			//fprintf(stderr, "Failed to open %s codec\n",
			//		av_get_media_type_string(type));
			return ret;
		}
		*stream_idx = stream_index;
	}

	return 0;
}

int AudioFile::get_format_from_sample_fmt( const char **fmt, enum AVSampleFormat sample_fmt )
{
	size_t i;
	struct sample_fmt_entry {
		enum AVSampleFormat sample_fmt; const char *fmt_be, *fmt_le;
	} sample_fmt_entries[] = {
		{ AV_SAMPLE_FMT_U8,  "u8",    "u8"    },
		{ AV_SAMPLE_FMT_S16, "s16be", "s16le" },
		{ AV_SAMPLE_FMT_S32, "s32be", "s32le" },
		{ AV_SAMPLE_FMT_FLT, "f32be", "f32le" },
		{ AV_SAMPLE_FMT_DBL, "f64be", "f64le" },
	};
	*fmt = NULL;

	for (i = 0; i < FF_ARRAY_ELEMS(sample_fmt_entries); i++) {
		struct sample_fmt_entry *entry = &sample_fmt_entries[i];
		if (sample_fmt == entry->sample_fmt) {
			*fmt = AV_NE(entry->fmt_be, entry->fmt_le);
			return 0;
		}
	}

	//fprintf(stderr,
	//		"sample format %s is not supported as output format\n",
	//		av_get_sample_fmt_name(sample_fmt));
	return -1;
}
#endif
