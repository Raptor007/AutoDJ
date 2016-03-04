#define __STDC_CONSTANT_MACROS
#define __STDC_FORMAT_MACROS

#ifdef USE_LIBAV
extern "C" {
#include <libavutil/imgutils.h>
#include <libavutil/samplefmt.h>
#include <libavutil/timestamp.h>
#include <libavformat/avformat.h>
}
#endif

#ifdef USE_SDL_MIXER
#include <SDL/SDL_audio.h>
#endif

class AudioFile
{
public:
	uint8_t *Data;
	size_t Allocated, Size;
	size_t Channels, SampleRate, BytesPerSample;
	const char *Error;
	
	AudioFile( const char *filename = NULL );
	~AudioFile();
	
	void Clear( void );
	void AddData( uint8_t *add_data, size_t add_size );
	
#ifdef USE_SDL_MIXER
	bool LoadWithSDLMixer( const char *filename, const SDL_AudioSpec *spec, bool use_ffmpeg = false );
#endif

#ifdef USE_LIBAV
	bool LoadWithLibAV( const char *filename );

private:
	AVFormatContext *fmt_ctx;
	AVCodecContext *audio_dec_ctx;
	int width, height;
	enum AVPixelFormat pix_fmt;
	AVStream *audio_stream;
	const char *src_filename;
	
	int audio_stream_idx;
	AVFrame *frame;
	AVPacket pkt;
	int audio_frame_count;
	
	int refcount;
	
	int decode_packet( int *got_frame );
	int open_codec_context( int *stream_idx, AVFormatContext *fmt_ctx, enum AVMediaType type );
	int get_format_from_sample_fmt( const char **fmt, enum AVSampleFormat sample_fmt );
#endif
};
