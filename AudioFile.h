#define __STDC_CONSTANT_MACROS
#define __STDC_FORMAT_MACROS
extern "C" {
#include <libavutil/imgutils.h>
#include <libavutil/samplefmt.h>
#include <libavutil/timestamp.h>
#include <libavutil/opt.h>
#include <libavformat/avformat.h>
#include <libavresample/avresample.h>
}

class AudioFile
{
public:
	uint8_t *Data;
	size_t Allocated, Size;
	size_t Channels, SampleRate, BytesPerSample;
	
	AudioFile( const char *filename = NULL );
	~AudioFile();
	
	void Clear( void );
	bool AddData( uint8_t *add_data, size_t add_size );
	bool Load( const char *filename );

private:
	AVFormatContext *fmt_ctx;
	AVCodecContext *audio_dec_ctx;
	AVStream *audio_stream;
	int audio_stream_idx;
	AVFrame *frame;
	AVPacket pkt;
	int audio_frame_count;
	int decoded;
	int got_frame;
	AVAudioResampleContext *avr;
	
	bool decode_packet( void );
	bool open_codec_context( int *stream_idx, AVFormatContext *fmt_ctx, enum AVMediaType type );
};
