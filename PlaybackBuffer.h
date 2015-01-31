#include "Mutex.h"


class PlaybackBuffer
{
public:
	void (*Callback)( void *userdata, Uint8 *stream, int len );
	int BufferSize, Buffered, StartAt;
	Uint8 *Buffer;
	Mutex Lock;
	
	PlaybackBuffer( void );
	virtual ~PlaybackBuffer();
	
	void SetSize( int size );
	
	void FillStream( void *userdata, Uint8 *stream, int len, bool lock = true );
	void AddToBuffer( void *userdata, int len );
};
