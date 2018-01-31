#include "PlaybackBuffer.h"


PlaybackBuffer::PlaybackBuffer( void )
{
	BufferSize = 0;
	Buffered = 0;
	StartAt = 0;
	LastSent = 0;
	Buffer = NULL;
}


PlaybackBuffer::~PlaybackBuffer()
{
	free( Buffer );
	Buffer = NULL;
}


void PlaybackBuffer::SetSize( int size )
{
	if( size )
	{
		if( Buffer )
			Buffer = (Uint8*) realloc( Buffer, size );
		else
			Buffer = (Uint8*) malloc( size );
	}
	else if( Buffer )
	{
		free( Buffer );
		Buffer = NULL;
	}
	
	BufferSize = Buffer ? size : 0;
	
	if( Buffered > BufferSize )
		Buffered = BufferSize;
}


void PlaybackBuffer::FillStream( void *userdata, Uint8 *stream, int len, bool lock )
{
	if( lock )
		Lock.Lock();
	
	if( ! Buffered )
		Callback( userdata, stream, len );
	else if( Buffered < len )
	{
		int was_buffered = Buffered;
		FillStream( userdata, stream, was_buffered, false );
		Callback( userdata, stream + was_buffered, len - was_buffered );
	}
	else
	{
		LastSent = StartAt;
		
		if( StartAt + len > BufferSize )
		{
			int chunk1 = BufferSize - StartAt;
			memcpy( stream, Buffer + StartAt, chunk1 );
			memcpy( stream + chunk1, Buffer, len - chunk1 );
		}
		else
			memcpy( stream, Buffer + StartAt, len );
		
		Buffered -= len;
		StartAt += len;
		StartAt %= BufferSize;
	}
	
	if( lock )
		Lock.Unlock();
}


void PlaybackBuffer::AddToBuffer( void *userdata, int len )
{
	Lock.Lock();
	
	if( len > BufferSize - Buffered )
		len = BufferSize - Buffered;
	
	if( len )
	{
		if( StartAt + Buffered >= BufferSize )
			Callback( userdata, Buffer + ((StartAt + Buffered) % BufferSize), len );
		else if( StartAt + Buffered + len > BufferSize )
		{
			int chunk1 = BufferSize - (StartAt + Buffered);
			Callback( userdata, Buffer + StartAt + Buffered, chunk1 );
			Callback( userdata, Buffer, len - chunk1 );
		}
		else
			Callback( userdata, Buffer + StartAt + Buffered, len );
		
		Buffered += len;
	}
	
	Lock.Unlock();
}


int PlaybackBuffer::Unfilled( void ) const
{
	return BufferSize - Buffered;
}
