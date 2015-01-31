/*
 *  Mutex.cpp
 */

#include "Mutex.h"
#include <cstddef>


Mutex::Mutex( void )
{
	RawMutex = SDL_CreateMutex();
}


Mutex::~Mutex( void )
{
	SDL_DestroyMutex( RawMutex );
	RawMutex = NULL;
}


bool Mutex::Lock( void )
{
	if( ! RawMutex )
		return false;
	
	return (SDL_mutexP( RawMutex ) == 0);
}


bool Mutex::Unlock( void )
{
	if( ! RawMutex )
		return false;
	
	return (SDL_mutexV( RawMutex ) == 0);
}
