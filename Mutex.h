/*
 *  Mutex.h
 */

#pragma once
class Mutex;

#include <SDL/SDL_thread.h>


class Mutex
{
public:
	Mutex( void );
	~Mutex();
	
	bool Lock( void );
	bool Unlock( void );
	
private:
	SDL_mutex *RawMutex;
};
