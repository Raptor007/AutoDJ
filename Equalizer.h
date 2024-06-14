#pragma once
class EqualizerParam;
class EqualizerFFT;
class Equalizer;

#define __STDC_CONSTANT_MACROS
#define __STDC_FORMAT_MACROS
#include <cstdio>
#include <cstring>
#include <map>
extern "C" {
#include <libavcodec/avfft.h>
}


class EqualizerParam
{
public:
	std::map<float,float> FreqScale;
	
	EqualizerParam( void );
	~EqualizerParam();
	
	float LookupScale( float freq ) const;
	float GetScale( float freq ) const;
	float GetScale( size_t index, size_t frames, unsigned int rate ) const;
	float MaxScale( void ) const;
	float MinScale( void ) const;
	float AvgScale( void ) const;
};


class EqualizerFFT
{
private:
	FFTContext *Context1, *Context2;
	FFTComplex *Complex;
	
public:
	unsigned int Channels, Rate;
	size_t Frames;
	size_t BufferSize;
	float *InputCopy;
	
	EqualizerFFT( unsigned int channels, unsigned int rate, size_t frames );
	~EqualizerFFT();
	
	void Process( float *buffer, EqualizerParam *param );
};


class Equalizer
{
private:
	std::map<size_t,EqualizerFFT*> EqualizerFFTs;
	
public:
	Equalizer( void );
	~Equalizer();
	
	void Process( float *buffer, unsigned int channels, unsigned int rate, size_t frames, EqualizerParam *param );
};
