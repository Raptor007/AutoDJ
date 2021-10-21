#pragma once
class ReverbBounce;
class ReverbParam;
class Reverb;

#include <cstdio>
#include <vector>


class ReverbBounce
{
public:
	size_t FramesBack;
	float AmpScale;
	
	ReverbBounce( void );
	ReverbBounce( size_t f, float a );
	~ReverbBounce();
};


class ReverbParam
{
public:
	float SpeakerSide, SpeakerFront, SideWall, FrontWall, BackWall, Ceiling, Floor;
	float HeadWidth, BehindScale;
	float BounceEnergy;
	float TotalScale;
	std::vector<ReverbBounce> SameSide, OppoSide;
	
	ReverbParam( void );
	~ReverbParam();
	
	float SpeakerDist( void ) const;
	float BouncedDist( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite ) const;
	float AmpScale( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite ) const;
	void AddBounce( int x_bounces, int y_bounces, int z_bounces, bool up, bool opposite, float speaker, int rate );
	void Setup( unsigned int rate );
};


class Reverb
{
private:
	float *History;
	unsigned int Rate;
	
public:
	Reverb( void );
	~Reverb();
	
	void Process( float *buffer, unsigned int channels, unsigned int rate, size_t frames, ReverbParam *param );
};
