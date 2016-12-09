#include <stdint.h>
#include <cstddef>

class FontDraw
{
public:
	uint32_t Color;
	uint32_t *Pixels;
	unsigned int Width, Height;
	size_t CharW, CharH;
	bool *Char[ 256 ];
	
	FontDraw( const uint8_t *font_bin );
	~FontDraw();
	
	void SetTarget( uint32_t *pixels, unsigned int width, unsigned int height );
	void Draw( int x, int y, const char *str );
	void Draw( int x, int y, const char *str, uint32_t color );
};
