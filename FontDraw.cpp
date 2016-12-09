#include "FontDraw.h"
#include <cstdlib>

FontDraw::FontDraw( const uint8_t *font_bin )
{
	Color = 0xFFFFFFFF;
	Pixels = NULL;
	Width = 0;
	Height = 0;
	
	CharW = (font_bin[18] + 256L*font_bin[19] + 256L*256L*font_bin[20] + 256L*256L*256L*font_bin[21]) / 256;
	CharH =  font_bin[22] + 256L*font_bin[23] + 256L*256L*font_bin[24] + 256L*256L*256L*font_bin[25];
	const uint8_t *bmp_data = font_bin + (font_bin[10] + 256L*font_bin[11] + 256L*256L*font_bin[12] + 256L*256L*256L*font_bin[13]);
	
	for( size_t i = 0; i < 256; i ++ )
	{
		Char[ i ] = (bool*) malloc( CharW * CharH );
		
		for( size_t y = 0; y < CharH; y ++ )
		{
			const uint8_t *line = bmp_data + (CharH - y - 1) * CharW * 32;
			
			for( size_t x = 0; x < CharW; x ++ )
			{
				uint8_t byte = line[ ((i * CharW) + x) / 8 ];
				size_t bit = ((i * CharW) + x) % 8;
				Char[ i ][ y * CharW + x ] = byte & (1 << (7 - bit));
			}
		}
	}
}

FontDraw::~FontDraw()
{
	for( size_t i = 0; i < 256; i ++ )
	{
		free( Char[ i ] );
		Char[ i ] = NULL;
	}
}

void FontDraw::SetTarget( uint32_t *pixels, unsigned int width, unsigned int height )
{
	Pixels = pixels;
	Width = width;
	Height = height;
}

void FontDraw::Draw( int x, int y, const char *str )
{
	Draw( x, y, str, Color );
}

void FontDraw::Draw( int x, int y, const char *str, uint32_t color )
{
	int left = x;
	
	while( size_t this_char = *((unsigned char*)( str ++ )) )
	{
		switch( this_char )
		{
			case '\n':
			{
				x = left;
				y += CharH;
				break;
			}
			case '\t':
			{
				int chars_so_far = (x - left) / CharW;
				x += (4 - (chars_so_far % 4)) * CharW;
				break;
			}
			default:
			{
				for( size_t cy = 0; cy < CharH; cy ++ )
				{
					for( size_t cx = 0; cx < CharW; cx ++ )
					{
						if( Char[ this_char ][ cy * CharW + cx ] && (x + cx > 0) && (y + cy > 0) && (x + cx < Width) && (y + cy < Height) )
							Pixels[ (y + cy) * Width + x + cx ] = color;
					}
				}
				
				x += CharW;
			}
		}
	}
}
