PREFIX = /opt/local
CC = $(PREFIX)/bin/g++
O = 2
CFLAGS = -O$(O) -march=native -mfpmath=sse -fomit-frame-pointer -ftree-vectorize -fno-strict-aliasing -flto -Wall -Wextra -pedantic -Wno-narrowing
INC = $(PREFIX)/include $(PREFIX)/include/SDL
LIB = libSDLmain.a libSDL_sound.a libSDL.a libspeex.a libmodplug.dylib libmikmod.a libsmpeg.a libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a
LIBDIR = $(PREFIX)/lib
FW = Cocoa OpenGL AudioUnit AudioToolbox IOKit Carbon
ARCH = i386

AutoDJ.app: AutoDJ Info.plist
	mkdir -p AutoDJ.app/Contents/MacOS
	rsync -ax AutoDJ AutoDJ.app/Contents/MacOS/
	rsync -ax Info.plist AutoDJ.app/Contents/

AutoDJ: $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	$(CC) $(foreach arch,$(ARCH),-arch $(arch)) $(CFLAGS) $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(foreach lib,$(LIB),$(LIBDIR)/$(lib)) $(foreach fw,$(FW),-framework $(fw)) -o AutoDJ
	chmod ugo+rx AutoDJ

%.o: %.cpp
	$(CC) $(foreach arch,$(ARCH),-arch $(arch)) $(CFLAGS) -c $< $(foreach inc,$(INC),-I$(inc)) -o $@

clean:
	rm -rf *.o AutoDJ AutoDJ.app

install:
	rsync -ax AutoDJ.app /Applications/

linux:
	make AutoDJ PREFIX=/usr LIBDIR=/usr/lib64 FW= ARCH= LIB="libSDLmain.a libSDL_sound.so libSDL.so"
	rm -rf AutoDJ.app
