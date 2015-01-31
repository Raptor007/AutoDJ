PREFIX = /opt/local
CC = $(PREFIX)/bin/g++
O = 2
INC = $(PREFIX)/include $(PREFIX)/include/SDL
LIB = libSDLmain.a libSDL_sound.a libSDL.a libspeex.a libmodplug.dylib libmikmod.a libsmpeg.a libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a
LIBDIR = $(PREFIX)/lib
FW = Cocoa OpenGL AudioUnit AudioToolbox IOKit Carbon

AutoDJ.app: AutoDJ Info.plist
	mkdir -p AutoDJ.app/Contents/MacOS
	rsync -ax AutoDJ AutoDJ.app/Contents/MacOS/
	rsync -ax Info.plist AutoDJ.app/Contents/

AutoDJ: AutoDJ.o PlaybackBuffer.o Mutex.o
	$(CC) -O$(O) AutoDJ.o PlaybackBuffer.o Mutex.o $(foreach lib,$(LIB),$(LIBDIR)/$(lib)) $(foreach fw,$(FW),-framework $(fw)) -o AutoDJ
	chmod ugo+rx AutoDJ

AutoDJ.o: AutoDJ.cpp
	$(CC) -O$(O) -c AutoDJ.cpp $(foreach inc,$(INC),-I$(inc)) -o AutoDJ.o

Mutex.o: Mutex.cpp
	$(CC) -O$(O) -c Mutex.cpp $(foreach inc,$(INC),-I$(inc)) -o Mutex.o

PlaybackBuffer.o: PlaybackBuffer.cpp
	$(CC) -O$(O) -c PlaybackBuffer.cpp $(foreach inc,$(INC),-I$(inc)) -o PlaybackBuffer.o

clean:
	rm -rf AutoDJ.o PlaybackBuffer.o Mutex.o AutoDJ AutoDJ.app

install:
	rsync -ax AutoDJ.app /Applications/

linux:
	make AutoDJ PREFIX=/usr LIBDIR=/usr/lib64 FW= LIB="libSDLmain.a libSDL_sound.so libSDL.so"
	rm -rf AutoDJ.app
