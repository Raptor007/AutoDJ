PREFIX = /usr
CC = $(PREFIX)/bin/g++
O = 2
CFLAGS = -O$(O) -march=native -mfpmath=sse -fomit-frame-pointer -ftree-vectorize -fno-strict-aliasing -flto -Wall -Wextra -pedantic -Wno-narrowing
INC = $(PREFIX)/include $(PREFIX)/include/SDL
LIB = libSDLmain.a libSDL_mixer.so libSDL.so
LIBDIR = $(PREFIX)/lib64
LIBRARIES = $(foreach lib,$(LIB),$(LIBDIR)/$(lib))
FRAMEWORKS = $(foreach fw,$(FW),-framework $(fw))
TARGET = autodj
INSTALL_DIR = /usr/local/bin
UNAME = $(shell uname)

ifeq ($(UNAME), Darwin)
PREFIX = /opt/local
LIBDIR = $(PREFIX)/lib
LIB = libSDLmain.a libSDL_mixer.a libSDL.a libspeex.a libmodplug.dylib libmikmod.a libsmpeg.a libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a
FW = Cocoa OpenGL AudioUnit AudioToolbox IOKit Carbon
ARCH = i386
CFLAGS += $(foreach arch,$(ARCH),-arch $(arch))
TARGET = AutoDJ.app
INSTALL_DIR = /Applications
endif

ifneq (,$(findstring MINGW,$(UNAME))$(findstring CYGWIN,$(UNAME)))
PREFIX = /mingw
CFLAGS += -mwindows -static-libgcc -static-libstdc++
LIBRARIES = -lmingw32 -lSDLmain -lSDL -lSDL_mixer
TARGET = AutoDJ.exe
INSTALL_DIR = /usr/bin
endif

.PHONY: default clean install

default: $(TARGET)

AutoDJ.app: AutoDJ Info.plist
	mkdir -p AutoDJ.app/Contents/MacOS
	rsync -ax AutoDJ AutoDJ.app/Contents/MacOS/
	rsync -ax Info.plist AutoDJ.app/Contents/

autodj AutoDJ AutoDJ.exe: $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	$(CC) $(CFLAGS) $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(LIBRARIES) $(FRAMEWORKS) -o $@
	chmod ugo+rx $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< $(foreach inc,$(INC),-I$(inc)) -o $@

clean:
	rm -rf *.o autodj AutoDJ AutoDJ.app

install:
	rsync -ax $(TARGET) $(INSTALL_DIR)/
