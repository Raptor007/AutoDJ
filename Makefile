PREFIX = /usr
CC = $(PREFIX)/bin/g++
O = 2
CFLAGS = -DUSE_LIBAV -DUSE_SDL_MIXER -DUSE_EXTERNAL_FFMPEG -O$(O) -march=native -mfpmath=sse -fomit-frame-pointer -ftree-vectorize -fno-strict-aliasing -flto -Wall -Wextra -pedantic -Wno-narrowing
INC = $(PREFIX)/include $(PREFIX)/include/SDL
LIB = libSDLmain.a libSDL_mixer.so libSDL.so ../local/lib/libavdevice.a ../local/lib/libavformat.a ../local/lib/libavfilter.a ../local/lib/libavcodec.a ../local/lib/libswresample.a ../local/lib/libswscale.a ../local/lib/libavutil.a libz.a liblzma.so libbz2.a
LIBDIR = $(PREFIX)/lib
ifeq ($(shell uname -m), x86_64)
LIBDIR = $(PREFIX)/lib64
endif
LIBRARIES = $(foreach lib,$(LIB),$(LIBDIR)/$(lib))
FRAMEWORKS = $(foreach fw,$(FW),-framework $(fw))
TARGET = autodj
INSTALL_DIR = /usr/local/bin
UNAME = $(shell uname)

ifeq ($(UNAME), Darwin)
PREFIX = /opt/local
LIBDIR = $(PREFIX)/lib
LIB = libSDLmain.a libSDL_mixer.a libSDL.a libavdevice.a libavformat.a libavfilter.a libavcodec.a libswresample.a libswscale.a libavutil.a libbluray.a libfreetype.a libpng.a libxml2.a libiconv.a libfaac.a libfdk-aac.a libmp3lame.a libopenjpeg.a libopus.a libschroedinger-1.0.a libtheoradec.a libtheora.a libvorbisenc.a libvorbis.a libvpx.a libx264.a libxvidcore.a libgmp.a libspeex.a libmodplug.dylib libmikmod.a libsmpeg.a libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a liborc-0.4.0.dylib libgnutls.dylib
FW = Cocoa VideoDecodeAcceleration OpenGL CoreVideo CoreFoundation AudioUnit AudioToolbox IOKit Carbon
ARCH = i386
CFLAGS += $(foreach arch,$(ARCH),-arch $(arch))
TARGET = AutoDJ.app
INSTALL_DIR = /Applications
endif

ifneq (,$(findstring MINGW,$(UNAME))$(findstring CYGWIN,$(UNAME)))
PREFIX = /mingw
CFLAGS += -mwindows -static-libgcc -static-libstdc++
LIBDIR = $(PREFIX)/lib
LIB = libiconv.a
LIBRARIES := -lmingw32 -lSDLmain -lSDL -lSDL_mixer -lavdevice -lavformat -lavfilter -lavcodec -lswresample -lswscale -lavutil $(LIBRARIES) -lws2_32 -lsecur32
TARGET = AutoDJ.exe
INSTALL_DIR = /usr/bin
endif

.PHONY: default clean install

default: $(TARGET)

AutoDJ.app: AutoDJ Info.plist
	mkdir -p $@/Contents/MacOS
	cp AutoDJ $@/Contents/MacOS/
	rsync -ax Info.plist $@/Contents/
	-cp $(PREFIX)/lib/libgcc/libstdc++.6.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libstdc++.6.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libstdc++.6.dylib "@executable_path/libstdc++.6.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libgcc/libgcc_s.1.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libgcc_s.1.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libgcc_s.1.dylib "@executable_path/libgcc_s.1.dylib" $@/Contents/MacOS/AutoDJ

autodj AutoDJ AutoDJ.exe: $(patsubst %.cpp,%.o,$(wildcard *.cpp))
	$(CC) $(CFLAGS) $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(LIBRARIES) $(FRAMEWORKS) -o $@
	chmod ugo+rx $@

%.o: %.cpp $(wildcard *.h)
	$(CC) $(CFLAGS) -c $< $(foreach inc,$(INC),-I$(inc)) -o $@

clean:
	rm -rf *.o autodj AutoDJ AutoDJ.app

install:
	rsync -ax $(TARGET) $(INSTALL_DIR)/
