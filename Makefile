PREFIX = /usr
CC = $(PREFIX)/bin/g++
O = 2
CFLAGS = -O$(O) -march=native -mfpmath=sse -ftree-vectorize -fno-strict-aliasing -flto -Wall -Wextra -pedantic -Wno-narrowing -Wno-deprecated-declarations
INC = $(PREFIX)/include $(PREFIX)/include/SDL
LIB = libSDLmain.a libSDL.so ../local/lib/libavdevice.a ../local/lib/libavformat.a ../local/lib/libavfilter.a ../local/lib/libavcodec.a ../local/lib/libavresample.a ../local/lib/libswscale.a ../local/lib/libavutil.a ../local/lib/libswresample.a libdl.so liblzma.so libbz2.a libz.a
LIBDIR = $(PREFIX)/lib
ifeq ($(shell uname -m), x86_64)
LIBDIR = $(PREFIX)/lib64
endif
LIBRARIES = $(foreach lib,$(LIB),$(LIBDIR)/$(lib))
FRAMEWORKS = $(foreach fw,$(FW),-framework $(fw))
MISC_OBJ = FontBin.o
TARGET = autodj
INSTALL_DIR = /usr/local/bin
UNAME = $(shell uname)

ifeq ($(UNAME), Darwin)
PREFIX = /opt/local
LIBDIR = $(PREFIX)/lib
LIB = libSDLmain.a libSDL.a libavdevice.a libavformat.a libavfilter.a libavcodec.a libavresample.a libswscale.a libavutil.a libswresample.a libbluray.a libfreetype.a libpng.a libxml2.a libiconv.a libfaac.a libfdk-aac.a libmp3lame.a libopus.a libschroedinger-1.0.a libtheoradec.a libtheora.a libvorbisenc.a libvorbis.a libvpx.a libx264.a libxvidcore.a libgmp.a libspeex.a libmodplug.dylib libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a libopenjp2.dylib libsoxr.dylib liborc-0.4.0.dylib libgnutls.dylib
FW = Cocoa VideoDecodeAcceleration OpenGL CoreVideo CoreFoundation AudioUnit AudioToolbox IOKit Carbon
ARCH = i386
CFLAGS += $(foreach arch,$(ARCH),-arch $(arch))
TARGET = AutoDJ.app
INSTALL_DIR = /Applications
endif

ifneq (,$(findstring MINGW,$(UNAME))$(findstring CYGWIN,$(UNAME)))
PREFIX = /mingw
CFLAGS += -DWAVEOUT -Wno-unused-parameter -mwindows -static-libgcc -static-libstdc++
LIBDIR = $(PREFIX)/lib
LIB = libiconv.a
LIBRARIES := -lmingw32 -lSDLmain -lSDL -lavdevice -lavformat -lavfilter -lavcodec -lavresample -lswscale -lavutil -lswresample $(LIBRARIES) -lws2_32 -lsecur32 -lwinmm
TARGET = AutoDJ.exe
MISC_OBJ += AutoDJ.res
INSTALL_DIR = /usr/bin
endif

.PHONY: default clean install

default: $(TARGET)

AutoDJ.app: AutoDJ Info.plist
	mkdir -p $@/Contents/MacOS
	cp AutoDJ $@/Contents/MacOS/
	rsync -ax Info.plist $@/Contents/
	mkdir -p $@/Contents/MacOS/Resources
	rsync -ax AutoDJ.icns $@/Contents/Resources/
	-cp $(PREFIX)/lib/libgcc/libstdc++.6.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libstdc++.6.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libstdc++.6.dylib "@executable_path/libstdc++.6.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libgcc/libgcc_s.1.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libgcc_s.1.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libgcc_s.1.dylib "@executable_path/libgcc_s.1.dylib" $@/Contents/MacOS/AutoDJ
	-codesign -s $(shell whoami) $@

AutoDJ.res: AutoDJ.rc
	windres $< -O coff -o $@

FontBin.o: FontBin.s Font1.bmp Font2.bmp
	$(CC) $(CFLAGS) -g -c $< -o $@

autodj AutoDJ AutoDJ.exe: $(patsubst %.cpp,%.o,$(wildcard *.cpp)) $(MISC_OBJ)
	$(CC) $(CFLAGS) $^ $(LIBRARIES) $(FRAMEWORKS) -o $@
	chmod ugo+rx $@

%.o: %.cpp $(wildcard *.h)
	$(CC) $(CFLAGS) -g -c $< $(foreach inc,$(INC),-I$(inc)) -o $@

clean:
	rm -rf *.o autodj AutoDJ AutoDJ.app

install:
	rsync -ax $(TARGET) $(INSTALL_DIR)/
