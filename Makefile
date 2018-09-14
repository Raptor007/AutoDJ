PREFIX = /usr
CC = $(PREFIX)/bin/g++
O = 2
CFLAGS = -O$(O) -march=native -mfpmath=sse -ftree-vectorize -fno-strict-aliasing -flto -Wall -Wextra -Wno-narrowing -Wno-deprecated-declarations
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
LIB = libSDLmain.a libSDL.a libavdevice.a libavformat.a libavfilter.a libavcodec.a libavresample.a libswscale.a libavutil.a libswresample.a libbluray.a libfreetype.a libpng.a libxml2.a libiconv.a libfaac.a libfdk-aac.a libmp3lame.a libopus.a libschroedinger-1.0.a libtheoradec.a libtheora.a libvorbisenc.a libvorbis.a libvpx.a libx264.a libx265.a libxvidcore.a libgmp.a libspeex.a libmodplug.dylib libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a libopenjp2.dylib libsoxr.dylib liborc-0.4.0.dylib libgnutls.dylib
FW = Cocoa VideoDecodeAcceleration OpenGL CoreVideo CoreFoundation AudioUnit AudioToolbox IOKit Carbon
MISC_OBJ += FileDrop.om
ARCH = i386
CFLAGS += $(foreach arch,$(ARCH),-arch $(arch)) -mmacosx-version-min=10.4 -dead_strip
TARGET = AutoDJ.app
INSTALL_DIR = /Applications
GCC_VERSION = $(shell $(CC) -dumpversion)
GCC_VERSION_MAJOR = $(word 1,$(subst ., ,$(GCC_VERSION)))
GCC_VERSION_MINOR = $(word 2,$(subst ., ,$(GCC_VERSION)))
ifneq (,$(findstring 4.9.,$(GCC_VERSION)))
CFLAGS += -Wl,-no_compact_unwind
else
ifneq (4,$(GCC_VERSION_MAJOR))
CFLAGS += -Wl,-no_compact_unwind
endif
endif
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
	-cp $(PREFIX)/lib/libgcc/libstdc++.6.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libstdc++.6.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libstdc++.6.dylib "@loader_path/libstdc++.6.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libgcc/libgcc_s.1.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libgcc_s.1.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libgcc_s.1.dylib "@loader_path/libgcc_s.1.dylib" $@/Contents/MacOS/AutoDJ && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgcc/libgcc_s.1.dylib "@loader_path/libgcc_s.1.dylib" $@/Contents/MacOS/libstdc++.6.dylib
	-cp $(PREFIX)/lib/libmodplug.1.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libmodplug.1.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libmodplug.1.dylib "@loader_path/libmodplug.1.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libopenjp2.7.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libopenjp2.7.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libopenjp2.7.dylib "@loader_path/libopenjp2.7.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libsoxr.0.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libsoxr.0.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libsoxr.0.dylib "@loader_path/libsoxr.0.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/liborc-0.4.0.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/liborc-0.4.0.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/liborc-0.4.0.dylib "@loader_path/liborc-0.4.0.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libgnutls.30.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libgnutls.30.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgnutls.30.dylib "@loader_path/libgnutls.30.dylib" $@/Contents/MacOS/AutoDJ
	-cp $(PREFIX)/lib/libsnowleopardfixes.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libsnowleopardfixes.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libsnowleopardfixes.dylib "@loader_path/libsnowleopardfixes.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libp11-kit.0.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libp11-kit.0.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libp11-kit.0.dylib "@loader_path/libp11-kit.0.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libz.1.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libz.1.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libz.1.dylib "@loader_path/libz.1.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libidn.11.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libidn.11.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libidn.11.dylib "@loader_path/libidn.11.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libunistring.2.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libunistring.2.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libunistring.2.dylib "@loader_path/libunistring.2.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libtasn1.6.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libtasn1.6.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libtasn1.6.dylib "@loader_path/libtasn1.6.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libhogweed.4.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libhogweed.4.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libhogweed.4.dylib "@loader_path/libhogweed.4.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libnettle.6.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libnettle.6.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libnettle.6.dylib "@loader_path/libnettle.6.dylib" $@/Contents/MacOS/libgnutls.30.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libnettle.6.dylib "@loader_path/libnettle.6.dylib" $@/Contents/MacOS/libhogweed.4.dylib
	-cp $(PREFIX)/lib/libgmp.10.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libgmp.10.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgmp.10.dylib "@loader_path/libgmp.10.dylib" $@/Contents/MacOS/libgnutls.30.dylib
	-cp $(PREFIX)/lib/libintl.8.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libintl.8.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libintl.8.dylib "@loader_path/libintl.8.dylib" $@/Contents/MacOS/libgnutls.30.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libintl.8.dylib "@loader_path/libintl.8.dylib" $@/Contents/MacOS/libidn.11.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libintl.8.dylib "@loader_path/libintl.8.dylib" $@/Contents/MacOS/libp11-kit.0.dylib
	-cp $(PREFIX)/lib/libffi.6.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libffi.6.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libffi.6.dylib "@loader_path/libffi.6.dylib" $@/Contents/MacOS/libp11-kit.0.dylib
	-cp $(PREFIX)/lib/libgmp.10.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libgmp.10.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libgmp.10.dylib "@loader_path/libgmp.10.dylib" $@/Contents/MacOS/libhogweed.4.dylib
	-cp $(PREFIX)/lib/libiconv.2.dylib $@/Contents/MacOS/ && chmod ugo+rx $@/Contents/MacOS/libiconv.2.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libiconv.2.dylib "@loader_path/libiconv.2.dylib" $@/Contents/MacOS/libidn.11.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libiconv.2.dylib "@loader_path/libiconv.2.dylib" $@/Contents/MacOS/libintl.8.dylib && $(PREFIX)/bin/install_name_tool -change $(PREFIX)/lib/libiconv.2.dylib "@loader_path/libiconv.2.dylib" $@/Contents/MacOS/libunistring.2.dylib
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
%.om: %.mm $(wildcard *.h)
	$(CC) $(CFLAGS) -Wno-unused-parameter -g -c $< $(foreach inc,$(INC),-I$(inc)) -o $@

clean:
	rm -rf *.o *.om autodj AutoDJ AutoDJ.app

install:
	rsync -ax $(TARGET) $(INSTALL_DIR)/
