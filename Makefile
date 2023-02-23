PREFIX = /usr
CC = $(PREFIX)/bin/g++
O = 2
MARCH = nocona
CFLAGS = -O$(O) -march=$(MARCH) -mfpmath=sse -ftree-vectorize -fno-strict-aliasing -flto -Wall -Wextra -Wno-narrowing -Wno-deprecated-declarations
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
MAC_CODESIGN = $(shell whoami)

ifeq ($(UNAME), Darwin)
O = 3
MARCH = native
PREFIX = /opt/local
LIBDIR = $(PREFIX)/lib
LIB = libSDLmain.a libSDL.a libavdevice.a libavformat.a libavfilter.a libavcodec.a libavresample.a libswscale.a libavutil.a libswresample.a libbluray.a libfreetype.a libbrotlidec-static.a libbrotlicommon-static.a libpng.a libxml2.a libicuuc.a libicudata.a libiconv.a libfaac.a libfdk-aac.a libmp3lame.a libopus.a libschroedinger-1.0.a libtheoradec.a libtheora.a libvorbisenc.a libvorbis.a libvpx.a libx264.a libx265.a libxvidcore.a libgmp.a libspeex.a libmodplug.dylib libflac.a libvorbisfile.a libvorbis.a libogg.a libXrandr.a libXrender.a libXext.a libX11.a libxcb.a libXdmcp.a libXau.a libbz2.a liblzma.a libz.a libopenjp2.dylib libsoxr.dylib liborc-0.4.0.dylib libgnutls.dylib
DYLIB = libgcc/libstdc++.6.dylib libgcc/libgcc_s.1.dylib libmodplug.1.dylib libopenjp2.7.dylib libsoxr.0.dylib liborc-0.4.0.dylib libgnutls.30.dylib libsnowleopardfixes.dylib libp11-kit.0.dylib libz.1.dylib libidn.11.dylib libunistring.2.dylib libtasn1.6.dylib libhogweed.4.dylib libnettle.6.dylib libgmp.10.dylib libintl.8.dylib libffi.6.dylib libgmp.10.dylib libiconv.2.dylib
FW = Cocoa VideoDecodeAcceleration OpenGL CoreVideo CoreFoundation AudioUnit AudioToolbox IOKit Carbon
LIBRARIES += -lc++
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
CFLAGS += -DWAVEOUT -Wno-unused-parameter -mwindows -static-libgcc -static-libstdc++ -Wl,--large-address-aware -Wl,-Bstatic
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
	$(foreach lib,$(DYLIB),cp -p "$(PREFIX)/lib/$(lib)" "$@/Contents/MacOS/"; chmod 644 "$@/Contents/MacOS/$(notdir $(lib))";)
	$(foreach lib,$(DYLIB),$(PREFIX)/bin/install_name_tool -change "$(PREFIX)/lib/$(lib)" "@loader_path/$(notdir $(lib))" "$@/Contents/MacOS/AutoDJ";)
	$(foreach lib1,$(DYLIB),$(foreach lib2,$(DYLIB),$(PREFIX)/bin/install_name_tool -change "$(PREFIX)/lib/$(lib1)" "@loader_path/$(notdir $(lib1))" "$@/Contents/MacOS/$(notdir $(lib2))";))
	-codesign -s "$(MAC_CODESIGN)" $@
	-if [ -L "$@/Contents/CodeResources" ]; then rm "$@/Contents/CodeResources"; rsync -ax "$@/Contents/_CodeSignature/CodeResources" "$@/Contents/"; fi

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
