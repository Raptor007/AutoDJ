#include <Foundation/NSObject.h>
#include <Foundation/NSString.h>
#include <Foundation/NSArray.h>
#include <AppKit/NSApplication.h>
#include <SDL/SDL.h>
#include <string>


@interface SDLMain : NSObject <NSApplicationDelegate>
@end

@interface FileDrop : SDLMain
- (id)init;
- (BOOL)application:(NSApplication *)theApplication openFile:(NSString *)filename;
@end


@implementation FileDrop : SDLMain

- (id)init
{
	self = [super init];
	return self;
}

- (BOOL)application:(NSApplication *)theApplication openFile:(NSString *)filename
{
	SDL_Event event = {0};
	event.type = SDL_USEREVENT;
	event.user.data1 = new std::string([filename UTF8String]);
	SDL_PushEvent( &event );
	return TRUE;
}

@end


void FileDropEnable( void )
{
	[[NSApp delegate] release];
	[NSApp setDelegate:[[FileDrop alloc] init]];
}
