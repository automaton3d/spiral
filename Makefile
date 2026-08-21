CC=gcc
CFLAGS=-std=c99 -O2 -Wall -Wextra

# SDL3 flags (adjust if SDL3 is installed elsewhere)
SDL_CFLAGS=$(shell pkg-config --cflags sdl3 2>/dev/null)
SDL_LIBS=$(shell pkg-config --libs sdl3 2>/dev/null)
ifeq ($(strip $(SDL_LIBS)),)
SDL_CFLAGS=-I/usr/include/SDL3
SDL_LIBS=-lSDL3
endif

all: spiral_axis test_axis spiral_3d

spiral_axis: spiral_axis.c
	$(CC) $(CFLAGS) -o spiral_axis spiral_axis.c

test_axis: test_axis.c spiral.c spiral.h
	$(CC) $(CFLAGS) -DNO_SDL -o test_axis test_axis.c spiral.c -lm

# Interactive 3D viewer.
# The CA core (spiral.c) is compiled headless (-DNO_SDL) so its own
# main()/render_frame are excluded, then linked with spiral_3d.c (viewer).
spiral_3d: spiral_3d.o spiral_nosdl.o
	$(CC) $(CFLAGS) -o spiral_3d spiral_3d.o spiral_nosdl.o $(SDL_LIBS) -lm

spiral_3d.o: spiral_3d.c spiral.h
	$(CC) $(CFLAGS) $(SDL_CFLAGS) -c spiral_3d.c

spiral_nosdl.o: spiral.c spiral.h
	$(CC) $(CFLAGS) -DNO_SDL -c -o spiral_nosdl.o spiral.c

clean:
	rm -f spiral_axis test_axis spiral_3d spiral_axis.exe test_axis.exe spiral_3d.exe *.obj *.o

run: spiral_axis
	./spiral_axis $(AX) $(AY) $(AZ)

run_test: test_axis
	./test_axis

run_3d: spiral_3d
	./spiral_3d $(AX) $(AY) $(AZ)