CC      = gcc
CFLAGS  = -O2 -Wall -Wextra
LDFLAGS = -lm

.PHONY: all test test_axis clean

all: test_axis

test: test_axis
	./test_axis

test_axis: spiral.c test_axis.c spiral.h
	$(CC) -DNO_SDL $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Build the 3D viewer when SDL3 is available.
# Adjust SDL3_CFLAGS and SDL3_LIBS if pkg-config is not installed.
SDL3_CFLAGS = $(shell pkg-config --cflags sdl3 2>/dev/null)
SDL3_LIBS   = $(shell pkg-config --libs sdl3 2>/dev/null)

spiral_auto: spiral_auto.c spiral.c spiral.h
	$(CC) $(CFLAGS) $(SDL3_CFLAGS) $^ -o $@ $(SDL3_LIBS) $(LDFLAGS)

clean:
	rm -f test_axis spiral_auto *.o
