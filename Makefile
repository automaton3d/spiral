CC=gcc
CFLAGS=-std=c99 -O2 -Wall -Wextra

all: spiral_axis test_axis

spiral_axis: spiral_axis.c
	$(CC) $(CFLAGS) -o spiral_axis spiral_axis.c

test_axis: test_axis.c spiral.c spiral.h
	$(CC) $(CFLAGS) -DNO_SDL -o test_axis test_axis.c spiral.c -lm

clean:
	rm -f spiral_axis test_axis spiral_axis.exe test_axis.exe *.obj

run: spiral_axis
	./spiral_axis $(AX) $(AY) $(AZ)

run_test: test_axis
	./test_axis
