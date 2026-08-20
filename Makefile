CC=gcc
CFLAGS=-std=c99 -O2 -Wall -Wextra

spiral_axis: spiral_axis.c
	$(CC) $(CFLAGS) -o spiral_axis spiral_axis.c

clean:
	rm -f spiral_axis spiral_axis.exe spiral_axis.obj

run: spiral_axis
	./spiral_axis $(AX) $(AY) $(AZ)
