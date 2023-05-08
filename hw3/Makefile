CC = mpicc
CFLAGS = -O2 -Wall -g
LFLAGS = -ljpeg -lm
sobel:	sobel.o
		$(CC) -o sobel sobel.o $(LFLAGS)

sobel.o:	sobel.c
		$(CC) $(CFLAGS) -c sobel.c


clean:
	       rm -f *~ *.o sobel out.jpg

