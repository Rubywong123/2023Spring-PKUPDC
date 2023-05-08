CC = nvcc
CFLAGS = -O3 -g
LFLAGS = -lm
fd:	fd.o
		$(CC) -o fd fd.o $(LFLAGS)

fd.o:	fd.c
		$(CC) $(CFLAGS) -c fd.cu


clean:
	       rm -f *~ *.o fd

