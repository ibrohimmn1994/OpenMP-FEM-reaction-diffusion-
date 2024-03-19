
CFLAGS= -fopenmp -Wall
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm -fopenmp

galsim: FEM.o graphics.o
	gcc -o FEM FEM.o graphics.o $(LDFLAGS)

galsim.o: FEM.c graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c FEM.c

graphics.o: graphics.c graphics.h
	gcc $(CFLAGS) $(INCLUDES) -c graphics.c

clean:
	rm -f ./FEM *.o
