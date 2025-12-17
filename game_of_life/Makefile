CFLAGS= -O3  -funroll-loops -march=native -fopenmp -std=c99
GFLAGS= -g
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm -fopenmp
C99FLAG = -std=c99

make: Game_of_Life.o graphics.o
	mpicc Game_of_Life.o graphics.o -o gol $(LDFLAGS)

Game_of_Life.o: Game_of_Life.c
	mpicc $(CFLAGS) $(INCLUDES) -c Game_of_Life_old.c

graphics.o: graphics.c
	gcc $(CFLAGS) $(INCLUDES) -c graphics.c -o graphics.o
clean: 
	rm -f ./gol *.o

test: test.o
	mpicc test.o -o gol $(LDFLAGS)

test.o:
	mpicc $(CFLAGS) $(INCLUDES) -c test2.c