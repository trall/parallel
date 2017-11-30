all: mpiRandomWalk.o run.c
	mpicc -std=gnu99 run.c mpiRandomWalk.o -o run -lm

mpiRandomWalk.o: mpiRandomWalk.h mpiRandomWalk.c
	mpicc -std=gnu99 mpiRandomWalk.c -c -lm
