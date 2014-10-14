
DEBUG=2
CC=mpicc

#CFLAGS=-Wall -g
CFLAGS=-Wall -O -DDEBUG=$(DEBUG)

#linker flags for libraries
LDFLAGS= -lprand
PROGS=parallel_bucketsort


all: $(PROGS)

parallel_bucketsort:	parallel_bucketsort.o
	$(CC) $(CFLAGS) -o $@ parallel_bucketsort.o $(LDFLAGS)

clean:
	/bin/rm --force *.o a.out $(PROGS)
