
##### Configurable options:

## Compiler:
CC=gcc
#CC=cc

## Compiler flags:

# GCC:  (also -march=pentium etc, for machine-dependent optimizing)
CFLAGS=-Wall -O3 -fPIC -fomit-frame-pointer -funroll-loops

# GCC w/ debugging:
#CFLAGS=-Wall -g -DINLINE=

# Compaq C / Digital C:
#CFLAGS=-arch=host -fast

# SunOS:
#CFLAGS=-fast

## Program options:

# Enable long options for cl (eg. "cl --help"), comment out to disable.
# Requires getopt_long(3)  (a GNU extension)
LONGOPTS = -DENABLE_LONG_OPTIONS

LDFLAGS = -shared
##### End of configurable options


all: cl


testcases: testcases.o cliquer.o graph.o reorder.o
	$(CC) $(LDFLAGS) -o $@ testcases.o cliquer.o graph.o reorder.o

cl: cl.o cliquer.o graph.o reorder.o
	$(CC) -o $@ cl.o cliquer.o graph.o reorder.o

cl.o: cl.c
	$(CC) $(CFLAGS) $(LONGOPTS) -c $<

cliquer.o: cliquer.c
	$(CC) $(CFLAGS) $(LONGOPTS) -c $<

reorder.o: reorder.c
	$(CC) $(CFLAGS) $(LONGOPTS) -c $<

graph.o: graph.c
	$(CC) $(CFLAGS) $(LONGOPTS) -c $<

clean:
	rm -f *.o *~ cl testcases

backup:
	mkdir "`date "+backup-%Y-%m-%d-%H-%M"`" 2>/dev/null || true
	cp * "`date "+backup-%Y-%m-%d-%H-%M"`"  2>/dev/null || true

test: testcases
	./testcases

libcliquer.so: cl.o cliquer.o graph.o reorder.o
	$(CC) $(LDFLAGS) -o $@ cl.o cliquer.o graph.o reorder.o

libcliquer.a: cl.o cliquer.o graph.o reorder.o
	$(CC) -o $@ cl.o cliquer.o graph.o reorder.o


install: libcliquer.so libcliquer.a
	mv libcliquer.* lib/
