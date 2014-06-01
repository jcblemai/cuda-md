#
# Choose optimization suitable for your architecture by modifying the 
# OPT settings: 
OPT = -O2 -funit-at-a-time -ffast-math -finline-functions

#OPT = -O2 -funit-at-a-time -ffast-math -march=athlon-xp -finline-functions
# -malign-double

# Compiler and linker to use: 
CC=g++
LD=g++

.cc.o:
	$(CC) $(OPT) -W -ffast-math -Wall -fno-rtti -fno-exceptions -c $<

all: fluid

run: all
	/usr/bin/time ./fluid | gnuplot
	#/usr/bin/time ./fluid

fluid: fluid.o alloc.o
	$(LD) fluid.o alloc.o -o fluid -lm

clean:
	-rm fluid.o alloc.o

dist:
	mkdir ljfluid && \
		cp alloc.cc alloc.h fluid.cc Makefile ljfluid/ && \
		tar chfz fluid.tgz ljfluid/ && \
		rm -r ljfluid/

alloc.o: alloc.cc alloc.h
fluid.o: fluid.cc alloc.h
