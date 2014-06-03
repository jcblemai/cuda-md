#
# Choose optimization suitable for your architecture by modifying the 
# OPT settings: 
OPT = -O2 -funit-at-a-time -ffast-math -finline-functions

#OPT = -O2 -funit-at-a-time -ffast-math -march=athlon-xp -finline-functions
# -malign-double

# Compiler and linker to use: 
CC=nvcc
LD=nvcc



.cc.o:
	$(CC) $(OPT) -W -ffast-math -Wall -fno-rtti -fno-exceptions -c $<

all: fluid fluid-cuda

run: all
	/usr/bin/time ./fluid | gnuplot
	#/usr/bin/time ./fluid

fluid: fluid.o alloc.o
	$(LD) fluid.o alloc.o -o fluid -lm

clean:
	-rm *.o

dist:
	mkdir ljfluid && \
		cp alloc.cc alloc.h fluid.cc Makefile ljfluid/ && \
		tar chfz fluid.tgz ljfluid/ && \
		rm -r ljfluid/

alloc.o: alloc.cc alloc.h
fluid.o: fluid.cc alloc.h

#CUDA part

CUDA_INC_PATH=/usr/local/cuda/include
CUDA_LIB_PATH=/usr/local/cuda/lib64


NVCC=nvcc
NVCFLAGS=-c -g -O3 -I$(CUDA_INC_PATH)  
LDFLAGS= -L$(CUDA_LIB_PATH) -lcudart -lcuda

fluid-cuda :
	$(NVCC) $(NVCFLAGS) fluid.cu -o fluid-cuda.o
	$(NVCC) $(LDFLAGS) fluid-cuda.o alloc.o -o fluid-cuda




