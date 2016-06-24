# Makefile
# GNU makefile for trijunction drag grain growth code
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../../mmsp/include

# compilers/flags
compiler = g++ -O3 -Wall
pcompiler = /usr/bin/mpic++.openmpi -O3 -Wall
flags = -I../ -I$(incdir) -I$(incdir)/../algorithms

# IBM compiler for AMOS
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_XL)/include
BG_LIB = -L$(BG_XL)/lib
qcompiler = $(BG_XL)/bin/mpixlcxx -O3 -qflag=w -qstrict -qmaxmem=-1

qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib
#qflags = $(BG_INC) $(BG_LIB) $(flags) $(CFLAGS) $(LDFLAGS)

# the program
vertex.out: vertex.cpp
	$(compiler) $(flags) $< -o $@ -lz

parallel: vertex.cpp
	$(pcompiler) $(flags) -include mpi.h $< -o parallel_GG.out -lz

bgq: vertex.cpp
	$(qcompiler) -DBGQ $(qflags) $< -o q_GG.out -lz

clean:
	rm -rf vertex.out parallel_GG.out q_GG.out
