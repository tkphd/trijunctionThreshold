# Makefile
# GNU Makefile for Voronoi tessellation and sparse phase-field grain growth
# Questions/comments to kellet@rpi.edu (Trevor Keller)

# includes
incdir = ../mmsp/include
algodir = ../mmsp/algorithms

# compilers/flags
compiler = g++ -O2 -g -Wall -ansi -pedantic
flags = -I$(incdir) -I$(algodir) -I$(algodir)/topology

mmsp2graph: mmsp2graph.cpp distance.cpp
	$(compiler) $(flags) $< -o $@ -lz -fopenmp

extract_volume: extract_volume.cpp
	$(compiler) $(flags) $< -o $@ -lz

clean:
	rm -rf extract_volume mmsp2graph
