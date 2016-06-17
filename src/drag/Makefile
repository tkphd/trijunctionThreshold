# Makefile
# GNU makefile for trijunction drag grain growth code
# Questions/comments to trevor.keller@gmail.com (Trevor Keller)

# includes
incdir = ../mmsp/include

# compilers/flags
compiler = g++ -O3 -Wall
flags = -I$(incdir)
links = -lz -lpng -lgsl -lgslcblas

mobility2png: mobility2png.cpp
	$(compiler) $(flags) $< -o $@ $(links)

clean:
	rm -rf mobility2png
