# Makefile
# GNU makefile for grid conversion programs using MMSP
# Questions/comments to gruberja@gmail.com (Jason Gruber)

# compilers/flags
compiler = g++
flags = -O2 -Wall -I ../include/

topo = ../algorithms/topology

# conversion programs
utility : mmsp2vti mmsp2pvd mmsp2png mmsp2tsv mmsp2xyz \
          mmsp2topo wrongendian

#
#    Utilities for standard visualization formats
#

# convert MMSP grid file to PNG image file type
mmsp2png : mmsp2png.cpp
	$(compiler) $(flags) $< -o $@ -lz -lpng

# convert MMSP grid file to ParaView Data file type
mmsp2pvd : mmsp2pvd.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to tab-delimited ASCII (TSV) file type
mmsp2tsv : mmsp2tsv.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to VTK Image file type
mmsp2vti : mmsp2vti.cpp
	$(compiler) $(flags) $< -o $@ -lz

# convert MMSP grid file to XYZ point cloud file type
mmsp2xyz : mmsp2xyz.cpp
	$(compiler) $(flags) $< -o $@ -lz

#
#    Utilities for niche applications
#

# extract grain boundary topology from MMSP grid file
mmsp2topo : mmsp2topo.cpp
	$(compiler) $(flags) -I $(topo) $< -o $@ -lz

# convert MMSP grid file with Little Endian data to Big Endian grid data, and vice versa
# (migrates data between machines with opposite endianness: if you're not sure, you don't need it)
wrongendian : wrongendian.cpp
	$(compiler) $(flags) $< -o $@ -lz -pthread

#
#    Cleanup
#

clean :
	rm -f mmsp2png mmsp2pvd mmsp2tsv mmsp2vti \
	      mmsp2topo mmsp2xyz wrongendian
