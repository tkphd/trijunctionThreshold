// tricrystal.hpp
// Anisotropic sparse phase field (sparsePF) grain growth example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

std::string PROGRAM = "tricrystal";
std::string MESSAGE = "Sparse phase field (sparsePF) grain growth code with vertex drag";

typedef double phi_type;

typedef MMSP::grid<1,MMSP::sparse<phi_type> > GRID1D;
typedef MMSP::grid<2,MMSP::sparse<phi_type> > GRID2D;
typedef MMSP::grid<3,MMSP::sparse<phi_type> > GRID3D;


// Anisotropic energy and mobility functions
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ANISOTROPY
#define ANISOTROPY
#include<iostream>
#include<cmath>

phi_type energy(int i, int j)
{
	// trivial case: no boundary
	if (i==j) return 0.0;
	return 1.0;
}

phi_type mobility(int i, int j)
{
	// trivial case: no boundary
	if (i==j) return 0.0;
	return 1.0;
}

#endif
