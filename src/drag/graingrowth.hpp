// graingrowth.hpp
// Anisotropic sparse phase field (sparsePF) grain growth example code
// Questions/comments to gruberja@gmail.com (Jason Gruber)

std::string PROGRAM = "graingrowth";
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

template <int dim, class T> T mobility(const T mu_lo, const T mu_hi, const T phi)
{
	// Create an annulus of low mobility near vertices

	T mean = 0.6435; // empirically determined
	T devn = 0.0125;

	/* TJ Drag: *only* on triple junctions; vertices unaffected!
	if (dim==3) {
		mean = 0.56; // empirically determined
		devn = 0.005;
	}
	*/

	// Inverted, shifted Gaussian curve
	T gauss = std::exp( -1.0*(phi - mean)*(phi - mean)/(2.0*devn*devn) );
	T mobility = mu_hi - (mu_hi - mu_lo) * gauss;

	if (mobility > mu_hi + 0.001) {
		std::cerr<<"Illegal mobility: too high! "<<mobility<<std::endl;
		std::exit(-1);
	} else if (mobility < mu_lo - 0.001) {
		std::cerr<<"Illegal mobility: too low! "<<mobility<<std::endl;
		std::exit(-1);
	}
	return mobility;
}

#endif
