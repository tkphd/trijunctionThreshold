// graingrowth.hpp
// Anisotropic coarsening algorithms for 2D and 3D sparse phase field (sparsePF) methods
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE

#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#ifdef BGQ
#include<mpi.h>
#endif
#include"MMSP.hpp"
#include"graingrowth.hpp"
#include"tessellate.hpp"

bool isLittleEndian() {
	short int number = 0x1;
	char *numPtr = (char*)&number;
	return (numPtr[0] == 1);
}

namespace MMSP {

void generate(int dim, char* filename) {
	#if (!defined MPI_VERSION || MPI_VERSION<2) && (defined BGQ)
	std::cerr<<"Error: MPI-2 is required for CCNI."<<std::endl;
	exit(1);
	#endif
	static unsigned long tstart;

 	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	#ifdef DEBUG
	if (rank==0) {
		if (!isLittleEndian()) std::cout<<"Byte order is big-endian."<<std::endl;
	}
	#endif
	if (dim == 2)	{
		int edge = 8192;
		int number_of_fields = 5120000;
		grid<2,sparse<phi_type> > initGrid(0, 0, edge, 0, edge);
		if (rank==0) std::cout<<"Grid origin: ("<<g0(initGrid,0)<<','<<g0(initGrid,1)<<"),"
		                        <<" dimensions: "<<g1(initGrid,0)-g0(initGrid,0)<<" × "<<g1(initGrid,1)-g0(initGrid,1)
		                        <<" with "<<number_of_fields<<" seeds"<<std::flush;
		#ifdef MPI_VERSION
		number_of_fields /= np;
		if (rank==0) std::cout<<", "<<number_of_fields<<" per rank"<<std::flush;
		if (rank==0 && number_of_fields % np != 0)
			std::cerr<<"\nWarning: Tessellation may hang with uneven distribution of seeds per thread."<<std::endl;
		#endif
		if (rank==0) std::cout<<"."<<std::endl;

		#if (!defined MPI_VERSION) && (defined BGQ)
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		tstart = time(NULL);
		tessellate<2,phi_type>(initGrid, number_of_fields);
		if (rank==0) std::cout<<"Tessellation complete ("<<time(NULL)-tstart<<" sec)."<<std::endl;
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		tstart=time(NULL);
		output(initGrid, filename);
		if (rank==0) std::cout<<"Voronoi tessellation written to "<<filename<<" ("<<time(NULL)-tstart<<" sec)."<<std::endl;
	}
}

template <int dim> void update(grid<dim, sparse<phi_type> >& oldGrid, int steps) {
	#if (!defined MPI_VERSION) && (defined BGQ)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
 	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	const phi_type dt = 0.01;
	const phi_type width = 10.0;
	const phi_type epsilon = 1.0e-8;

	static int iterations = 1;

	for (int step = 0; step < steps; step++) {
		if (rank==0) print_progress(step, steps);
		// update grid must be overwritten each time
		ghostswap(oldGrid);
		grid<dim, sparse<phi_type> > newGrid(oldGrid);

		for (int i = 0; i < nodes(oldGrid); i++) {
			vector<int> x = position(oldGrid, i);

			// determine nonzero fields within
			// the neighborhood of this node
			// (2 adjacent voxels along each cardinal direction)
			sparse<int> s;
			for (int j = 0; j < dim; j++)
				for (int k = -1; k <= 1; k++) {
				  x[j] += k;
				  for (int h = 0; h < length(oldGrid(x)); h++) {
				    int sindex = index(oldGrid(x), h);
				    set(s, sindex) = 1;
				  }
				  x[j] -= k;
				}
			phi_type S = phi_type(length(s));

			// if only one field is nonzero,
			// then copy this node to newGrid
			if (S < 2.0) newGrid(i) = oldGrid(i);
			else {
				// compute laplacian of each field
				sparse<phi_type> lap = laplacian(oldGrid, i);

				// compute variational derivatives
				sparse<phi_type> dFdp;
				for (int h = 0; h < length(s); h++) {
				  int hindex = index(s, h);
				  for (int j = h + 1; j < length(s); j++) {
				    int jindex = index(s, j);
				    phi_type gamma = energy(hindex, jindex);
				    phi_type eps = 4.0 / acos(-1.0) * sqrt(0.5 * gamma * width);
				    phi_type w = 4.0 * gamma / width;
				    // Update dFdp_h and dFdp_j, so the inner loop can be over j>h instead of j≠h
				    set(dFdp, hindex) += 0.5 * eps * eps * lap[jindex] + w * oldGrid(i)[jindex];
				    set(dFdp, jindex) += 0.5 * eps * eps * lap[hindex] + w * oldGrid(i)[hindex];
				  }
				}

				// compute time derivatives
				sparse<phi_type> dpdt;
				for (int h = 0; h < length(s); h++) {
				  int hindex = index(s, h);
				  for (int j = h + 1; j < length(s); j++) {
				    int jindex = index(s, j);
				    phi_type mu = mobility(hindex, jindex);
				    set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
				    set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
				  }
				}

				// compute update values
				phi_type sum = 0.0;
				for (int h = 0; h < length(s); h++) {
				  int sindex = index(s, h);
				  phi_type value = oldGrid(i)[sindex] + dt * (2.0 / S) * dpdt[sindex]; // Extraneous factor of 2?
				  if (value > 1.0) value = 1.0;
				  if (value < 0.0) value = 0.0;
				  if (value > epsilon) set(newGrid(i), sindex) = value;
				  sum += newGrid(i)[sindex];
				}

				// project onto Gibbs simplex (enforce Σφ=1)
				phi_type rsum = 0.0;
				if (fabs(sum) > 0.0) rsum = 1.0 / sum;
				for (int h = 0; h < length(newGrid(i)); h++) {
				  int sindex = index(newGrid(i), h);
				  set(newGrid(i), sindex) *= rsum;
				}
			}
		} // Loop over nodes(oldGrid)
		swap(oldGrid, newGrid);
	} // Loop over steps
	ghostswap(oldGrid);
	++iterations;
}

template <class T> std::ostream& operator<<(std::ostream& o, sparse<T>& s) {
	o<<"    Index  Value\n";
	for (int i=0; i<length(s); ++i) {
		int sindex = index(s, i);
		o<<"    "<<std::setw(5)<<std::right<<sindex<<"  "<<s[sindex]<<'\n';
	}
	return o;
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

