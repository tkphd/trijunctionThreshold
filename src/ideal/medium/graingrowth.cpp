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
#if (defined CCNI) || (defined BGQ)
#include<mpi.h>
#endif
#include"MMSP.hpp"
#include"anisotropy.hpp"
#include"graingrowth.hpp"
#include"tessellate.hpp"

typedef float phi_type;

void print_progress(const int step, const int steps, const int iterations);

bool isLittleEndian() {
	short int number = 0x1;
	char *numPtr = (char*)&number;
	return (numPtr[0] == 1);
}

namespace MMSP {

void generate(int dim, char* filename) {
	#if (!defined MPI_VERSION || MPI_VERSION<2) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI-2 is required for CCNI."<<std::endl;
	exit(1);
	#endif
  static unsigned long tstart;

 	int id=0;
	#ifdef MPI_VERSION
	id = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	#ifdef DEBUG
	if (id==0) {
		if (!isLittleEndian()) std::cout<<"Byte order is big-endian."<<std::endl;
	}
	#endif
	if (dim == 2)	{
		const int edge = 16384;
		int number_of_fields = edge*10000/2048; //static_cast<int>(pow(edge/10.0,2.0)/4); //edge*10000/2048; // average grain is a disk of radius 10
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		MMSP::grid<2,MMSP::sparse<phi_type> > grid(0, 0, edge, 0, edge);
		if (id==0) std::cout<<"Grid origin: ("<<g0(grid,0)<<','<<g0(grid,1)<<"),"
				                <<" dimensions: "<<g1(grid,0)-g0(grid,0)<<" × "<<g1(grid,1)-g0(grid,1)
				                <<" with "<<number_of_fields<<" seeds"<<std::flush;
		#ifdef MPI_VERSION
		number_of_fields /= np;
		if (id==0) std::cout<<", "<<number_of_fields<<" per rank"<<std::flush;
		#endif
		if (id==0) std::cout<<"."<<std::endl;

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
	  tstart = time(NULL);
		tessellate<2,phi_type>(grid, number_of_fields);
		if (id==0) std::cout<<"Tessellation complete ("<<time(NULL)-tstart<<" sec)."<<std::endl;
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		tstart=time(NULL);
		output(grid, filename);
		if (id==0) std::cout<<"Voronoi tessellation written to "<<filename<<" ("<<time(NULL)-tstart<<" sec)."<<std::endl;
	}

	if (dim == 3)	{
		const int edge = 512;
		/*
		//int number_of_fields = static_cast<int>(pow(edge/8.0,3.0)/8); // Average grain is a sphere of radius 10 voxels
		int number_of_fields = static_cast<int>(pow(edge,3.0)/(64*pow(8,3.0)/sqrt(27))); // BCC packing with R=8 vox, Vc=(4R/sqrt(3))^3.
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		*/
		int number_of_fields = 8192;
		MMSP::grid<3,MMSP::sparse<phi_type> > grid(0,0,edge,0,edge,0,edge);
		if (id==0) std::cout<<"Grid origin: ("<<g0(grid,0)<<','<<g0(grid,1)<<','<<g0(grid,2)<<"),"
				                <<" dimensions: "<<g1(grid,0)-g0(grid,0)<<" × "<<g1(grid,1)-g0(grid,1)<<" × "<<g1(grid,2)-g0(grid,2)
				                <<" with "<<number_of_fields<<" seeds"<<std::flush;
		#ifdef MPI_VERSION
		number_of_fields /= np;
		if (id==0) std::cout<<", "<<number_of_fields<<" per rank"<<std::flush;
		#endif
		if (id==0) std::cout<<"."<<std::endl;

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
	  tstart = time(NULL);
		tessellate<3,phi_type>(grid, number_of_fields);
		//const std::string seedfile("points2.txt");
		//tessellate<3,phi_type>(grid, seedfile);
		if (id==0) std::cout<<"Tessellation complete ("<<time(NULL)-tstart<<" sec)."<<std::endl;
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		tstart=time(NULL);
		output(grid, filename);
		if (id==0) std::cout<<"Voronoi tessellation written to "<<filename<<" ("<<time(NULL)-tstart<<" sec)."<<std::endl;
	}
}

template <int dim> void update(MMSP::grid<dim, sparse<phi_type> >& grid, int steps) {
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int id=0;
	#ifdef MPI_VERSION
 	id=MPI::COMM_WORLD.Get_rank();
	#endif
	const phi_type dt = 0.01;
	const phi_type width = 10.0;
	const phi_type epsilon = 1.0e-8;

	static int iterations = 1;

	for (int step = 0; step < steps; step++) {
		if (id==0) print_progress(step, steps, iterations);
		// update grid must be overwritten each time
		ghostswap(grid);
		MMSP::grid<dim, sparse<phi_type> > update(grid);

		for (int i = 0; i < nodes(grid); i++) {
			vector<int> x = position(grid, i);

			// determine nonzero fields within
			// the neighborhood of this node
			// (2 adjacent voxels along each cardinal direction)
			sparse<int> s;
			for (int j = 0; j < dim; j++)
				for (int k = -1; k <= 1; k++) {
				  x[j] += k;
				  for (int h = 0; h < length(grid(x)); h++) {
				    int index = MMSP::index(grid(x), h);
				    set(s, index) = 1;
				  }
				  x[j] -= k;
				}
			phi_type S = phi_type(length(s));

			// if only one field is nonzero,
			// then copy this node to update
			if (S < 2.0) update(i) = grid(i);
			else {
				// compute laplacian of each field
				sparse<phi_type> lap = laplacian(grid, i);

				// compute variational derivatives
				sparse<phi_type> dFdp;
				for (int h = 0; h < length(s); h++) {
				  int hindex = MMSP::index(s, h);
				  for (int j = h + 1; j < length(s); j++) {
				    int jindex = MMSP::index(s, j);
				    phi_type gamma = energy(hindex, jindex);
				    phi_type eps = 4.0 / acos(-1.0) * sqrt(0.5 * gamma * width);
				    phi_type w = 4.0 * gamma / width;
				    // Update dFdp_h
				    set(dFdp, hindex) += 0.5 * eps * eps * lap[jindex] + w * grid(i)[jindex];
				    // Update dFdp_j, so the inner loop can be over j>h instead of j≠h
				    set(dFdp, jindex) += 0.5 * eps * eps * lap[hindex] + w * grid(i)[hindex];
				  }
				}

				// compute time derivatives
				sparse<phi_type> dpdt;
				for (int h = 0; h < length(s); h++) {
				  int hindex = MMSP::index(s, h);
				  for (int j = h + 1; j < length(s); j++) {
				    int jindex = MMSP::index(s, j);
				    phi_type mu = mobility(hindex, jindex);
				    set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
				    set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
				  }
				}

				// compute update values
				phi_type sum = 0.0;
				for (int h = 0; h < length(s); h++) {
				  int index = MMSP::index(s, h);
				  phi_type value = grid(i)[index] + dt * (2.0 / S) * dpdt[index]; // Extraneous factor of 2?
				  if (value > 1.0) value = 1.0;
				  if (value < 0.0) value = 0.0;
				  if (value > epsilon) set(update(i), index) = value;
				  sum += update(i)[index];
				}

				// project onto Gibbs simplex (enforce Σφ=1)
				phi_type rsum = 0.0;
				if (fabs(sum) > 0.0) rsum = 1.0 / sum;
				for (int h = 0; h < length(update(i)); h++) {
				  int index = MMSP::index(update(i), h);
				  set(update(i), index) *= rsum;
				}
			}
		} // Loop over nodes(grid)
		swap(grid, update);
	} // Loop over steps
	ghostswap(grid);
	++iterations;
}

template <class T> std::ostream& operator<<(std::ostream& o, sparse<T>& s) {
	o<<"    Index  Value\n";
	for (int i=0; i<length(s); ++i) {
		int index = MMSP::index(s, i);
		o<<"    "<<std::setw(5)<<std::right<<index<<"  "<<s[index]<<'\n';
	}
	return o;
}

} // namespace MMSP

void print_progress(const int step, const int steps, const int iterations) {
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step==0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
	} else if (step==steps-1) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout<<"•] "<<std::setw(2)<<std::right<<deltat/3600<<"h:"
										<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
										<<std::setw(2)<<std::right<<deltat%60<<"s"
										<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}
#endif

#include"MMSP.main.hpp"

