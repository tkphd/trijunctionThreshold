// tricrystal.hpp
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
#include"tricrystal.hpp"

bool isLittleEndian() {
	short int number = 0x1;
	char *numPtr = (char*)&number;
	return (numPtr[0] == 1);
}

namespace MMSP {

void generate(int dim, char* filename)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	if (rank==0) {
		std::ofstream vfile("v.log",std::ofstream::out);
		vfile.close();
	}


	if (dim == 2)	{
		const int Lx = 1024;
		const int Ly = 128;

		grid<2,sparse<phi_type> > initGrid(0, 0,Lx, 0,Ly);

		for (int d=0; d<dim; d++) {
			if (x0(initGrid, d) == g0(initGrid,d))
				b0(initGrid,d) = Dirichlet;
			else if (x1(initGrid,d) == g1(initGrid,d))
				b1(initGrid,d) = Dirichlet;
		}

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			if (x[0]>Ly && x[1]>0.25*Ly && x[1]<0.75*Ly)
				set(initGrid(n), 0) = 1.0;
			else if (x[1]<Ly/2)
				set(initGrid(n), 1) = 1.0;
			else
				set(initGrid(n), 2) = 1.0;
		}

		output(initGrid, filename);
	}

	if (dim == 3)	{
		const int Lx = 512;
		const int Ly = 256;
		const int Lz = 64;

		grid<3,sparse<phi_type> > initGrid(0, 0,Lx, 0,Ly, 0,Lz);

		for (int d=0; d<dim; d++) {
			if (x0(initGrid, d) == g0(initGrid,d))
				b0(initGrid,d) = Dirichlet;
			else if (x1(initGrid,d) == g1(initGrid,d))
				b1(initGrid,d) = Dirichlet;
		}

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			if (x[0]>Ly && x[1]>0.25*Ly && x[1]<0.75*Ly)
				set(initGrid(n), 0) = 1.0;
			else if (x[1]<Ly/2)
				set(initGrid(n), 1) = 1.0;
			else
				set(initGrid(n), 2) = 1.0;
		}

		output(initGrid, filename);
	}
}

template <int dim> void update(grid<dim, sparse<phi_type> >& oldGrid, int steps)
{
	int rank=0;
	#ifdef MPI_VERSION
 	rank=MPI::COMM_WORLD.Get_rank();
	#endif
	const phi_type dt = 0.01;
	const phi_type width = 10.0;
	const phi_type epsilon = 1.0e-8;
	const phi_type mu_hi = 1.00;
	const phi_type mu_lo = 0.01;

	for (int step = 0; step < steps; step++) {
		if (rank==0) print_progress(step, steps);
		// newGrid grid must be overwritten each time
		ghostswap(oldGrid);
		grid<dim, sparse<phi_type> > newGrid(oldGrid);

		for (int d=0; d<dim; d++) {
			if (x0(oldGrid, d) == g0(oldGrid,d)) {
				b0(oldGrid,d) = Dirichlet;
				b0(newGrid,d) = Dirichlet;
			} else if (x1(oldGrid,d) == g1(oldGrid,d)) {
				b1(oldGrid,d) = Dirichlet;
				b1(newGrid,d) = Dirichlet;
			}
		}

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
					int pindex = index(oldGrid(x), h);
					set(s, pindex) = 1;
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
					// Update dFdp_h
					set(dFdp, hindex) += 0.5 * eps * eps * lap[jindex] + w * oldGrid(i)[jindex];
					// Update dFdp_j, so the inner loop can be over j>h instead of j≠h
					set(dFdp, jindex) += 0.5 * eps * eps * lap[hindex] + w * oldGrid(i)[hindex];
				  }
				}

				// compute time derivatives
				sparse<phi_type> dpdt;
				phi_type mag_phi = oldGrid(x).getMagPhi();
				phi_type mu=1.0;
				if (dim==2) mu = mobility<2>(mu_lo, mu_hi, mag_phi);
				else if (dim==3) mu = mobility<3>(mu_lo, mu_hi, mag_phi);
				for (int h = 0; h < length(s); h++) {
					int hindex = index(s, h);
					for (int j = h + 1; j < length(s); j++) {
						int jindex = index(s, j);
						set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
						set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
					}
				}

				// compute update values
				phi_type sum = 0.0;
				for (int h = 0; h < length(s); h++) {
					int pindex = index(s, h);
					phi_type value = oldGrid(i)[pindex] + dt * (2.0 / S) * dpdt[pindex]; // Extraneous factor of 2?
					if (value > 1.0) value = 1.0;
					if (value < 0.0) value = 0.0;
					if (value > epsilon) set(newGrid(i), pindex) = value;
					sum += newGrid(i)[pindex];
				}

				// project onto Gibbs simplex (enforce Σφ=1)
				phi_type rsum = 0.0;
				if (fabs(sum) > 0.0) rsum = 1.0 / sum;
				for (int h = 0; h < length(newGrid(i)); h++) {
					int pindex = index(newGrid(i), h);
					set(newGrid(i), pindex) *= rsum;
				}
			}

		} // Loop over nodes(oldGrid)

		if (step % 10 == 0) {
			// Scan along just above the mid-line for the grain boundary.
			// When found, determine its angle.
			vector<int> x(dim, 0);
			int delta = 3;
			phi_type contour = 0.505;

			x[0] = x0(newGrid,0);
			x[1] = (g1(newGrid,1) - g0(newGrid,1))/2;
			while (x[0]<x1(newGrid) && x[1]>=y0(newGrid) && x[1]<y1(newGrid) && newGrid(x)[2]>contour)
				x[0]++;
			if (x[0] == x0(newGrid))
				x[0] = g0(newGrid,0);
			int v0 = x[0];
			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Allreduce(&x[0], &v0, 1, MPI_INT, MPI_MAX);
			#endif

			x[1] += delta;
			while (x[0]>= x0(newGrid) && x[0]<x1(newGrid) && x[1]>=y0(newGrid) && x[1]<y1(newGrid) && newGrid(x)[2]>contour)
				x[0]++;
			if (x[0] == x0(newGrid))
				x[0] = g0(newGrid,0);
			int v1 = x[0];
			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Allreduce(&x[0], &v1, 1, MPI_INT, MPI_MAX);
			#endif

			x[1] += delta;
			while (x[0]>= x0(newGrid) && x[0]<x1(newGrid) && x[1]>=y0(newGrid) && x[1]<y1(newGrid) && newGrid(x)[2]>contour)
				x[0]++;
			if (x[0] == x0(newGrid))
				x[0] = g0(newGrid,0);
			int v2 = x[0];
			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Allreduce(&x[0], &v2, 1, MPI_INT, MPI_MAX);
			#endif

			double theta = 180.0/M_PI * std::atan2(dx(newGrid,1)*delta, dx(newGrid,0)*(v2-v1));
			std::ofstream vfile;
			if (rank==0) {
				vfile.open("v.log",std::ofstream::out | std::ofstream::app);
				vfile << dx(newGrid,0)*v0 << '\t' << dx(newGrid,0)*v1 << '\t' << dx(newGrid,0)*v2 << '\t' << theta << '\n';
				vfile.close();
			}
		}

		swap(oldGrid, newGrid);
	} // Loop over steps
	ghostswap(oldGrid);
}

template <class T> std::ostream& operator<<(std::ostream& o, sparse<T>& s) {
	o<<"	Index  Value\n";
	for (int i=0; i<length(s); ++i) {
		int pindex = index(s, i);
		o<<"	"<<std::setw(5)<<std::right<<pindex<<"  "<<s[pindex]<<'\n';
	}
	return o;
}

} // namespace MMSP

#endif

#include"MMSP.main.hpp"

