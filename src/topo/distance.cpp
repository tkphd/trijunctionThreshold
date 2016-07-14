// distance.cpp
// Implementations of methods to propagate vertex distances in MMSP grids
// Note: Parallelization is provided by OpenMP, not MPI.

#ifndef _DISTANCE_CPP_
#define _DISTANCE_CPP_
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <ctime>
#include <limits>
#include <cassert>
#include <omp.h>

#include "MMSP.hpp"
#include "sparse_priority_queue.h"
#include "distance.h"

template <int dim, typename T> double radius(const MMSP::vector<T>& a, const MMSP::vector<T>& b)
{
	double radius=0.0;
	for (int d=0; d<dim; d++)
		radius += (b[d]-a[d]) * (b[d]-a[d]);
	return sqrt(radius);
}


MMSP::vector<int> getPosition(const SparseDistanceVoxel& dv)
{
	MMSP::vector<int> x(3);
	for (int d=0; d<3; d++)
		x[d]=dv.getPosition(d);
	return x;
}


template<int dim, typename T, typename U>
void propagate_cost(const SparseDistanceVoxel* core_voxel,
                    const int index,
                    const MMSP::grid<dim,MMSP::sparse<U> >& grid,
                    const MMSP::grid<dim,T>& cost,
                    const std::set<int>& neighbor_ids,
                    MMSP::grid<dim,SparseDistanceVoxel>& distimage,
                    SparsePriorityQueue& queue)
{
	if (dim == 2) {
		const MMSP::vector<int> position = getPosition(*core_voxel);
		const int core_distance = distimage(position).getValue(index);
		MMSP::vector<int> x(position);
		// Loop over 9 neighboring Voxels
		for (x[1]=position[1]-1; x[1]<position[1]+2; x[1]++) {
			for (x[0]=position[0]-1; x[0]<position[0]+2; x[0]++) {
				if (x[0]==position[0] && x[1]==position[1]) continue; // know thyself
				T step=radius<2,int>(position, x);
				MMSP::vector<int> p(x);
				for (int d=0; d<dim; d++)
					MMSP::check_boundary(p[d],MMSP::x0(grid,d),MMSP::x1(grid,d),b0(grid,d),b1(grid,d));

				if (grid(p).length()<2) continue;
				// Make sure test point has at least one phase in common with source vertex
				int matching_ids = 0;
				const MMSP::sparse<U>* s = &(grid(p));
				for (int i=0; i<s->length(); i++) {
					int a=s->index(i);
					if (neighbor_ids.find(a) != neighbor_ids.end())
						++matching_ids; // match found :)
				}
				if (matching_ids>0) {
					// Enqueue valid voxel
					SparseDistanceVoxel* Voxel = &( distimage(p) );
					T distance = core_distance + step + cost(p);
					if (Voxel->getValue(index) > distance) {
						Voxel->setValue(index,distance); // This part requires mutex locks for OpenMP
						if ( !queue.in_heap( Voxel ) ) queue.push(index, Voxel); // Add new Voxel to heap
						else queue.update_position(index, Voxel); // Update position of Voxel in the heap based on new value
					}
				}
			}
		}
	} else if (dim == 3) {
		const MMSP::vector<int> position = getPosition(*core_voxel);
		const int core_distance = distimage(position).getValue(index);
		MMSP::vector<int> x(position);
		// Loop over 27 neighboring Voxels
		for (x[2]=position[2]-1; x[2]<position[2]+2; x[2]++) {
			for (x[1]=position[1]-1; x[1]<position[1]+2; x[1]++) {
				for (x[0]=position[0]-1; x[0]<position[0]+2; x[0]++) {
					if (x[0]==position[0] && x[1]==position[1] && x[2]==position[2]) continue; // know thyself
					T step=radius<3,int>(position, x);
					MMSP::vector<int> p(x);
					for (int d=0; d<dim; d++)
						MMSP::check_boundary(p[d],MMSP::x0(grid,d),MMSP::x1(grid,d),b0(grid,d),b1(grid,d));

					if (grid(p).length()<3) continue;
					// Make sure test point has at least one phase in common with source vertex
					int matching_ids = 0;
					const MMSP::sparse<U>* s = &(grid(p));
					for (int i=0; i<s->length(); i++) {
						int a=s->index(i);
						if (neighbor_ids.find(a) != neighbor_ids.end())
							++matching_ids; // match found :)
					}
					if (matching_ids>1) {
						// Enqueue valid voxel
						SparseDistanceVoxel* Voxel = &( distimage(p) );
						T distance = core_distance + step + cost(p);
						if (Voxel->getValue(index) > distance) {
							Voxel->setValue(index,distance);
							if ( !queue.in_heap( Voxel ) ) queue.push(index, Voxel); // Add new Voxel to heap
							else queue.update_position(index, Voxel); // Update position of Voxel in the heap based on new value
						}
					}
				}
			}
		}
	}
} // propagate_cost


template<int dim, typename T>
bool isLocalMin(const MMSP::grid<dim,MMSP::sparse<T> >& grid, const MMSP::vector<int> position)
{
	MMSP::vector<int> x(position);
	const double center=grid(x).getMagPhi();

	if (grid(x).length() <= dim) return false;

	if (dim==2) {
		for (x[1]=position[1]-1; x[1]<=position[1]+1; x[1]++) {
			for (x[0]=position[0]-1; x[0]<=position[0]+1; x[0]++) {
				if (x[0]==position[0] && x[1]==position[1]) continue; // know thyself
				MMSP::vector<int> p(x);
				for (int d=0; d<dim; d++)
					MMSP::check_boundary(p[d],MMSP::x0(grid,d),MMSP::x1(grid,d),MMSP::b0(grid,d),MMSP::b1(grid,d));
				if (center>grid(p).getMagPhi()) return false;
			}
		}
	} else if (dim==3) {
		for (x[2]=position[2]-1; x[2]<=position[2]+1; x[2]++) {
			for (x[1]=position[1]-1; x[1]<=position[1]+1; x[1]++) {
				for (x[0]=position[0]-1; x[0]<=position[0]+1; x[0]++) {
					if (x[0]==position[0] && x[1]==position[1] && x[2]==position[2]) continue; // know thyself
					MMSP::vector<int> p(x);
					for (int d=0; d<dim; d++)
						MMSP::check_boundary(p[d],MMSP::x0(grid,d),MMSP::x1(grid,d),MMSP::b0(grid,d),MMSP::b1(grid,d));
					if (center>grid(p).getMagPhi()) return false;
				}
			}
		}
	}
	return true;
} // isLocalMin


template<int dim, typename T>
void locate_vertices(const MMSP::grid<dim,MMSP::sparse<T> >& grid, std::map<int,double>& grainsizes, std::vector<MMSP::vector<int> >& global_vertices)
{
	global_vertices.clear();

	// Associate grain IDs with each vertex. If duplicate sets exist, retain the best candidate.
	std::map<std::set<int>,MMSP::vector<int> > phase_ids; // map set of phases to vertex position
	MMSP::vector<int> vertex(3,0);
	for (int n=0; n<nodes(grid); n++) {
		MMSP::vector<int> x=MMSP::position(grid,n);
		for (int i=0; i<grid(x).length(); i++)
			grainsizes[grid(x).index(i)]+=grid(x).value(i);
		if (isLocalMin(grid,x)) {
			for (int d=0; d<dim; d++)
				vertex[d]=x[d];
			// Associate a set of neighbor-ids with this queue
			std::set<int> local_ids;
			for (int j=0; j<grid(vertex).length(); j++)
				local_ids.insert( grid(vertex).index(j) );
			std::map<std::set<int>,MMSP::vector<int> >::iterator insert_set = phase_ids.find(local_ids);
			if (insert_set==phase_ids.end()) {
				// Vertex is unique. Rejoice!
				phase_ids[local_ids]=vertex; // map set of phases to this vertex for comparison
				global_vertices.push_back(vertex); // store vertex in global collection
			}
		}
	}
} // locate_vertices


template<int dim, typename T>
void determine_weights(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                       const std::vector<MMSP::vector<int> >& global_vertices,
                       MMSP::grid<dim,double>& cost)
{
	double min=M_PI; // pi will not naturally occur: if you see it, something's wrong
	double max=-M_PI;
	double bulk_penalty = 1000.0;
	double bulk_threshold = 0.7125; //1.01 / sqrt(2.0); // must be above 0.70711, the magnitude of an edge

	double epsilonsq = 4.5e-5;

	if (dim==2) {
		// Cost of moving through the grain boundary equals the curvature of |phi| at
		// the point. Mean curvature 2H = -div.normal=-div.(grad |phi|/|grad|phi||).
		// If 2H>0, take cost as |phi| instead.

		MMSP::grid<2,double> magnitude(0,MMSP::g0(grid,0),MMSP::g1(grid,0),MMSP::g0(grid,1),MMSP::g1(grid,1));
		for (int d=0; d<dim; d++)
			dx(magnitude,d) = dx(grid,d);

		#pragma omp parallel for
		for (int i=0; i<MMSP::nodes(grid); i++) {
			MMSP::vector<int> x = MMSP::position(grid,i);
			magnitude(x) = grid(i).getMagPhi();
		}

		#pragma omp parallel
		{
			double myMin = min;
			double myMax = max;

			#pragma omp for nowait
			for (int i=0; i<MMSP::nodes(magnitude); i++) {
				MMSP::vector<int> x = MMSP::position(magnitude,i);
				double H = 0.0;
				double lap = 0.0;

				// Get central values
				MMSP::vector<T> grad = MMSP::gradient(magnitude, x);
				double maggradsq = grad*grad;
				const T& pc = magnitude(x);
				const T  Mc = (maggradsq>epsilonsq) ? 1.0 / sqrt(maggradsq) : 0.0;

				for (int d=0; d<dim; d++) {
					// Central second-order difference
					// Get low values
					x[d] -= 1;
					grad = MMSP::gradient(magnitude, x);
					maggradsq = grad*grad;
					const T& pl = magnitude(x);
					const T  Ml = (maggradsq>epsilonsq) ? 1.0 / sqrt(maggradsq) : 0.0;

					// Get high values
					x[d] += 2;
					grad = MMSP::gradient(magnitude, x);
					maggradsq = grad*grad;
					const T& ph = magnitude(x);
					const T  Mh = (maggradsq>epsilonsq) ? 1.0 / sqrt(maggradsq) : 0.0;

					// Return to center
					x[d] -= 1;

					// Put 'em all together (note -=, since it's the negative curvature)
					double weight = 1.0/(dx(grid,d) * dx(grid,d));
					H -= 0.5*weight*( (Mh+Mc)*(ph-pc) - (Mc+Ml)*(pc-pl) );
					lap -= weight*(ph - 2.0*pc + pl);
				}

				if (H>0.0) {
					if (lap<0.0)
						H = lap;
					else
						H = magnitude(i);
				}
				if (H < myMin)
					myMin = H;
				if (H > myMax)
					myMax = H;
				cost(x) = H;
			}

			#pragma omp critical
			{
				if (myMin < min)
					min = myMin;
				if (myMax > max)
					max = myMax;
			}

		}

		#pragma omp parallel for
		for (int i=0; i<MMSP::nodes(cost); i++) {
			// Set high cost in bulk, and bump grain boundary costs positive
			MMSP::vector<int> x = position(cost, i);
			if (magnitude(x) > bulk_threshold)
				cost(i) = bulk_penalty * max;
			else
				cost(i) = cost(i) - min;
		}
	} else if (dim==3) {
		// Cost of moving through the grain boundary equals the curvature of |phi| at
		// the point. Mean curvature 2H = -div.normal=-div.(grad |phi|/|grad|phi||).
		// If 2H>0, take cost as |phi| instead.

		MMSP::grid<3,double> magnitude(0,MMSP::g0(grid,0),MMSP::g1(grid,0),MMSP::g0(grid,1),MMSP::g1(grid,1),MMSP::g0(grid,2),MMSP::g1(grid,2));
		for (int d=0; d<dim; d++)
			dx(magnitude,d) = dx(grid,d);

		#pragma omp parallel for
		for (int i=0; i<MMSP::nodes(grid); i++) {
			MMSP::vector<int> x = MMSP::position(grid,i);
			magnitude(x) = grid(i).getMagPhi();
		}

		#pragma omp parallel
		{
			double myMin = min;
			double myMax = max;

			#pragma omp for nowait
			for (int i=0; i<MMSP::nodes(magnitude); i++) {
				MMSP::vector<int> x = MMSP::position(magnitude,i);
				double H = 0.0;
				double lap = 0.0;

				// Get central values
				MMSP::vector<T> grad = MMSP::gradient(magnitude, x);
				double maggradsq = grad*grad;
				const T& pc = magnitude(x);
				const T  Mc = (maggradsq>epsilonsq) ? 1.0 / sqrt(maggradsq) : 0.0;

				for (int d=0; d<dim; d++) {
					// Central second-order difference
					// Get low values
					x[d] -= 1;
					grad = MMSP::gradient(magnitude, x);
					maggradsq = grad*grad;
					const T& pl = magnitude(x);
					const T  Ml = (maggradsq>epsilonsq) ? 1.0 / sqrt(maggradsq) : 0.0;

					// Get high values
					x[d] += 2;
					grad = MMSP::gradient(magnitude, x);
					maggradsq = grad*grad;
					const T& ph = magnitude(x);
					const T  Mh = (maggradsq>epsilonsq) ? 1.0 / sqrt(maggradsq) : 0.0;

					// Return to center
					x[d] -= 1;

					// Put 'em all together (note -=, since it's the negative curvature)
					double weight = 1.0/(dx(grid,d) * dx(grid,d));
					H -= 0.5*weight*( (Mh+Mc)*(ph-pc) - (Mc+Ml)*(pc-pl) );
					lap -= weight*(ph - 2.0*pc + pl);
				}

				if (H>0.0) {
					if (lap<0.0)
						H = lap;
					else
						H = magnitude(i);
				}
				if (H < myMin)
					myMin = H;
				if (H > myMax)
					myMax = H;
				cost(x) = H;
			}

			#pragma omp critical
			{
				if (myMin < min)
					min = myMin;
				if (myMax > max)
					max = myMax;
			}

		}

		#pragma omp parallel for
		for (int i=0; i<MMSP::nodes(cost); i++) {
			// Set high cost in bulk, and bump grain boundary costs positive
			MMSP::vector<int> x = position(cost, i);
			if (magnitude(x) > bulk_threshold)
				cost(i) = bulk_penalty * max;
			else
				cost(i) = cost(i) - min;
		}
	}
} // determine_weights


template<int dim, typename T, typename U>
void approximate_cost(const MMSP::grid<dim,MMSP::sparse<U> >& phi,
                      const MMSP::grid<dim,T>& cost,
                      std::vector<MMSP::vector<int> >& global_vertices,
                      MMSP::grid<dim,SparseDistanceVoxel>& distance_grid,
                      MMSP::grid<dim,MMSP::sparse<T> >& grid)
{
	// Implements a fast marching algorithm to generate the distance map
	// Based on code written by Prof. Barbara Cutler, RPI Comp. Sci. Dept., for CSCI-1200.
	if (global_vertices.size()==0) return;

	MMSP::vector<int> x;

	// Initialize distance_grid; explicitly scope variables for OpenMP
	for (int i=0; i<MMSP::nodes(distance_grid); i++) {
		x=MMSP::position(grid,i);
		for (int d=0; d<dim; d++)
			distance_grid(x).setPosition(d,x[d]);
	}

	// Associate grain IDs with each vertex.
	std::vector<std::set<int> > phase_ids;
	phase_ids.resize(global_vertices.size());

	for (int i=0; i<int(global_vertices.size()); i++) {
		MMSP::sparse<U>* s = &(phi(global_vertices[i]));
		// Associate a set of neighbor-ids with this queue
		for (int j=0; j<s->length(); j++)
			phase_ids[i].insert( s->index(j) );
	}

	// create the voxel Heap
	std::vector<SparsePriorityQueue> queues;
	if (1) {
		SparsePriorityQueue queue;
		while (queues.size()<global_vertices.size())
			queues.push_back(queue);
	}

	// Propagating the distance to each vertex in parallel using OpenMP.
	// Before execution,      export OMP_NUM_THREADS=xxx
	// If adjacent vertices are being marched, this loop WILL CRASH
	// when two threads modify the same SparseDistanceVoxel, since
	// adding an index triggers regeneration of the sparse guts.
	// Needs mutual exclusion locks.

	// Each priority queue should be independent: operating on index i only.

	for (int i=0; i<int(global_vertices.size()); i++) {
		x = global_vertices[i];

		// Enqueue this node's vertices (local)
		SparseDistanceVoxel* p = &( distance_grid(x) );
		p->setValue(i, 0.0);
		propagate_cost(p, i, phi, cost, phase_ids[i], distance_grid, queues[i]);
	}

	for (int i=0; i<int(global_vertices.size()); i++) {
		x = global_vertices[i];

		// Fast-march grid points associated with this vertex
		while ( !(queues[i].empty()) ) {
			const SparseDistanceVoxel* p = queues[i].top();
			queues[i].pop(i);
			propagate_cost(p, i, phi, cost, phase_ids[i], distance_grid, queues[i]);
		}

		assert(queues[i].empty());

	}

	// Copy result from distance_grid to phase-field grid
	for (int i=0; i<MMSP::nodes(grid); i++) {
		x=MMSP::position(grid, i);
		SparseDistanceVoxel* Voxel=&(distance_grid(x));
		for (int j=0; j<Voxel->length(); j++) {
			int index = Voxel->index(j);
			grid(x).set(index) = Voxel->getValue(index);
		}
	}
} // approximate_cost


template <int dim, typename T>
void trace_pathway(const MMSP::grid<dim,MMSP::sparse<T> >& phase_grid,
                   const MMSP::grid<dim,SparseDistanceVoxel>& dist_grid,
                   const int id,
                   const MMSP::vector<int>& start,
                   const MMSP::vector<int>& finish,
                   std::map<int,std::vector<MMSP::vector<int> > >& path)
{
	// Pursue a Greedy algorithm (Method of Steepest Descents) to reach the Finish from the Start
	// ... unless a more attractive Finish is found
	const double epsilon=sqrt(std::numeric_limits<double>::epsilon());
	MMSP::vector<int> x(start);
	MMSP::vector<int> y(finish);
	for (int d=0; d<dim; d++) {
		MMSP::check_boundary(x[d],MMSP::x0(dist_grid,d),MMSP::x1(dist_grid,d),MMSP::b0(dist_grid,d),MMSP::b1(dist_grid,d));
		MMSP::check_boundary(y[d],MMSP::x0(dist_grid,d),MMSP::x1(dist_grid,d),MMSP::b0(dist_grid,d),MMSP::b1(dist_grid,d));
	}
	std::vector<MMSP::vector<int> > trial;

	// Construct sets of phases known at the start and end points.
	// Each point on a valid edge should share dim common grains.
	std::set<grainid> start_set;
	for (int i=0; i<phase_grid(x).length(); i++)
		if (phase_grid(x).value(i)>epsilon)
			start_set.insert(phase_grid(x).index(i));
	std::set<grainid> finish_set;
	for (int i=0; i<phase_grid(y).length(); i++)
		if (phase_grid(y).value(i)>epsilon)
			finish_set.insert(phase_grid(y).index(i));
	// Compare common phases in start and finish vertex locations.
	int common_phases=0;
	for (std::set<grainid>::const_iterator it_g=start_set.begin(); it_g!=start_set.end(); it_g++)
		if (finish_set.find(*it_g)!=finish_set.end())
			common_phases++;

	if (common_phases < dim) return; // dim cuts double edges, but may go too far

	if (dim==2) {
		while (radius<2,double>(x,y) > epsilon) {
			MMSP::vector<int> v(x);
			MMSP::vector<int> best(x);
			int closestVertex=id;
			for (v[1]=x[1]-1; v[1]<x[1]+2; v[1]++) {
				for (v[0]=x[0]-1; v[0]<x[0]+2; v[0]++) {
					MMSP::vector<int> p(v);
					for (int d=0; d<dim; d++)
						MMSP::check_boundary(p[d],MMSP::x0(dist_grid,d),MMSP::x1(dist_grid,d),MMSP::b0(dist_grid,d),MMSP::b1(dist_grid,d));

					const SparseDistanceVoxel* TestVoxel = &(dist_grid(p));
					const SparseDistanceVoxel* BestVoxel = &(dist_grid(best));

					bool isCloserToTarget = (  TestVoxel->getValue(id)            < BestVoxel->getValue(id) );
					bool isStillClosest = (    TestVoxel->getValue(closestVertex) < BestVoxel->getValue(id) + sqrt(epsilon));
					bool isCloserToClosest = ( TestVoxel->getValue(closestVertex) < BestVoxel->getValue(closestVertex));
					bool isAlongEdge = isCloserToTarget && isStillClosest && isCloserToClosest;

					if (isAlongEdge) {
						// Check whether another vertex is closer
						for (int i=0; i<TestVoxel->length(); i++) {
							int vid = TestVoxel->index(i);
							if (vid==closestVertex) continue;

							bool isCloserThanClosest = (  TestVoxel->getValue(vid) < TestVoxel->getValue(closestVertex) );
							bool isCloserThanTarget = (   TestVoxel->getValue(vid) < BestVoxel->getValue(vid) );
							bool isBetweenStartFinish = ( TestVoxel->getValue(vid) < dist_grid(start).getValue(vid) );

							if (isBetweenStartFinish && isCloserThanClosest && isCloserThanTarget)
								closestVertex=vid;
						}
						bool isSteppingCloser = ( TestVoxel->getValue(closestVertex) < BestVoxel->getValue(closestVertex) );
						if (isSteppingCloser)
							best=p;
					}
				}
			}
			if (best[0]==x[0] && best[1]==x[1]) {
				trial.clear(); // Walker got stuck: must be a vertex here.
				return;
			}
			x=best;
			trial.push_back(x);
		} // while radius>epsilon
	} else if (dim==3) {
		while (radius<3,double>(x,y) > epsilon) {
			MMSP::vector<int> v(x);
			MMSP::vector<int> best(x);
			int closestVertex=id;
			for (v[2]=x[2]-1; v[2]<x[2]+2; v[2]++) {
				for (v[1]=x[1]-1; v[1]<x[1]+2; v[1]++) {
					for (v[0]=x[0]-1; v[0]<x[0]+2; v[0]++) {
						MMSP::vector<int> p(v);
						for (int d=0; d<dim; d++)
							MMSP::check_boundary(p[d],MMSP::x0(dist_grid,d),MMSP::x1(dist_grid,d),MMSP::b0(dist_grid,d),MMSP::b1(dist_grid,d));

						const SparseDistanceVoxel* TestVoxel = &(dist_grid(p));
						const SparseDistanceVoxel* BestVoxel = &(dist_grid(best));
						bool isCloserToTarget = (  TestVoxel->getValue(id)            < BestVoxel->getValue(id) );
						bool isStillClosest = (    TestVoxel->getValue(closestVertex) < BestVoxel->getValue(id) + sqrt(epsilon));
						bool isCloserToClosest = ( TestVoxel->getValue(closestVertex) < BestVoxel->getValue(closestVertex));
						bool isAlongEdge = isCloserToTarget && isStillClosest && isCloserToClosest;
						if (isAlongEdge) {
							// Check whether another vertex is closer
							for (int i=0; i<TestVoxel->length(); i++) {
								int vid = TestVoxel->index(i);
								if (vid==closestVertex) continue;
								bool isCloserThanClosest = (  TestVoxel->getValue(vid) < TestVoxel->getValue(closestVertex) );
								bool isCloserThanTarget = (   TestVoxel->getValue(vid) < TestVoxel->getValue(id) );
								bool isBetweenStartFinish = ( TestVoxel->getValue(vid) < dist_grid(start).getValue(vid) );
								if (isCloserThanClosest && isCloserThanTarget && isBetweenStartFinish)
									closestVertex=vid;
							}
							bool isSteppingCloser = ( TestVoxel->getValue(closestVertex) < BestVoxel->getValue(closestVertex) );
							if (isSteppingCloser)
								best=p;
						}
					}
				}
			}
			if (best[0]==x[0] && best[1]==x[1] && best[2]==x[2]) {
				trial.clear(); // Walker got stuck: must be a vertex here.
				return;
			}
			x=best;
			trial.push_back(x);
		} // while radius>epsilon
	} else {
		std::cerr<<"Error: "<<dim<<"-D grid not supported."<<std::endl;
		exit(-1);
	}
	if (trial.size()==0) return;
	path[id] = trial;
} // trace_pathway


template <int dim, typename T>
void locate_edges(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                  const MMSP::grid<dim,SparseDistanceVoxel>& distance,
                  const std::vector<MMSP::vector<int> >& global_vertices,
                  std::map<grainid,std::map<grainid,std::vector<MMSP::vector<int> > > >& global_edges)
{
	global_edges.clear();
	// Choose a starting vertex
	for (int v0=0; v0<int(global_vertices.size()); v0++) {
		// Choose a destination vertex
		for (int j=0; j<distance(global_vertices[v0]).length(); j++) {
			int v1=distance(global_vertices[v0]).index(j);

			// Attempt walking from v0 to v1. Store valid walks in global_edges.

			if (v1==v0)
				continue; // no loops
			else if (global_edges[v0].find(v1)!=global_edges[v0].end() // edge exists
			         && global_edges[v0].find(v1)->second.size()>0)    // edge length is non-zero
				continue; // edge already exists
			else {
				trace_pathway(grid,distance,v1,global_vertices[v0],global_vertices[v1],global_edges[v0]);
				// Check whether path (v0,v1) was constructed.
				std::map<grainid,std::vector<MMSP::vector<int> > >::iterator m=global_edges[v0].find(v1);
				if (m!=global_edges[v0].end() && m->second.size()>0)
					global_edges[v1][v0]=global_edges[v0][v1]; // Copy the result to path (v0,v1).
			}

			// Check whether (v0,v1) is a valid edge. Erase zero-length edges.
			std::map<grainid,std::vector<MMSP::vector<int> > >::iterator it_edg=global_edges[v0].find(v1);
			if (it_edg!=global_edges[v0].end() && it_edg->second.size()==0)
				global_edges[v0].erase(it_edg);
		}
	}

	// Erase vertices with no valid edges.
	for (std::map<grainid,std::map<int,std::vector<MMSP::vector<int> > > >::iterator it_vtx=global_edges.begin(); it_vtx!=global_edges.end(); it_vtx++) {
		if (it_vtx->second.size()==0) {
			std::map<grainid,std::map<grainid,std::vector<MMSP::vector<int> > > >::iterator badv(it_vtx);
			it_vtx--;
			global_edges.erase(badv);
		}
	}
} // locate_edges


template <int dim, typename T>
void search_cycles(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                   const std::vector<MMSP::vector<int> >& global_vertices,
                   const std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >& global_edges,
                   std::set<std::vector<int> >& global_cycles)
{
	// Breadth-first search for chordless cycles of the graph. These are faces.
	// Note: Some cycles span multiple grains, and should not be taken seriously.

	// Source code modified from http://stackoverflow.com/questions/4022662/find-all-chordless-cycles-in-an-undirected-graph
	// to accept an adjacency-list, rather than an adjacency-matrix, and to limit cycle length.

	// The code generates triangular cycles IJKI and longer cycles JIK...MJ
	// Adjacent vertices must share dim common phases. This property is used to limit cycle length.
	// Every vertex in a cycle must share the same dim-1 phases: these are the grains forming the "face."

	for (int i=0; i<int(global_vertices.size())-2; i++) {
		std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >::const_iterator it_i=global_edges.find(i);
		if (it_i==global_edges.end())
			continue;
		// Store set of phases present at vertex I
		std::set<grainid> phases_i;
		for (int p=0; p<grid(global_vertices[i]).length(); p++)
			phases_i.insert(grid(global_vertices[i]).index(p));

		for (int j=i+1; j<int(global_vertices.size())-1; j++) {
			std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >::const_iterator it_j=global_edges.find(j);
			if (it_j==global_edges.end())
				continue;
			else if (it_i->second.find(j)==it_i->second.end() || (it_i->second.find(j))->second.size()==0) //if (!adjacency[i+j*global_vertices.size()]) -> !adj(i,j)
				continue;
			// Store set of phases present at vertex J
			std::set<grainid> phases_j;
			for (int p=0; p<grid(global_vertices[j]).length(); p++)
				phases_j.insert(grid(global_vertices[j]).index(p));
			int match_j=0;
			for (std::set<int>::const_iterator it_phi=phases_j.begin(); it_phi!=phases_j.end(); it_phi++)
				if (phases_i.find(*it_phi)!=phases_i.end())
					match_j++; // J shares phi with I
			if (match_j<dim)
				continue; // J cannot be adjacent to I, even if an edge exists, which it should not

			std::list<std::vector<int> > candidates;
			for (int k=j+1; k<int(global_vertices.size()); k++) {
				if (it_i->second.find(k)==it_i->second.end() || (it_i->second.find(k))->second.size()==0) //if (!adjacency[i+k*global_vertices.size()]) -> !adj(i,k)
					continue;
				// Compare phases in k to those in i. There must be at least dim matches.
				std::set<grainid> phases_k;
				for (int p=0; p<grid(global_vertices[k]).length(); p++)
					phases_k.insert(grid(global_vertices[k]).index(p));
				int match_ki=0, match_kj=0;
				for (std::set<int>::const_iterator it_phi=phases_k.begin(); it_phi!=phases_k.end(); it_phi++) {
					if (phases_i.find(*it_phi)!=phases_i.end())
						match_ki++; // K shares phi with I and J
					if (phases_j.find(*it_phi)!=phases_j.end())
						match_kj++;
				}
				// adj(i,j) is already established. The closed cycle IJKI must have adj(i,k) AND adj(j,k) as well.
				if (match_ki>=dim && match_kj>=dim && (it_j->second.find(k)!=it_j->second.end()) && ((it_j->second.find(k))->second.size()>0)) { // if (adjacency[j+k*global_vertices.size()]) -> adj(j,k)
					// k is adjacent to both i and j: legitimate triangular cycle found!
					std::vector<int> cycle;
					cycle.push_back(i);
					cycle.push_back(j);
					cycle.push_back(k);
					cycle.push_back(i);
					global_cycles.insert(cycle);
					continue;
				}
				// JIK...LMJ may be a cycle.
				std::vector<int> v;
				v.resize(3);
				v[0]=j;
				v[1]=i;
				v[2]=k;
				candidates.push_back(v);
			}
			while (!candidates.empty()) {
				// adj(i,j) is known.
				std::vector<int> v = candidates.front();
				candidates.pop_front();
				int l = v.back();
				// Store set of phases present at vertex L
				std::set<grainid> phases_l;
				for (int p=0; p<grid(global_vertices[l]).length(); p++)
					phases_l.insert(grid(global_vertices[l]).index(p));
				// Store set of phases present at vertex K (v[2])
				std::set<grainid> phases_k;
				for (int p=0; p<grid(global_vertices[v[2]]).length(); p++)
					phases_k.insert(grid(global_vertices[v[2]]).index(p));
				std::set<grainid> phases_jik;
				if (1) {
					std::set<grainid> phases;
					// collect all possible grainids
					for (int n=0; n<3; n++)
						for (int p=0; p<grid(global_vertices[v[n]]).length(); p++)
							phases.insert(grid(global_vertices[v[n]]).index(p));
					// store only those ids present at all three vertices starting this cycle
					for (std::set<grainid>::const_iterator it_phi=phases.begin(); it_phi!=phases.end(); it_phi++)
						if (phases_i.find(*it_phi)!=phases_i.end() && phases_j.find(*it_phi)!=phases_j.end() && phases_k.find(*it_phi)!=phases_k.end())
							phases_jik.insert(*it_phi);
				}
				if (phases_jik.size() < dim-1)
					continue; // jik cannot be on the same face
				for (int m=i+2; m<int(global_vertices.size()); m++) {
					if (std::find(v.begin(), v.end(), m) != v.end())
						continue;
					std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >::const_iterator it_m=global_edges.find(m);
					if (it_m==global_edges.end())
						continue;
					if (it_m->second.find(l)==it_m->second.end() || (it_m->second.find(l))->second.size()==0) //if (!adjacency[m+l*global_vertices.size()]) -> !adj(m,l)
						continue;
					std::set<grainid> phases_m;
					for (int p=0; p<grid(global_vertices[m]).length(); p++)
						phases_m.insert(grid(global_vertices[m]).index(p));
					// Compare phases in m to those in i, j, and k. At least dim-1 phases must be present for this to be a face.
					int match_m=0;
					for (std::set<int>::const_iterator it_phi=phases_m.begin(); it_phi!=phases_m.end(); it_phi++)
						if (phases_jik.find(*it_phi)!=phases_jik.end())
							match_m++;
					if (match_m<dim-1)
						continue;
					// Compare phases in m to those in l. At least dim phases must match for these to be adjacent.
					bool chord = false;
					for (int n=1; n<int(v.size())-1; n++)
						if (it_m->second.find(v[n])!=it_m->second.end() && (it_m->second.find(v[n]))->second.size()>0) // if (adjacency[m+v[n]*global_vertices.size()])
							chord = true;
					if (chord)
						continue;
					// If m is the penultimate vertex in the cycle, then m must share dim phases with j.
					match_m=0;
					for (std::set<int>::const_iterator it_phi=phases_m.begin(); it_phi!=phases_m.end(); it_phi++)
						if (phases_j.find(*it_phi)!=phases_j.end())
							match_m++;
					if (match_m>=dim && it_m->second.find(j)!=it_m->second.end() && (it_m->second.find(j))->second.size()>0) { //if (adjacency[m+j*global_vertices.size()])
						std::vector<int> cycle;
						for (int n=0; n<int(v.size()); n++)
							cycle.push_back(v[n]);
						cycle.push_back(m);
						cycle.push_back(v[0]);
						global_cycles.insert(cycle);
						continue;
					}
					// If m shares fewer than dim phases with j, that's OK- it may be an earlier link in the chain.
					std::vector<int> w(v);
					w.push_back(m);
					candidates.push_back(w);
				}
			}
		}
	}
	if (dim==2) printf("Found %4lu chordless cycles. ", global_cycles.size());
	else printf("Found %5lu chordless cycles. ", global_cycles.size());
}


template <int dim, typename T>
int assign_features(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                    const std::vector<MMSP::vector<int> >& global_vertices,
                    const std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >& global_edges,
                    const std::set<std::vector<int> >& global_cycles,
                    std::map<grainid,std::set<int> >& grainverts,
                    std::map<grainid,std::set<std::set<int> > >& grainedges,
                    std::map<grainid,std::set<std::vector<int> > >& graincycles,
                    std::map<grainid,std::map<int,std::vector<int> > >& grainfaces)
{
	// Assign vertex IDs (v) to constituent grains (using sparse data)
	for (unsigned int v=0; v<global_vertices.size(); v++) {
		const MMSP::sparse<T>* sp=&(grid(global_vertices[v]));
		for (int i=0; i<sp->length(); i++)
			grainverts[sp->index(i)].insert(v);
	} // assign vertices

	// Assign edges to grains based on the local set of vertices and the global set of edges
	//	For each vertex-pair v0 and v1 assigned to a grain, if there exists
	//	a global edge (v0,v1), then this grain shares that edge.
	std::vector<int> vertex_hits(global_vertices.size(),0); // keep track of vertices contributing to edges
	for (std::map<grainid,std::set<int> >::iterator it_grn=grainverts.begin(); it_grn!=grainverts.end(); it_grn++) {
		for (std::set<int>::const_iterator it_v0=it_grn->second.begin(); it_v0!=it_grn->second.end(); it_v0++) {
			// See if v0 has any associated edges
			const std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >::const_iterator it_edg_v0=global_edges.find(*it_v0);
			if (it_edg_v0==global_edges.end()) continue; // if no, move on (Should remove this vertex -- check if it happens!)
			for (std::set<int>::const_iterator it_v1=it_grn->second.begin(); it_v1!=it_grn->second.end(); it_v1++) {
				if (it_v0==it_v1) continue;
				// If the edge (v0,v1) exists, add it to this grain
				std::map<int,std::vector<MMSP::vector<int> > >::const_iterator it_edg_v1=it_edg_v0->second.find(*it_v1);
				if (it_edg_v1==it_edg_v0->second.end()) continue; // if no, move on
				std::set<int> edg;
				edg.insert(*it_v0);
				edg.insert(*it_v1);
				grainedges[it_grn->first].insert(edg);
				//vertex_hits[*it_v0]++;
			}
		}
	} // assign edges

	// Assign cycles to grains based on the local set of vertices
	for (std::map<grainid,std::set<int> >::iterator it_grn=grainverts.begin(); it_grn!=grainverts.end(); it_grn++) {
		for (std::set<std::vector<int> >::const_iterator it_cyc=global_cycles.begin(); it_cyc!=global_cycles.end(); it_cyc++) {
			bool isMyCycle=true;
			for (std::vector<int>::const_iterator it_vtx=it_cyc->begin(); isMyCycle && it_vtx!=it_cyc->end(); it_vtx++)
				if (it_grn->second.find(*it_vtx)==it_grn->second.end())
					isMyCycle=false;
			if (isMyCycle)
				graincycles[it_grn->first].insert(*it_cyc);
		}
	} // assign cycles

	// Assign shared cycles ("faces") to grains. Count the number of grains each cycle contributes to. Erase unshared cycles.
	for (std::map<grainid,std::set<std::vector<int> > >::iterator it_g0=graincycles.begin(); it_g0!=graincycles.end(); it_g0++) {
		for (std::set<std::vector<int> >::iterator it_c0=it_g0->second.begin(); it_c0!=it_g0->second.end(); it_c0++) {
			int cyclecount=0;
			std::vector<int> c0(*it_c0);
			c0.pop_back();
			std::sort(c0.begin(),c0.end());
			for (std::map<grainid,std::set<std::vector<int> > >::iterator it_g1=graincycles.begin(); it_g1!=graincycles.end(); it_g1++) {
				for (std::set<std::vector<int> >::iterator it_c1=it_g1->second.begin(); it_c1!=it_g1->second.end(); it_c1++) {
					if (it_g0->first==it_g1->first) continue;
					std::vector<int> c1(*it_c1);
					c1.pop_back();
					std::sort(c1.begin(),c1.end());
					if (c0==c1) {
						cyclecount++;
						if (dim==2) {// This cycle belongs to more than one grain. Erase it from the other, and keep searching.
							std::set<std::vector<int> >::iterator badc(it_c1);
							it_c1--;
							it_g1->second.erase(badc);
						}
						if (dim==3) // This cycle is shared with another grain. Rejoice!
							grainfaces[it_g0->first][it_g1->first]=*it_c0;
					}
				}
			}

			/*
			if (dim==3 && cyclecount==0) {
				// This cycle only contributed to this grain. Erase it.
				std::set<std::vector<int> >::iterator badc(it_c0);
				it_c0--;
				it_g0->second.erase(badc);
			}
			*/
			if (dim==2 && cyclecount>0) {
				// This cycle contributed to more than one grain. Erase it.
				std::set<std::vector<int> >::iterator badc(it_c0);
				it_c0--;
				it_g0->second.erase(badc);
			}
		}
	} // assign faces

	// Edges may exist that contribute to no cycles within this grain. Erase them.
	for (std::map<grainid,std::set<std::set<int> > >::iterator it_grn_edg=grainedges.begin(); it_grn_edg!=grainedges.end(); it_grn_edg++) {
		std::map<grainid,std::set<std::vector<int> > >::const_iterator it_grn_cyc=graincycles.find(it_grn_edg->first);
		if (it_grn_cyc==graincycles.end()) { // This grain does not exist: it has no cycles. Erase the whole thing.
			if (dim==3) {
				std::map<grainid,std::set<std::set<int> > >::iterator badg(it_grn_edg);
				it_grn_edg--;
				grainedges.erase(badg);
			}
			continue;
		}
		for (std::set<std::set<int> >::iterator it_edg=it_grn_edg->second.begin(); it_edg!=it_grn_edg->second.end(); it_edg++) {
			int edgecount=0; // track how many cycles this edge contributes to
			for (std::set<std::vector<int> >::const_iterator it_cyc=it_grn_cyc->second.begin(); it_cyc!=it_grn_cyc->second.end(); it_cyc++) {
				std::set<int>::const_iterator it_v = it_edg->begin();
				std::vector<int>::const_iterator it_v0 = std::find(it_cyc->begin(), it_cyc->end(), *it_v);
				it_v++;
				std::vector<int>::const_iterator it_v1 = std::find(it_cyc->begin(), it_cyc->end(), *it_v);
				if (it_v0!=it_cyc->end() && it_v1!=it_cyc->end())
					edgecount++; // both vertices on this edge belong to this cycle
			} // loop over cycles of this grain
			// DANGER
			//if (dim==3 && edgecount==0) {// Edge contributed to one cycle. Prune it.
			if (edgecount==0) {// Edge contributed to one cycle. Prune it.
				std::set<std::set<int> >::iterator bade(it_edg);
				it_edg--;
				it_grn_edg->second.erase(bade);
			}
		} // loop over edges of this grain
	} // prune edges

	// Count the number of grains each edge contributes to. If none, erase the edge.
	for (std::map<grainid,std::set<std::set<int> > >::iterator it_g0=grainedges.begin(); it_g0!=grainedges.end(); it_g0++) {
		for (std::set<std::set<int> >::iterator it_edg=it_g0->second.begin(); it_edg!=it_g0->second.end(); it_edg++) {
			int edgecount=0; // keep track of how many times this edge contributes to a grain
			for (std::map<grainid,std::set<std::set<int> > >::const_iterator it_g1=grainedges.begin(); it_g1!=grainedges.end(); it_g1++) {
				if (it_g0==it_g1) continue;
				std::set<std::set<int> >::const_iterator it_e1=it_g1->second.find(*it_edg);
				if (it_e1!=it_g1->second.end())
					edgecount++;
			}
			//if ((dim==2 && edgecount==0) || (dim==3 && edgecount<3)) {// Edge only contributed to one grain. Prune it.
			if (edgecount==0) {// Edge only contributed to one grain. Prune it.
				std::set<std::set<int> >::iterator bade(it_edg);
				it_edg--;
				it_g0->second.erase(bade);
			}
		}
	}

	// Vertices may exist that contribute to no edges within this grain. Erase them.
	for (std::map<grainid,std::set<int> >::iterator it_grn_vtx=grainverts.begin(); it_grn_vtx!=grainverts.end(); it_grn_vtx++) {
		std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn_edg=grainedges.find(it_grn_vtx->first);
		if (it_grn_edg==grainedges.end()) { // This grain does not exist: it has no cycles. Erase the whole thing.
			std::map<grainid,std::set<std::vector<int> > >::iterator it_cyc=graincycles.find(it_grn_vtx->first);
			if (it_cyc!=graincycles.end())
				graincycles.erase(it_cyc);
			std::map<grainid,std::map<int,std::vector<int> > >::iterator it_fac=grainfaces.find(it_grn_vtx->first);
			if (it_fac!=grainfaces.end())
				grainfaces.erase(it_fac);
			std::map<grainid,std::set<int> >::iterator badg(it_grn_vtx);
			it_grn_vtx--;
			grainverts.erase(badg);
			continue;
		}
		for (std::set<int>::iterator it_vtx=it_grn_vtx->second.begin(); it_vtx!=it_grn_vtx->second.end(); it_vtx++) {
			int vertcount=0; // track how many edges this edge contributes to
			for (std::set<std::set<int> >::const_iterator it_edg=it_grn_edg->second.begin(); it_edg!=it_grn_edg->second.end(); it_edg++) {
				std::set<int>::const_iterator it_v = it_edg->find(*it_vtx);
				if (it_v!=it_edg->end())
					vertcount++;
			} // loop over edges on this grain
			if (vertcount==0) {// Vertex contributed to no edges. Prune it.
				std::set<int>::iterator badv(it_vtx);
				it_vtx--;
				it_grn_vtx->second.erase(badv);
			}
		} // loop over vertices on this grain
	}

	// Count the number of grains each vertex contributes to. If none, erase the vertex.
	for (std::map<grainid,std::set<int> >::iterator it_g0=grainverts.begin(); it_g0!=grainverts.end(); it_g0++) {
		for (std::set<int>::iterator it_vtx=it_g0->second.begin(); it_vtx!=it_g0->second.end(); it_vtx++) {
			int vertcount=0; // keep track of how many times this vertex contributes to a grain
			for (std::map<grainid,std::set<int> >::const_iterator it_g1=grainverts.begin(); it_g1!=grainverts.end(); it_g1++) {
				if (it_g0==it_g1) continue;
				std::set<int>::const_iterator it_v1=it_g1->second.find(*it_vtx);
				if (it_v1!=it_g1->second.end()) {
					vertcount++;
				}
			}
			if (vertcount==0) {// Vertex only contributed to one grain. Prune it.
				std::set<int>::iterator badv(it_vtx);
				it_vtx--;
				it_g0->second.erase(badv);
			} else
				vertex_hits[*it_vtx]++;
		}
	}

	int plateaunic=0;
	for (unsigned int i=0; i<vertex_hits.size(); i++) {
		if (dim==2 && vertex_hits[i]==3) // each vertex should be touched by each of 3 incident cycles
			plateaunic++;
		else if (dim==3 && vertex_hits[i]==4) // each vertex should be touched by each of 4 incident cycles
			plateaunic++;
	}

	return plateaunic;
} // assign_features

std::map<int,int> regenerate_features(const std::map<grainid,std::map<int,std::vector<int> > >& grainfaces,
                                      std::map<grainid,std::set<int> >& grainverts,
                                      std::map<grainid,std::set<std::set<int> > >& grainedges,
                                      std::map<grainid,std::set<std::vector<int> > >& graincycles)
{
	grainverts.clear();
	grainedges.clear();
	//graincycles.clear();
	std::map<int,int> plateaunic;
	//for (std::map<grainid,std::map<int,std::vector<int> > >::const_iterator it_grn=grainfaces.begin(); it_grn!=grainfaces.end(); it_grn++) {
	//for (std::map<int,std::vector<int> >::const_iterator it_cyc=it_grn->second.begin(); it_cyc!=it_grn->second.end(); it_cyc++) {
	for (std::map<int,std::set<std::vector<int> > >::const_iterator it_fac=graincycles.begin(); it_fac!=graincycles.end(); it_fac++) {
		for (std::set<std::vector<int> >::const_iterator it_cyc=it_fac->second.begin(); it_cyc!=it_fac->second.end(); it_cyc++) {
			//graincycles[it_grn->first].insert(it_cyc->second);
			std::vector<int>::const_iterator it_v0=it_cyc->begin();
			std::vector<int>::const_iterator it_v1(it_v0);
			it_v1++;
			while (it_v1!=it_cyc->end()) {
				grainverts[it_fac->first].insert(*it_v1);
				std::set<int> edg;
				edg.insert(*it_v0);
				edg.insert(*it_v1);
				plateaunic[*it_v1]++;
				grainedges[it_fac->first].insert(edg);
				it_v0++;
				it_v1++;
				edg.clear();
			}
		}
	}
	//}

	return plateaunic;
} // regenerate_features

std::string progressbar(const int n, const int N, const char c)
{
	/*
		Usage: Call once before your for loop (over i) with n=0,
					 then at the end of each iteration with n=i+1
	*/
	std::stringstream bar;
	static unsigned long tstart;
	if (n==0) {
		tstart = time(NULL);
		bar<<" ["<<std::flush;
	} else if (n==N) {
		unsigned long deltat = time(NULL)-tstart;
		bar<<c<<"] "<<std::setw(2)<<std::right<<deltat/3600<<"h:"
		   <<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
		   <<std::setw(2)<<std::right<<deltat%60<<"s"
		   <<".\n";
	} else if ((20 * n) % (N-N%5) == 0) bar<<c;
	return bar.str();
} // progressbar


#endif // _DISTANCE_CPP_
