// distance.h
// Definitions of methods to propagate vertex distances in MMSP grids
#ifndef _DISTANCE_H_
#define _DISTANCE_H_

typedef int grainid;
struct flt2 { float x; float y; };

template <int dim, typename T>
double radius(const MMSP::vector<T>& a, const MMSP::vector<T>& b);

template <int dim, typename T>
int radiussq(const MMSP::vector<T>& a, const MMSP::vector<T>& b);

MMSP::vector<int> getPosition(const SparseDistanceVoxel& dv);

template<int dim, typename T, typename U>
void propagate_cost(const SparseDistanceVoxel* core_voxel,
                    const int index,
                    const MMSP::grid<dim,MMSP::sparse<U> >& grid,
                    const MMSP::grid<dim,T>& cost,
                    MMSP::grid<dim,SparseDistanceVoxel>& distimage,
                    SparsePriorityQueue& queue);

template<int dim, typename T>
bool isLocalMin(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                const MMSP::vector<int> position);

template<int dim, typename T>
void locate_vertices(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                     std::vector<MMSP::vector<int> >& global_vertices);

template<int dim, typename T>
void determine_weights(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                       const std::vector<MMSP::vector<int> >& global_vertices,
                       MMSP::grid<dim,double>& cost);

template<int dim, typename T, typename U>
void approximate_cost(const MMSP::grid<dim,MMSP::sparse<U> >& phi,
                      const MMSP::grid<dim,T>& cost,
                      const std::vector<MMSP::vector<int> >& global_vertices,
                      MMSP::grid<dim,SparseDistanceVoxel>& distance_grid,
                      MMSP::grid<dim,MMSP::sparse<T> >& grid);

template <int dim, typename T>
void trace_pathway(const MMSP::grid<dim,MMSP::sparse<T> >& phase_grid,
                   const MMSP::grid<dim,SparseDistanceVoxel>& dist_grid,
                   const int id,
                   const MMSP::vector<int>& start,
                   const MMSP::vector<int>& finish,
                   std::map<int,std::vector<MMSP::vector<int> > >& path);

template <int dim, typename T>
void locate_edges(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                  const MMSP::grid<dim,SparseDistanceVoxel>& distance,
                  const std::vector<MMSP::vector<int> >& global_vertices,
                  std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >& global_edges);

template <int dim, typename T>
void search_cycles(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                   const std::vector<MMSP::vector<int> >& global_vertices,
                   const std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >& global_edges,
                   std::set<std::vector<int> >& global_cycles);

template <int dim, typename T>
int assign_features(const MMSP::grid<dim,MMSP::sparse<T> >& grid,
                    const std::vector<MMSP::vector<int> >& global_vertices,
                    const std::map<int,std::map<int,std::vector<MMSP::vector<int> > > >& global_edges,
                    const std::set<std::vector<int> >& global_cycles,
                    std::map<grainid,std::set<int> >& grainverts,
                    std::map<grainid,std::set<std::set<int> > >& grainedges,
                    std::map<grainid,std::set<std::vector<int> > >& graincycles,
                    std::map<grainid,std::map<grainid,std::vector<int> > >& grainfaces);

std::map<int,int> regenerate_features(const std::map<grainid,std::map<int,std::vector<int> > >& grainfaces,
                                      std::map<grainid,std::set<int> >& grainverts,
                                      std::map<grainid,std::set<std::set<int> > >& grainedges,
                                      std::map<grainid,std::set<std::vector<int> > >& graincycles);

void export_schlegels(const std::map<grainid,std::set<int> >& grainverts,
                      const std::map<grainid,std::set<std::set<int> > >& grainedges,
                      const std::map<grainid,std::set<std::vector<int> > >& graincycles,
                      const std::string& timestep);

std::string progressbar(const int n, const int N, const char c);

#endif // _DISTANCE_H_
