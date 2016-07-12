// File:    mmsp2graph.cpp
// Purpose: reads MMSP grid containing sparse vector of order parameters
//          writes grain-graph data using fast marching algorithm:
//          a better approach than the flattened representation of mmsp2topo

// Questions/Comments to trevor.keller@gmail.com (Trevor Keller)


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <map>
#include <utility>
#include <vector>
#include <ctime>
#include <zlib.h>

#include "MMSP.hpp"
#include "distance.cpp"
#include "schlegel.hpp"

int main(int argc, char* argv[])
{
	if ( argc != 3 ) {
		std::cout << "Usage: " << argv[0] << " input.dat output.tsv\n";
		return ( 1 );
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr<<"File input error: could not open "<<argv[1]<<".\n\n";
		exit(-1);
	}

	// read data type
	std::string type;
	getline(input,type,'\n');

	// grid type error check
	if (type.substr(0,4)!="grid") {
		std::cerr<<"File input error: "<<argv[1]<<" does not contain grid data."<<std::endl;
		exit(-1);
	}

	// parse data type
	bool float_type = (type.find("float") != std::string::npos);
	bool double_type = (type.find("double") != std::string::npos);
	bool sparse_type = (type.find("sparse") != std::string::npos);

	if (not float_type  and  not double_type  and  not sparse_type) {
		std::cerr << "File input error: support for "<<type<<" is not implemented." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input>>dim;
	if (dim!=2&&dim!=3) {
		std::cerr<<"ERROR: "<<dim<<"-D input not implemented.\n"<<std::endl;
		exit(-1);
	}

	// read number of fields
	int fields;
	input >> fields;

	// read grid sizes
	int x0[3] = {0, 0, 0};
	int x1[3] = {0, 0, 0};
	for (int i = 0; i < dim; i++)
		input >> x0[i] >> x1[i];

	// read cell spacing
	float dx[3] = {1.0, 1.0, 1.0};
	for (int i = 0; i < dim; i++)
		input >> dx[i];

	// ignore trailing endlines
	input.ignore(10, '\n');

	if (sparse_type) {
		static unsigned long tstart=time(NULL);
		if (float_type) {
			if (dim==2) {
				const std::string gridname(argv[1]);
				std::cout<<gridname.substr(gridname.find_last_of("/")+1,gridname.length())<<": "<<std::flush;
				MMSP::grid<2,MMSP::sparse<float> > grid(argv[1]);
				MMSP::grid<2,double> curvature(0, x0[0], x1[0], x0[1], x1[1]);
				MMSP::grid<2,SparseDistanceVoxel> distance(0, x0[0], x1[0], x0[1], x1[1]);
				MMSP::grid<2,MMSP::sparse<double> > distgrid(0, x0[0], x1[0], x0[1], x1[1]);
				MMSP::grid<2,unsigned char> pathways(0, x0[0], x1[0], x0[1], x1[1]);

				std::map<grainid,double> grainsizes;
				std::vector<MMSP::vector<int> > global_vertices;
				std::map<int,std::map<int,std::vector<MMSP::vector<int> > > > global_edges;
				std::set<std::vector<int> > global_cycles;

				locate_vertices(grid, grainsizes, global_vertices);
				determine_weights(grid,global_vertices,curvature);
				approximate_cost(grid, curvature, global_vertices, distance, distgrid);
				locate_edges(grid, distance, global_vertices, global_edges);
				search_cycles(grid, global_vertices, global_edges, global_cycles);
				std::map<grainid,std::set<int> > grainverts; // (grainid, vertexid)
				std::map<grainid,std::set<std::set<int> > > grainedges; // [grainid, {(v,v),(v,v),...})]
				std::map<grainid,std::set<std::vector<int> > > graincycles; // [grainid, {(v,v,v,v),(v,v,v,v),...}]
				std::map<grainid,std::map<grainid,std::vector<int> > > grainfaces; // [grainid, (grainid, {v, v, v...})]    (blank for 2D)
				int plateaunic=assign_features<2>(grid, global_vertices, global_edges, global_cycles, grainverts, grainedges, graincycles, grainfaces);



				// Export intermediate grids and report topology
				std::string altfile(argv[2]);
				altfile.resize(altfile.length()-4);
				altfile.append("_path.dat");
				for (int n=0; n<MMSP::nodes(pathways); ++n)
					pathways(n)=255;
				// Instead of iterating over global_edges, iterate over the grains. Paint edges and track hits on vertices.
				std::vector<int> vertex_hits;
				while (vertex_hits.size()<global_vertices.size()) vertex_hits.push_back(0);
				for (std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn=grainedges.begin(); it_grn!=grainedges.end(); it_grn++) {
					for (std::set<std::set<int> >::const_iterator it_edg=it_grn->second.begin(); it_edg!=it_grn->second.end(); it_edg++) {
						std::set<int>::const_iterator it_v0=it_edg->begin();
						std::set<int>::const_iterator it_v1(it_v0);
						it_v1++;
						//while (it_v1!=it_edg->end()) {
						// v0 and v1 are the endpoints of an edge. Count 'em.
						++vertex_hits[*it_v0];
						++vertex_hits[*it_v1];
						for (std::vector<MMSP::vector<int> >::const_iterator it_path=global_edges[*it_v0][*it_v1].begin(); it_path!=global_edges[*it_v0][*it_v1].end(); it_path++) {
							unsigned char c=pathways(*(it_path));
							if (c < static_cast<unsigned char>(65)) continue;
							else if (c > static_cast<unsigned char>(128)) pathways(*(it_path))=128;
							else pathways(*(it_path))=c-static_cast<unsigned char>(32);
						}
						//it_v0++;
						//it_v1++;
						//}
					}
				}
				for (int v=0; v<int(global_vertices.size()); v++) {
					pathways(global_vertices[v])=0;
					if (vertex_hits[v]/2 != 3) { // Make a box around problem vertices
						const MMSP::vector<int> x=global_vertices[v];
						MMSP::vector<int> y(x);
						for (y[1]=x[1]-2; y[1]<x[1]+3; y[1]+=2)
							for (y[0]=x[0]-2; y[0]<x[0]+3; y[0]+=2)
								pathways(y)=64;
					}
				}
				printf("%4d stable of %4lu total vertices (%3.0f%%). ", plateaunic, global_vertices.size(),100.0*plateaunic/global_vertices.size());
				MMSP::output(pathways, altfile.c_str());

				altfile.resize(altfile.length()-8);
				altfile.append("curv.dat");
				MMSP::grid<2,double> capped(curvature);
				#pragma omp parallel for
				for (int i=0; i<nodes(capped); i++)
					capped(i) = (curvature(i)>2.5) ? 2.5 : curvature(i);
				MMSP::output(capped, altfile.c_str());

				altfile.resize(altfile.length()-8);
				altfile.append("topo.csv");
				std::ofstream of(altfile.c_str());
				of<<"ID,V,v,e,f,x,errors,cycles\n";

				// Calculate fraction of valid grains based on the set of grain vertices
				int invariant=0;
				for (std::map<grainid,std::set<int> >::const_iterator i=grainverts.begin(); i!=grainverts.end(); i++) {
					const int id=i->first;
					//std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					const int v=i->second.size(); //(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int x=v-e;
					if (x==0)
						invariant++;
				}
				printf("V: %4d of %4lu (%3.0f%%), ", invariant, grainverts.size(), 100.0*invariant/grainverts.size());

				// Calculate fraction of valid grains and Output based on the set of grain edges
				invariant=0;
				double unitsz=1.0;
				for (int d=0; d<dim; d++) unitsz*=MMSP::dx(grid,d);
				for (std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.begin(); j!=grainedges.end(); j++) {
					const int id=j->first;
					std::map<grainid,double>::const_iterator szitr=grainsizes.find(id);
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const double A=unitsz*(szitr!=grainsizes.end())?szitr->second:0.0;
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=j->second.size();
					const int f=(k!=grainfaces.end())?k->second.size():0;
					const int x=v-e;
					const int c=(l!=graincycles.end())?l->second.size():0;
					if (x==0) {
						invariant++;
						of<<id<<','<<A<<','<<v<<','<<e<<','<<f<<','<<x<<','<<'0'<<','<<c<<'\n';
					} else {
						of<<id<<','<<A<<','<<v<<','<<e<<','<<f<<','<<x<<','<<'1'<<','<<c<<'\n';
					}
				}
				of.close();
				printf("E: %4d of %4lu (%3.0f%%) grains valid. %4lu seconds.\n", invariant, grainedges.size(), 100.0*invariant/grainedges.size(), time(NULL)-tstart);


			} else if (dim==3) {


				std::vector<unsigned long> cyclemonitor;
				const std::string gridname(argv[1]);
				std::cout<<gridname.substr(gridname.find_last_of("/")+1,gridname.length())<<": "<<std::flush;
				MMSP::grid<3,MMSP::sparse<float> > grid(argv[1]);

				MMSP::grid<3,double> curvature(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);
				MMSP::grid<3,SparseDistanceVoxel> distance(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);
				MMSP::grid<3,MMSP::sparse<double> > distgrid(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);
				MMSP::grid<3,unsigned char> pathways(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);

				std::map<grainid,double> grainsizes;
				std::vector<MMSP::vector<int> > global_vertices;
				std::map<int,std::map<int,std::vector<MMSP::vector<int> > > > global_edges;
				std::set<std::vector<int> > global_cycles;

				locate_vertices(grid, grainsizes, global_vertices);
				determine_weights(grid,global_vertices,curvature);
				approximate_cost(grid, curvature, global_vertices, distance, distgrid);
				locate_edges(grid, distance, global_vertices, global_edges);
				search_cycles(grid, global_vertices, global_edges, global_cycles);
				std::map<grainid,std::set<int> > grainverts; // (grainid, vertexid)
				std::map<grainid,std::set<std::set<int> > > grainedges; // [grainid, {(v,v),(v,v),...})]
				std::map<grainid,std::set<std::vector<int> > > graincycles; // [grainid, {(v,v,v,v),(v,v,v,v),...}]
				std::map<grainid,std::map<grainid,std::vector<int> > > grainfaces; // [grainid, (grainid, {v, v, v...})]
				int plateaunic=assign_features<3>(grid, global_vertices, global_edges, global_cycles, grainverts, grainedges, graincycles, grainfaces);


				// Export intermediate grids and report topology
				//MMSP::output(distgrid, argv[2]); // takes too long :/

				std::string altfile(argv[2]);
				altfile.resize(altfile.length()-4);
				altfile.append("_path.dat");
				for (int n=0; n<MMSP::nodes(pathways); n++)
					pathways(n)=255;
				// Paint pathways based on grain cycles
				for (std::map<grainid,std::set<std::vector<int> > >::const_iterator it_grn=graincycles.begin(); it_grn!=graincycles.end(); it_grn++) {
					for (std::set<std::vector<int> >::const_iterator it_cyc=it_grn->second.begin(); it_cyc!=it_grn->second.end(); it_cyc++) {
						assert(it_cyc->size()>2);
						std::vector<int>::const_iterator it_v0 = it_cyc->begin();
						std::vector<int>::const_iterator it_v1(it_v0);
						it_v1++;
						while (it_v1!=it_cyc->end()) {
							for (std::vector<MMSP::vector<int> >::const_iterator it_path=global_edges[*it_v0][*it_v1].begin(); it_path!=global_edges[*it_v0][*it_v1].end(); it_path++) {
								unsigned char c=pathways(*(it_path));
								if (c < static_cast<unsigned char>(65)) continue;
								else if (c > static_cast<unsigned char>(128)) pathways(*(it_path))=128;
								else pathways(*(it_path))=c-static_cast<unsigned char>(32);
							}
							++it_v0;
							++it_v1;
						}
					}
				}
				for (std::map<grainid,std::set<int> >::const_iterator it_vtx=grainverts.begin(); it_vtx!=grainverts.end(); it_vtx++)
					for (std::set<int>::const_iterator it_v=it_vtx->second.begin(); it_v!=it_vtx->second.end(); it_v++)
						pathways(global_vertices[*it_v])=250;
				MMSP::output(pathways, altfile.c_str());
				printf("%4d stable of %4lu total vertices (%3.0f%%). ", plateaunic, global_vertices.size(),100.0*plateaunic/global_vertices.size());


				altfile.resize(altfile.length()-8);
				altfile.append("topo.tsv");
				std::ofstream of(altfile.c_str());
				of<<"ID\tV\tv\te\tf\tx\terrors\tcycles\tep\tfx\tde\tdf\tWeinberg\tsymmetry";
				for (int p=3; p<6; p++)	of<<"\tp"<<p;
				of<<"\tDispersion\n";

				// Calculate fraction of valid grains based on the set of grain vertices
				int invariant=0;
				for (std::map<grainid,std::set<int> >::const_iterator i=grainverts.begin(); i!=grainverts.end(); i++) {
					const int id=i->first;
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const int v=i->second.size(); //(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					if (x==2)
						invariant++;
				}
				printf("V: %4d of %4lu (%3.0f%%), ", invariant, grainverts.size(), 100.0*invariant/grainverts.size());

				// Calculate fraction of valid grains and Output based on the set of grain edges
				invariant=0;
				double unitsz=1.0;
				for (int d=0; d<dim; d++) unitsz*=MMSP::dx(grid,d);
				for (std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.begin(); j!=grainedges.end(); j++) {
					const int id=j->first;
					std::map<grainid,double>::const_iterator szitr=grainsizes.find(id);
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const double V=unitsz*(szitr!=grainsizes.end())?szitr->second:0.0;
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=j->second.size();
					const int f=(k!=grainfaces.end())?k->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					of<<id<<'\t'<<V<<'\t'<<v<<'\t'<<e<<'\t'<<f<<'\t'<<x;
					if (x==2) {
						invariant++;
						of<<"\t0\t";
					} else
						of<<"\t1\t";
					of<<c<<'\t'<<3*v/2<<'\t'<<2+e-v<<'\t'<<e-(3*v/2)<<'\t'<<f-(2+e-v)<<std::flush;
					if (2.0*e/v<3.0)
						of<<"\t0\t0"; // graph must contain 2-connected vertices
					else if (i!=grainverts.end() && k!=grainfaces.end() && l!=graincycles.end()) {
						Schlegel schlegel_diagram(i,j,l);
						std::string timestep(altfile);
						timestep.resize(timestep.length()-9);
						int dotpos=timestep.find_last_of(".");
						timestep=timestep.substr(dotpos+1,timestep.length());
						const std::string wv=schlegel_diagram.get_weinberg_vector(timestep);
						// if (x==2)
						schlegel_diagram.export_eps(j, x, wv);
						of<<'\t'<<wv<<'\t'<<schlegel_diagram.get_sym();
					} else {
						of<<"\t0\t0";
					}
					// output subset of p-vec and dispersion
					if (c>0) {
						std::vector<int> pvec;
						for (std::set<std::vector<int> >::const_iterator p=l->second.begin(); p!=l->second.end(); p++) {
							while (pvec.size()<p->size())
								pvec.push_back(0);
							if (p->size()>0) pvec[p->size()-1]++; // cycles include origin twice
						}
						for (int p=3; p<6; p++)
							of<<'\t'<<pvec[p];
						float sig=0.0f;
						const float pbar=6.0-12.0/c;
						for (unsigned int i=0; i<pvec.size(); i++)
							sig += float(pvec[i])*(float(i)-pbar)*(float(i)-pbar);
						of<<'\t'<<std::sqrt(sig/float(c))<<'\n';
					} else {
						of<<"\t0\t0\t0\t0\n";
					}
				}
				of.close();
				printf("E: %4d of %4lu (%3.0f%%), ", invariant, grainedges.size(), 100.0*invariant/grainedges.size());

				// Calculate fraction of valid grains based on the set of grain faces
				// Export "cell adjacency"
				altfile.resize(altfile.length()-9);
				altfile.append(".tsv");
				altfile.replace(altfile.find("topo_"),4,"cell");
				of.open(altfile.c_str());
				invariant=0;
				for (std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.begin(); k!=grainfaces.end(); k++) {
					const int id=k->first;
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					if (x==2)
						invariant++;

					of<<k->first; // output central grain
					for (std::map<grainid,std::vector<int> >::const_iterator neigh=k->second.begin(); neigh!=k->second.end(); neigh++)
						of<<'\t'<<neigh->first;
					of<<'\n';
				}
				of.close();
				printf("F: %4d of %4lu (%3.0f%%), ", invariant, grainfaces.size(), 100.0*invariant/grainfaces.size());

				// Calculate fraction of valid grains based on the set of grain cycles
				invariant=0;
				for (std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.begin(); l!=graincycles.end(); l++) {
					const int id=l->first;
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					if (x==2)
						invariant++;
				}
				printf("C: %4d of %4lu (%3.0f%%) grains valid. %4lu seconds.\n", invariant, graincycles.size(), 100.0*invariant/graincycles.size(), time(NULL)-tstart);

				// Export radial data for |phi| near a vertex on the largest grain
				altfile.replace(altfile.find("cell"),4,"vert");
				of.open(altfile.c_str());
				double largest=0.0;
				std::map<grainid,std::set<int> >::const_iterator largest_grain=grainverts.begin();
				for (std::map<grainid,std::set<int> >::const_iterator j=grainverts.begin(); j!=grainverts.end(); j++) {
					std::map<grainid,double>::const_iterator k = grainsizes.find(j->first);
					if (k->second > largest) {
						largest = k->second;
						largest_grain = j;
					}
				}
				of<<"x\ty\tz\tr\t|phi|\n";
				std::set<int>::const_iterator vertex = largest_grain->second.begin();
				vertex++; // why not?
				MMSP::vector<int> v = global_vertices[*vertex];
				MMSP::vector<int> x(v);
				for (x[2]=v[2]-25; x[2]<v[2]+26; x[2]++) {
					for (x[1]=v[1]-25; x[1]<v[1]+26; x[1]++) {
						for (x[0]=v[0]-25; x[0]<v[0]+26; x[0]++) {
							MMSP::vector<int> p(x);
							double r=0.0;
							for (int d=0; d<dim; d++) {
								r += 1.0*(x[d]-v[d])*(x[d]-v[d]);
								MMSP::check_boundary(p[d],MMSP::x0(grid,d),MMSP::x1(grid,d),MMSP::b0(grid,d),MMSP::b1(grid,d));
							}
							of<<x[0]-v[0]<<'\t'<<x[1]-v[1]<<'\t'<<x[2]-v[2]<<'\t'<<std::sqrt(r)<<'\t'<<grid(p).getMagPhi()<<'\n';
						}
					}
				}
				of.close();

			}
		}
		else if (double_type) {
			if (dim==2) {
				const std::string gridname(argv[1]);
				std::cout<<gridname.substr(gridname.find_last_of("/")+1,gridname.length())<<": "<<std::flush;
				MMSP::grid<2,MMSP::sparse<double> > grid(argv[1]);
				MMSP::grid<2,double> curvature(1, x0[0], x1[0], x0[1], x1[1]);
				MMSP::grid<2,SparseDistanceVoxel> distance(1, x0[0], x1[0], x0[1], x1[1]);
				MMSP::grid<2,MMSP::sparse<double> > distgrid(0, x0[0], x1[0], x0[1], x1[1]);
				MMSP::grid<2,unsigned char> pathways(1, x0[0], x1[0], x0[1], x1[1]);

				std::map<grainid,double> grainsizes;
				std::vector<MMSP::vector<int> > global_vertices;
				std::map<int,std::map<int,std::vector<MMSP::vector<int> > > > global_edges;
				std::set<std::vector<int> > global_cycles;

				locate_vertices(grid, grainsizes, global_vertices);
				determine_weights(grid,global_vertices,curvature);
				approximate_cost(grid, curvature, global_vertices, distance, distgrid);
				locate_edges(grid, distance, global_vertices, global_edges);
				search_cycles(grid, global_vertices, global_edges, global_cycles);
				std::map<grainid,std::set<int> > grainverts; // (grainid, vertexid)
				std::map<grainid,std::set<std::set<int> > > grainedges; // [grainid, {(v,v),(v,v),...})]
				std::map<grainid,std::set<std::vector<int> > > graincycles; // [grainid, {(v,v,v,v),(v,v,v,v),...}]
				std::map<grainid,std::map<grainid,std::vector<int> > > grainfaces; // [grainid, (grainid, {v, v, v...})]    (blank for 2D)
				int plateaunic=assign_features<2>(grid, global_vertices, global_edges, global_cycles, grainverts, grainedges, graincycles, grainfaces);



				// Export intermediate grids and report topology
				std::string altfile(argv[2]);
				altfile.resize(altfile.length()-4);
				altfile.append("_path.dat");
				for (int n=0; n<MMSP::nodes(pathways); ++n)
					pathways(n)=255;
				// Instead of iterating over global_edges, iterate over the grains. Paint edges and track hits on vertices.
				std::vector<int> vertex_hits;
				while (vertex_hits.size()<global_vertices.size()) vertex_hits.push_back(0);
				for (std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn=grainedges.begin(); it_grn!=grainedges.end(); it_grn++) {
					for (std::set<std::set<int> >::const_iterator it_edg=it_grn->second.begin(); it_edg!=it_grn->second.end(); it_edg++) {
						std::set<int>::const_iterator it_v0=it_edg->begin();
						std::set<int>::const_iterator it_v1(it_v0);
						it_v1++;
						//while (it_v1!=it_edg->end()) {
						// v0 and v1 are the endpoints of an edge. Count 'em.
						++vertex_hits[*it_v0];
						++vertex_hits[*it_v1];
						for (std::vector<MMSP::vector<int> >::const_iterator it_path=global_edges[*it_v0][*it_v1].begin(); it_path!=global_edges[*it_v0][*it_v1].end(); it_path++) {
							unsigned char c=pathways(*(it_path));
							if (c < static_cast<unsigned char>(65)) continue;
							else if (c > static_cast<unsigned char>(128)) pathways(*(it_path))=128;
							else pathways(*(it_path))=c-static_cast<unsigned char>(32);
						}
						//it_v0++;
						//it_v1++;
						//}
					}
				}
				for (int v=0; v<int(global_vertices.size()); v++) {
					pathways(global_vertices[v])=0;
					if (vertex_hits[v]/2 != 3) { // Make a box around problem vertices
						const MMSP::vector<int> x=global_vertices[v];
						MMSP::vector<int> y(x);
						for (y[1]=x[1]-2; y[1]<x[1]+3; y[1]+=2)
							for (y[0]=x[0]-2; y[0]<x[0]+3; y[0]+=2)
								pathways(y)=64;
					}
				}
				printf("%4d stable of %4lu total vertices (%3.0f%%). ", plateaunic, global_vertices.size(),100.0*plateaunic/global_vertices.size());
				MMSP::output(pathways, altfile.c_str());

				altfile.resize(altfile.length()-8);
				altfile.append("curv.dat");
				MMSP::grid<2,double> capped(curvature);
				#pragma omp parallel for
				for (int i=0; i<nodes(capped); i++)
					capped(i) = (curvature(i)>2.5) ? 2.5 : curvature(i);
				MMSP::output(capped, altfile.c_str());

				altfile.resize(altfile.length()-8);
				altfile.append("topo.csv");
				std::ofstream of(altfile.c_str());
				of<<"ID,V,v,e,f,x,errors,cycles\n";

				// Calculate fraction of valid grains based on the set of grain vertices
				int invariant=0;
				for (std::map<grainid,std::set<int> >::const_iterator i=grainverts.begin(); i!=grainverts.end(); i++) {
					const int id=i->first;
					//std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					const int v=i->second.size(); //(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int x=v-e;
					if (x==0)
						invariant++;
				}
				printf("V: %4d of %4lu (%3.0f%%), ", invariant, grainverts.size(), 100.0*invariant/grainverts.size());

				// Calculate fraction of valid grains and Output based on the set of grain edges
				invariant=0;
				double unitsz=1.0;
				for (int d=0; d<dim; d++) unitsz*=MMSP::dx(grid,d);
				for (std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.begin(); j!=grainedges.end(); j++) {
					const int id=j->first;
					std::map<grainid,double>::const_iterator szitr=grainsizes.find(id);
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const double A=unitsz*(szitr!=grainsizes.end())?szitr->second:0.0;
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=j->second.size();
					const int f=(k!=grainfaces.end())?k->second.size():0;
					const int x=v-e;
					const int c=(l!=graincycles.end())?l->second.size():0;
					if (x==0) {
						invariant++;
						of<<id<<','<<A<<','<<v<<','<<e<<','<<f<<','<<x<<','<<'0'<<','<<c<<'\n';
					} else {
						of<<id<<','<<A<<','<<v<<','<<e<<','<<f<<','<<x<<','<<'1'<<','<<c<<'\n';
					}
				}
				of.close();
				printf("E: %4d of %4lu (%3.0f%%) grains valid. %4lu seconds.\n", invariant, grainedges.size(), 100.0*invariant/grainedges.size(), time(NULL)-tstart);


			} else if (dim==3) {


				std::vector<unsigned long> cyclemonitor;
				const std::string gridname(argv[1]);
				std::cout<<gridname.substr(gridname.find_last_of("/")+1,gridname.length())<<": "<<std::flush;
				MMSP::grid<3,MMSP::sparse<double> > grid(argv[1]);

				MMSP::grid<3,double> curvature(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);
				MMSP::grid<3,SparseDistanceVoxel> distance(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);
				MMSP::grid<3,MMSP::sparse<double> > distgrid(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);
				MMSP::grid<3,unsigned char> pathways(0, x0[0], x1[0], x0[1], x1[1], x0[2], x1[2]);

				std::map<grainid,double> grainsizes;
				std::vector<MMSP::vector<int> > global_vertices;
				std::map<int,std::map<int,std::vector<MMSP::vector<int> > > > global_edges;
				std::set<std::vector<int> > global_cycles;

				locate_vertices(grid, grainsizes, global_vertices);
				determine_weights(grid,global_vertices,curvature);
				approximate_cost(grid, curvature, global_vertices, distance, distgrid);
				locate_edges(grid, distance, global_vertices, global_edges);
				search_cycles(grid, global_vertices, global_edges, global_cycles);
				std::map<grainid,std::set<int> > grainverts; // (grainid, vertexid)
				std::map<grainid,std::set<std::set<int> > > grainedges; // [grainid, {(v,v),(v,v),...})]
				std::map<grainid,std::set<std::vector<int> > > graincycles; // [grainid, {(v,v,v,v),(v,v,v,v),...}]
				std::map<grainid,std::map<grainid,std::vector<int> > > grainfaces; // [grainid, (grainid, {v, v, v...})]
				int plateaunic=assign_features<3>(grid, global_vertices, global_edges, global_cycles, grainverts, grainedges, graincycles, grainfaces);


				// Export intermediate grids and report topology
				//MMSP::output(distgrid, argv[2]); // takes too long :/

				std::string altfile(argv[2]);
				altfile.resize(altfile.length()-4);
				altfile.append("_path.dat");
				for (int n=0; n<MMSP::nodes(pathways); n++)
					pathways(n)=255;
				// Paint pathways based on grain cycles
				for (std::map<grainid,std::set<std::vector<int> > >::const_iterator it_grn=graincycles.begin(); it_grn!=graincycles.end(); it_grn++) {
					for (std::set<std::vector<int> >::const_iterator it_cyc=it_grn->second.begin(); it_cyc!=it_grn->second.end(); it_cyc++) {
						assert(it_cyc->size()>2);
						std::vector<int>::const_iterator it_v0 = it_cyc->begin();
						std::vector<int>::const_iterator it_v1(it_v0);
						it_v1++;
						while (it_v1!=it_cyc->end()) {
							for (std::vector<MMSP::vector<int> >::const_iterator it_path=global_edges[*it_v0][*it_v1].begin(); it_path!=global_edges[*it_v0][*it_v1].end(); it_path++) {
								unsigned char c=pathways(*(it_path));
								if (c < static_cast<unsigned char>(65)) continue;
								else if (c > static_cast<unsigned char>(128)) pathways(*(it_path))=128;
								else pathways(*(it_path))=c-static_cast<unsigned char>(32);
							}
							++it_v0;
							++it_v1;
						}
					}
				}
				for (std::map<grainid,std::set<int> >::const_iterator it_vtx=grainverts.begin(); it_vtx!=grainverts.end(); it_vtx++)
					for (std::set<int>::const_iterator it_v=it_vtx->second.begin(); it_v!=it_vtx->second.end(); it_v++)
						pathways(global_vertices[*it_v])=250;
				MMSP::output(pathways, altfile.c_str());
				printf("%4d stable of %4lu total vertices (%3.0f%%). ", plateaunic, global_vertices.size(),100.0*plateaunic/global_vertices.size());


				altfile.resize(altfile.length()-8);
				altfile.append("topo.tsv");
				std::ofstream of(altfile.c_str());
				of<<"ID\tV\tv\te\tf\tx\terrors\tcycles\tep\tfx\tde\tdf\tWeinberg\tsymmetry";
				for (int p=3; p<6; p++)	of<<"\tp"<<p;
				of<<"\tDispersion\n";

				// Calculate fraction of valid grains based on the set of grain vertices
				int invariant=0;
				for (std::map<grainid,std::set<int> >::const_iterator i=grainverts.begin(); i!=grainverts.end(); i++) {
					const int id=i->first;
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const int v=i->second.size(); //(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					if (x==2)
						invariant++;
				}
				printf("V: %4d of %4lu (%3.0f%%), ", invariant, grainverts.size(), 100.0*invariant/grainverts.size());

				// Calculate fraction of valid grains and Output based on the set of grain edges
				invariant=0;
				double unitsz=1.0;
				for (int d=0; d<dim; d++) unitsz*=MMSP::dx(grid,d);
				for (std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.begin(); j!=grainedges.end(); j++) {
					const int id=j->first;
					std::map<grainid,double>::const_iterator szitr=grainsizes.find(id);
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const double V=unitsz*(szitr!=grainsizes.end())?szitr->second:0.0;
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=j->second.size();
					const int f=(k!=grainfaces.end())?k->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					of<<id<<'\t'<<V<<'\t'<<v<<'\t'<<e<<'\t'<<f<<'\t'<<x;
					if (x==2) {
						invariant++;
						of<<"\t0\t";
					} else
						of<<"\t1\t";
					of<<c<<'\t'<<3*v/2<<'\t'<<2+e-v<<'\t'<<e-(3*v/2)<<'\t'<<f-(2+e-v)<<std::flush;
					if (2.0*e/v<3.0)
						of<<"\t0\t0"; // graph must contain 2-connected vertices
					else if (i!=grainverts.end() && k!=grainfaces.end() && l!=graincycles.end()) {
						Schlegel schlegel_diagram(i,j,l);
						std::string timestep(altfile);
						timestep.resize(timestep.length()-9);
						int dotpos=timestep.find_last_of(".");
						timestep=timestep.substr(dotpos+1,timestep.length());
						const std::string wv=schlegel_diagram.get_weinberg_vector(timestep);
						// if (x==2)
						schlegel_diagram.export_eps(j, x, wv);
						of<<'\t'<<wv<<'\t'<<schlegel_diagram.get_sym();
					} else {
						of<<"\t0\t0";
					}
					// output subset of p-vec and dispersion
					if (c>0) {
						std::vector<int> pvec;
						for (std::set<std::vector<int> >::const_iterator p=l->second.begin(); p!=l->second.end(); p++) {
							while (pvec.size()<p->size())
								pvec.push_back(0);
							if (p->size()>0) pvec[p->size()-1]++; // cycles include origin twice
						}
						for (int p=3; p<6; p++)
							of<<'\t'<<pvec[p];
						float sig=0.0f;
						const float pbar=6.0-12.0/c;
						for (unsigned int i=0; i<pvec.size(); i++)
							sig += float(pvec[i])*(float(i)-pbar)*(float(i)-pbar);
						of<<'\t'<<std::sqrt(sig/float(c))<<'\n';
					} else {
						of<<"\t0\t0\t0\t0\n";
					}
				}
				of.close();
				printf("E: %4d of %4lu (%3.0f%%), ", invariant, grainedges.size(), 100.0*invariant/grainedges.size());

				// Calculate fraction of valid grains based on the set of grain faces
				// Export "cell adjacency"
				altfile.resize(altfile.length()-9);
				altfile.append(".tsv");
				altfile.replace(altfile.find("topo_"),4,"cell");
				of.open(altfile.c_str());
				invariant=0;
				for (std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.begin(); k!=grainfaces.end(); k++) {
					const int id=k->first;
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.find(id);
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					if (x==2)
						invariant++;

					of<<k->first; // output central grain
					for (std::map<grainid,std::vector<int> >::const_iterator neigh=k->second.begin(); neigh!=k->second.end(); neigh++)
						of<<'\t'<<neigh->first;
					of<<'\n';
				}
				of.close();
				printf("F: %4d of %4lu (%3.0f%%), ", invariant, grainfaces.size(), 100.0*invariant/grainfaces.size());

				// Calculate fraction of valid grains based on the set of grain cycles
				invariant=0;
				for (std::map<grainid,std::set<std::vector<int> > >::const_iterator l=graincycles.begin(); l!=graincycles.end(); l++) {
					const int id=l->first;
					std::map<grainid,std::set<int> >::const_iterator i=grainverts.find(id);
					std::map<grainid,std::set<std::set<int> > >::const_iterator j=grainedges.find(id);
					std::map<grainid,std::map<grainid,std::vector<int> > >::const_iterator k=grainfaces.find(id);
					const int v=(i!=grainverts.end())?i->second.size():0;
					const int e=(j!=grainedges.end())?j->second.size():0;
					const int c=(l!=graincycles.end())?l->second.size():0;
					const int x=c-e+v;
					if (x==2)
						invariant++;
				}
				printf("C: %4d of %4lu (%3.0f%%) grains valid. %4lu seconds.\n", invariant, graincycles.size(), 100.0*invariant/graincycles.size(), time(NULL)-tstart);

				// Export radial data for |phi| near a vertex on the largest grain
				altfile.replace(altfile.find("cell"),4,"vert");
				of.open(altfile.c_str());
				double largest=0.0;
				std::map<grainid,std::set<int> >::const_iterator largest_grain=grainverts.begin();
				for (std::map<grainid,std::set<int> >::const_iterator j=grainverts.begin(); j!=grainverts.end(); j++) {
					std::map<grainid,double>::const_iterator k = grainsizes.find(j->first);
					if (k->second > largest) {
						largest = k->second;
						largest_grain = j;
					}
				}
				of<<"x\ty\tz\tr\t|phi|\n";
				std::set<int>::const_iterator vertex = largest_grain->second.begin();
				vertex++; // why not?
				MMSP::vector<int> v = global_vertices[*vertex];
				MMSP::vector<int> x(v);
				for (x[2]=v[2]-25; x[2]<v[2]+26; x[2]++) {
					for (x[1]=v[1]-25; x[1]<v[1]+26; x[1]++) {
						for (x[0]=v[0]-25; x[0]<v[0]+26; x[0]++) {
							MMSP::vector<int> p(x);
							double r=0.0;
							for (int d=0; d<dim; d++) {
								r += 1.0*(x[d]-v[d])*(x[d]-v[d]);
								MMSP::check_boundary(p[d],MMSP::x0(grid,d),MMSP::x1(grid,d),MMSP::b0(grid,d),MMSP::b1(grid,d));
							}
							of<<x[0]-v[0]<<'\t'<<x[1]-v[1]<<'\t'<<x[2]-v[2]<<'\t'<<std::sqrt(r)<<'\t'<<grid(p).getMagPhi()<<'\n';
						}
					}
				}
				of.close();

			}
		}

	} else {

		std::cerr<<"File input error: "<<type<<" not implemented."<<std::endl;
		exit(-1);
	}

	return 0;
}
