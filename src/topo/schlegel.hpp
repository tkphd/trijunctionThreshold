#ifndef _SCHLEGEL_HPP_
#define _SCHLEGEL_HPP_

#include<map>
#include<string>
#include<set>
#include<vector>
#include<sys/stat.h>

const float r=100.0f;
const int iter=16384;

bool fileExists(const std::string& file) {
	struct stat buf;
	return (stat(file.c_str(), &buf) == 0);
}

class Schlegel
{
public:
	Schlegel() {
		std::cerr<<"Error in Schlegel diagram: default constructor must not be used!"<<std::endl;
		exit(-1);
	}
	Schlegel(const std::map<grainid,std::set<int> >::const_iterator it_grn_vtx,
	         const std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn_edg,
	         const std::map<grainid,std::set<std::vector<int> > >::const_iterator it_grn_cyc);
	void export_eps(const std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn_edg, const int x, const std::string& wv);
	std::string get_weinberg_vector(const std::string& timestep);
	int get_sym(){return wv_symmetries;}

private:
	grainid id;
	bool two_edge;
	std::map<int,std::set<int> > adjacency;
	std::map<int,flt2> positions;
	std::vector<int> weinberg_vector;
	int wv_rotations, wv_symmetries;
	//std::map<int,int> weinberg_map; // map Weinberg labels (0-26) to vertex IDs
	void weinberg_recurse(int vert_from, int vert_to, const std::map<int,std::vector<int> >& edge_map, std::vector<int> &sequence, std::map<int,int>& verts_visited, std::map<int,std::map<int,int> >& edges_visited, bool mirror);
	void relabel_weinberg(std::vector<int>& sequence);
}; // class Schlegel

Schlegel::Schlegel(const std::map<grainid,std::set<int> >::const_iterator it_grn_vtx,
                   const std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn_edg,
                   const std::map<grainid,std::set<std::vector<int> > >::const_iterator it_grn_cyc)
{
	wv_rotations=0;
	wv_symmetries=0;

	/*
		Algorithm for generating Schlegel diagrams after E. Lazar (http://www.seas.upenn.edu/~mlazar/schlegels.html)

		1. Select a face, preferably with the most edges. Space its vertices equally around the unit circle. Lock positions.
		2. Place the remaining vertices at the origin.
		3. While things are moving, or fewer than 100 iterations, or both:
			a. Select a mobile vertex (ie, not on the unit circle).
			b. Find the positions of the vertex's neighbors (using adjacency matrix).
			c. Move the vertex to the middle of its neighbors.
		4. Export to Encapsulated PostScript.
	*/

	this->id=it_grn_edg->first;

	// Revise valid edges to "true"
	this->adjacency.clear();
	for (std::set<std::set<int> >::const_iterator it_edg=it_grn_edg->second.begin(); it_edg!=it_grn_edg->second.end(); it_edg++) {
		std::set<int>::const_iterator it_v0=it_edg->begin();
		std::set<int>::const_iterator it_v1(it_v0);
		it_v1++;
		assert (it_v1!=it_edg->end());
		this->adjacency[*it_v0].insert(*it_v1);
		this->adjacency[*it_v1].insert(*it_v0);
	}

	// Search for 2-valent vertices: assume they represent algorithm failures.
	two_edge=false;
	for (std::map<int,std::set<int> >::const_iterator it_adj=adjacency.begin(); !two_edge && it_adj!=adjacency.end(); it_adj++)
		if (it_adj->second.size()<3)
			two_edge=true;

	//	====================================================================
	/*
	// Experimental: Add cycle pairs to adjacencies
	for (std::set<std::vector<int> >::const_iterator it_cyc=it_grn_cyc->second.begin(); it_cyc!=it_grn_cyc->second.end(); it_cyc++) {
		std::vector<int>::const_iterator it_v0=it_cyc->begin();
		std::vector<int>::const_iterator it_v1(it_v0);
		it_v1++;
		while (it_v1!=it_cyc->end()) {
			this->adjacency[*it_v0].insert(*it_v1);
			this->adjacency[*it_v1].insert(*it_v0);
			it_v0++;
			it_v1++;
		}
	}
	*/
	//	====================================================================

	// Determine bounding face as the one with the most vertices
	std::set<std::vector<int> >::const_iterator max_face=it_grn_cyc->second.begin();
	for (std::set<std::vector<int> >::const_iterator it_cyc=it_grn_cyc->second.begin(); it_cyc!=it_grn_cyc->second.end(); it_cyc++)
		if (it_cyc->size() > max_face->size())
			max_face=it_cyc;

	this->positions.clear();
	std::map<int,bool> fixedpos;
	flt2 origin;
	origin.x=origin.y=0.0f;
	for (std::set<int>::const_iterator it_vtx=it_grn_vtx->second.begin(); it_vtx!=it_grn_vtx->second.end(); it_vtx++) {
		if (this->adjacency.find(*it_vtx)==this->adjacency.end()) continue;
		this->positions[*it_vtx]=origin;
		fixedpos[*it_vtx]=false;
	}
	float theta=0.0f;
	for (std::vector<int>::const_iterator it_vtx=max_face->begin(); it_vtx!=max_face->end(); it_vtx++) {
		if (this->adjacency.find(*it_vtx)==this->adjacency.end()) continue;
		fixedpos[*it_vtx]=true;
		this->positions[*it_vtx].x = r*sin(theta);
		this->positions[*it_vtx].y = r*cos(theta);
		theta += 2.0*M_PI/(max_face->size()-1);
	}

	// Move each free vertex to the midpoint of its neighbors
	for (int i=0; i<iter; i++) {
		for (std::map<int,flt2>::iterator it_pos=this->positions.begin(); it_pos!=this->positions.end(); it_pos++) {
			if (fixedpos[it_pos->first]) continue; // skip vertices on the border
			flt2 newpos;
			newpos.x=0.0f;
			newpos.y=0.0f;
			for (std::set<int>::const_iterator it_vtx=this->adjacency[it_pos->first].begin(); it_vtx!=this->adjacency[it_pos->first].end(); it_vtx++) {
				newpos.x += this->positions[*it_vtx].x/this->adjacency[it_pos->first].size();
				newpos.y += this->positions[*it_vtx].y/this->adjacency[it_pos->first].size();
			}
			it_pos->second.x=newpos.x;
			it_pos->second.y=newpos.y;
		}
	}
}

void Schlegel::relabel_weinberg(std::vector<int> &sequence)
{
	std::map<int,int> association;
	int next = 1;

	for (unsigned int k=0; k<sequence.size(); k++) {
		std::map<int,int>::iterator itr=association.find(sequence[k]);
		if (itr==association.end()) {
			association[sequence[k]]=next;
			sequence[k]=next;
			next++;
		} else {
			sequence[k]=itr->second;
		}
	}
}

std::string Schlegel::get_weinberg_vector(const std::string& timestep)
{
	// Construct Weinberg vector from vertex positions (Schlegel diagram)

	if (this->positions.size()==0 || two_edge)
		return std::string("0");

	//assert (this->adjacency.size() == this->positions.size());

	// construct maps of vertex sequences in counter-clockwise order using theta = arctan(y/x)
	std::map<int,std::vector<int> > edge_map;

	//for (std::map<int,flt2>::const_iterator it_core=this->positions.begin(); it_core!=this->positions.end(); it_core++) {
	//std::map<int,std::set<int> >::const_iterator it_adj=this->adjacency.find(it_core->first);
	for (std::map<int,std::set<int> >::const_iterator it_adj=this->adjacency.begin(); it_adj!=this->adjacency.end(); it_adj++) {
		std::map<int,flt2>::const_iterator it_core=this->positions.find(it_adj->first); // "core" is the cetral vertex
		std::map<double,int> thetas; // sequence: crossprod, from, next
		for (std::set<int>::const_iterator it_branch=it_adj->second.begin(); it_branch!=it_adj->second.end(); it_branch++) {
			// "branch" is one of the adjacent vertices
			if (*it_branch==it_core->first) continue;
			double y = this->positions[*it_branch].y-this->positions[it_core->first].y;
			double x = this->positions[*it_branch].x-this->positions[it_core->first].x;
			double theta = atan2(y, x); // implementation of arctan which avoids dividing by zero
			thetas[theta]=*it_branch;
		}
		// construct counter-clockwise sequence of vertices
		for (std::map<double,int>::const_iterator it_tht=thetas.begin(); it_tht!=thetas.end(); it_tht++)
			edge_map[it_core->first].push_back(it_tht->second);
	}

	// perform a walk through the vertices: starting from each vertex, walking first to each of its adjacencies
	std::vector<std::vector<int> > all;

	int required_weinberg_length=0; // the Weinberg vector _must_ traverse each edge once in each direction!
	for (std::map<int,std::set<int> >::const_iterator it_adj=this->adjacency.begin(); it_adj!=this->adjacency.end(); it_adj++)
		required_weinberg_length+=it_adj->second.size();

	//for (std::map<int,flt2>::const_iterator it_to=this->positions.begin(); it_to!=this->positions.end(); it_to++) {
	//	std::map<int,std::set<int> >::const_iterator it_adj=this->adjacency.find(it_to->first);
	for (int mirror=0; mirror<2; mirror++) {
		for (std::map<int,std::set<int> >::const_iterator it_adj=this->adjacency.begin(); it_adj!=this->adjacency.end(); it_adj++) {
			for (std::set<int>::const_iterator it_vtx=it_adj->second.begin(); it_vtx!=it_adj->second.end(); it_vtx++) {
				std::vector<int> sequence;
				std::map<int,int> verts_visited;
				std::map<int,std::map<int,int> > edges_visited;

				int from = it_adj->first;
				int to = *it_vtx;

				sequence.push_back(from);
				verts_visited[from]++;
				edges_visited[from][to]++;
				weinberg_recurse(from, to, edge_map, sequence, verts_visited, edges_visited, mirror);

				if (int(sequence.size()) != required_weinberg_length) continue;

				relabel_weinberg(sequence);

				if (sequence[2]==sequence[0]) return std::string("0"); // lonely vertex

				all.push_back(sequence);
			}
		}
	}

	if (all.size()==0) return std::string("0");

	std::map<std::vector<int>,int> symmetries;
	for (unsigned int i = 0; i < all.size(); i++) {
		symmetries[all[i]]++;
	}

	//assert (symmetries.size() > 0);
	if (symmetries.size()==0) return std::string("0");

	wv_rotations = all.size();
	wv_symmetries = (symmetries.size()==0)?0:all.size() / symmetries.size(); // protect against divide-by-zero

	// Export the trial Weinberg vectors to ./grains/weinberg_XXXX.tttt.log
	/*
	std::stringstream ss;
	ss<<"grains/weinberg_"<<std::setfill('0')<<std::setw(5)<<std::right<<this->id<<'.'<<timestep<<".log";
	std::ofstream wvfile(ss.str().c_str());
	wvfile<<"Weinberg_Vector";
	for (int i=15; i<required_weinberg_length; i++) wvfile<<' ';
	wvfile<<"\tfrequency\n";
	for (std::map<std::vector<int>,int>::const_iterator it_wv=symmetries.begin(); it_wv!=symmetries.end(); it_wv++) {
		for (unsigned int i=0; i<it_wv->first.size(); i++) {
			char c=(it_wv->first[i]<27)?char(it_wv->first[i]+64):char(it_wv->first[i]+70);
			wvfile<<c;
		}
		wvfile<<'\t'<<it_wv->second<<'\n';
	}
	wvfile.close();
	*/

	weinberg_vector = symmetries.begin()->first;
	std::string wv(weinberg_vector.size(),'Z');
	for (unsigned int i=0; i<weinberg_vector.size(); i++)
		wv[i]=(weinberg_vector[i]<27)?char(weinberg_vector[i]+64):char(weinberg_vector[i]+70);

	return wv;
}

void Schlegel::weinberg_recurse(int vert_from, int vert_to, const std::map<int,std::vector<int> >& edge_map, std::vector<int> &sequence, std::map<int,int>& verts_visited, std::map<int,std::map<int,int> >& edges_visited, bool mirror)
{
	// BASE CASE
	if (sequence.size()==0 || verts_visited[vert_to] == int(adjacency[vert_to].size())) {
		return;
	}

	// SANITY CHECKS
	/*
	assert (verts_visited[vert_from] > 0);
	assert (verts_visited[vert_to] >= 0);
	assert (edges_visited[vert_from][vert_to] > 0);
	assert (verts_visited[vert_from] <= int(adjacency[vert_from].size()));
	assert (verts_visited[vert_to] <= int(adjacency[vert_to].size()));
	assert (edges_visited[vert_from][vert_to] <= 2);
	*/
	if (verts_visited[vert_from] <= 0 ||
			verts_visited[vert_to] < 0 ||
			edges_visited[vert_from][vert_to] <= 0 ||
			verts_visited[vert_from] > int(adjacency[vert_from].size()) ||
			verts_visited[vert_to] > int(adjacency[vert_to].size()) ||
			edges_visited[vert_from][vert_to] > 2
	) {
		sequence.clear();
		return;
	}

	// NORMAL CASE: stay on this face, "TURN RIGHT"
	//const Face& my_face = p->getFaceWithEdge(vert_from, vert_to);
	//int turn_right_next_index = (my_face.whichVert(vert_to) + 1) % my_face.numVertices();
	//int turn_right_vert_next = my_face[turn_right_next_index];
	std::map<int,std::vector<int> >::const_iterator to_vtx=edge_map.find(vert_to);
	int turn_from_index = std::find(to_vtx->second.begin(), to_vtx->second.end(), vert_from) - to_vtx->second.begin();
	int turn_right_next_index = (turn_from_index + 1) % to_vtx->second.size();
	int turn_right_vert_next = to_vtx->second[turn_right_next_index];

	// ALTERNATE, "TURN LEFT"
	//const Face& opposite_face = p->getFaceWithEdge(vert_to, vert_from);
	//int turn_left_next_index = (opposite_face.whichVert(vert_to) + opposite_face.numVertices() - 1) % opposite_face.numVertices();
	//int turn_left_vert_next = opposite_face[turn_left_next_index];
	int turn_left_next_index = (turn_from_index + to_vtx->second.size() - 1) % to_vtx->second.size();
	int turn_left_vert_next = to_vtx->second[turn_left_next_index];


	if (mirror) {
		int tmp = turn_right_vert_next;
		turn_right_vert_next = turn_left_vert_next;
		turn_left_vert_next = tmp;
	}



	// if we haven't been to this vertex before, simply go around the face (turn right)
	if (verts_visited[vert_to] == 0) {
		sequence.push_back(vert_to);
		verts_visited[vert_to]++;
		edges_visited[vert_to][turn_right_vert_next]++;
		weinberg_recurse(vert_to, turn_right_vert_next, edge_map, sequence, verts_visited, edges_visited, mirror);
	} else {
		// assert (verts_visited[vert_to] == 1 || verts_visited[vert_to] == 2); // assumes 3-valent graph
		// if we haven't already traversed the backward edge, take that route (backwards)
		if (edges_visited[vert_to][vert_from] == 0) {
			sequence.push_back(vert_to);
			verts_visited[vert_to]++;
			edges_visited[vert_to][vert_from]++;
			weinberg_recurse(vert_to, vert_from, edge_map, sequence, verts_visited, edges_visited, mirror);
		}	else {
			// This assertion fails when 2-edged faces are present.
			//assert (edges_visited[std::make_pair(vert_to, vert_from)] == 1);

			// if turn right is available, do that (turn right)
			if (edges_visited[vert_to][turn_right_vert_next] == 0) {
				sequence.push_back(vert_to);
				verts_visited[vert_to]++;
				edges_visited[vert_to][turn_right_vert_next]++;
				weinberg_recurse(vert_to, turn_right_vert_next, edge_map, sequence, verts_visited, edges_visited, mirror);
			} else { // otherwise, "turn left" (if 3-valent; otherwise, cycle CCW)
				// This assertion fails when 2-edged faces are present.
				// assert (edges_visited[std::make_pair(vert_to, turn_right_vert_next)] == 1);
				//sequence.push_back(vert_to);
				//verts_visited[vert_to]++;
				//edges_visited[vert_to][turn_left_vert_next]++;
				//weinberg_recurse(vert_to, turn_left_vert_next, edge_map, sequence, verts_visited, edges_visited, mirror);

				// don't turn left, increment right until it becomes available
				if (!mirror) {
					int skip_right_next_index = (turn_right_next_index+1) % to_vtx->second.size();
					int skip_right_vert_next = to_vtx->second[skip_right_next_index];
					while ((edges_visited[vert_to][skip_right_vert_next] != 0) && (skip_right_next_index!=turn_right_next_index)) {
						skip_right_next_index = (skip_right_next_index+1) % to_vtx->second.size();
						skip_right_vert_next = to_vtx->second[skip_right_next_index];
					}
					sequence.push_back(vert_to);
					verts_visited[vert_to]++;
					edges_visited[vert_to][skip_right_vert_next]++;
					weinberg_recurse(vert_to, skip_right_vert_next, edge_map, sequence, verts_visited, edges_visited, mirror);
				} else {
					int skip_left_next_index = (turn_left_next_index + to_vtx->second.size() - 1) % to_vtx->second.size();
					int skip_left_vert_next = to_vtx->second[skip_left_next_index];
					while ((edges_visited[vert_to][skip_left_vert_next] != 0) && (skip_left_next_index!=turn_left_next_index)) {
						skip_left_next_index = (skip_left_next_index + to_vtx->second.size() - 1) % to_vtx->second.size();
						skip_left_vert_next = to_vtx->second[skip_left_next_index];
					}
					sequence.push_back(vert_to);
					verts_visited[vert_to]++;
					edges_visited[vert_to][skip_left_vert_next]++;
					weinberg_recurse(vert_to, skip_left_vert_next, edge_map, sequence, verts_visited, edges_visited, mirror);
				}
			}
		}
	}
}

//void Schlegel::export_eps(const std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn_edg, const int x, const std::string& timestep)
void Schlegel::export_eps(const std::map<grainid,std::set<std::set<int> > >::const_iterator it_grn_edg, const int x, const std::string& wv)
{
	// Export the result as an EPS file named ./grains/schlegel-XXXX.eps
	std::stringstream ss;
	//ss<<"grains/grain_"<<std::setfill('0')<<std::setw(5)<<std::right<<it_grn_edg->first<<'.'<<timestep<<".eps";
	ss<<"../grains/schlegel-"<<wv<<".eps";
	// Check whether the file exists. If it doesn't, or if it does and this topology is complete (x==2), write to disk.
	if (!fileExists(ss.str()) || x==2) {
		std::ofstream imgfile(ss.str().c_str());
		imgfile<<"%!PS-Adobe-3.0 EPSF-3.0\n"
	    	   <<"%%BoundingBox: "<<int(-r-2)<<' '<<int(-r-2)<<' '<<int(r+2)<<' '<<int(r+2)<<'\n'
	  	     <<"1.0 setlinewidth\n";

		if (x==2) imgfile<<"0  setgray\n\n";
		else imgfile<<"0.5 setgray\n\n"; // Euler characteristic violated: this grain is suspect! Shade it funny colors.

		for (std::map<int,flt2>::iterator it_pos=this->positions.begin(); it_pos!=this->positions.end(); it_pos++) {
			for (std::set<int>::const_iterator it_vtx=this->adjacency[it_pos->first].begin(); it_vtx!=this->adjacency[it_pos->first].end(); it_vtx++) {
				if (*it_vtx<=it_pos->first) continue; // skip vertices with smaller ID: edge already drawn
				imgfile<< it_pos->second.x <<'\t'<< it_pos->second.y <<"\tnewpath moveto\t"<< this->positions[*it_vtx].x <<'\t'<< this->positions[*it_vtx].y<<"\tlineto closepath stroke\n";
			}
		}
		imgfile<<"showpage\n";
		imgfile.close();
	}
	ss.str("");
}

#endif // _SCHLEGEL_HPP_
