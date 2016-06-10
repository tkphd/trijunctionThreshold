// File:    extract_volume.cpp
// Purpose: reads MMSP grid containing sparse vector of order parameters
//          writes volume of specified grain in the specified data file

// Questions/Comments to kellet@rpi.edu (Trevor Keller)


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

int main(int argc, char* argv[])
{
  if ( argc != 3 )
  {
    std::cout << "Usage: " << argv[0] << " input.dat grainid\n";
    std::cout << "   eg, " << argv[0] << " BGQ.00100.dat 11\n";
    exit(-1);
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
  bool sparse_type = (type.find("sparse") != std::string::npos);

  if (not float_type   and 	not sparse_type){
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
	  if (float_type) {
	  	const int gid = std::atoi(argv[2]);
	  	if (dim==2) {
		  	MMSP::grid<2,MMSP::sparse<float> > grid(argv[1]);
		  	double pxarea=1.0;
		  	for (int d=0; d<dim; d++) pxarea*=MMSP::dx(grid,d);

		  	double area=0.0;
		  	for (int n=0; n<nodes(grid); n++){
		  		MMSP::vector<int> x=MMSP::position(grid, n);
		  		area += grid(x)[gid];
		  	}
				std::cout<<area*pxarea<<'\n';
	  	} else if (dim==3) {
		  	MMSP::grid<3,MMSP::sparse<float> > grid(argv[1]);
		  	double vxvol=1.0;
		  	for (int d=0; d<dim; d++) vxvol*=MMSP::dx(grid,d);

		  	double volume=0.0;
		  	for (int n=0; n<nodes(grid); n++){
		  		MMSP::vector<int> x=MMSP::position(grid, n);
		  		volume += grid(x)[gid];
		  	}
				std::cout<<volume*vxvol<<'\n';
		  } else {
		    std::cerr<<"File input error: "<<dim<<"-D not implemented."<<std::endl;
    		exit(-1);
		  }
		}
  } else {
    std::cerr<<"File input error: "<<type<<" not implemented."<<std::endl;
    exit(-1);
  }

  return 0;
}
