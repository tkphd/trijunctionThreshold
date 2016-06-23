// mmsp2png.cpp
// Convert MMSP grid data to grayscale PNG image format
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include<iostream>
#include<cstdlib>
#include<string>
#include<sstream>
#include<fstream>
#include<vector>
#include<set>
#include<cmath>
#include<png.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_multifit_nlin.h>

#include"MMSP.hpp"
#include"graingrowth.hpp"
#include"curvefit.hpp"

double trijunction_profile(const double& a, const double& theta, const double& x)
{
	// y = a*acos(exp(-2*theta*x/a)*sin(theta))/(2*theta) - a*(-theta + pi/2)/(2*theta)
	return a*acos(exp(-2.0*theta*x/a)*sin(theta))/(2.0*theta) - a*(M_PI/2.0 - theta)/(2.0*theta);
}

int writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename);

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,T>& GRID,
                                                    const int& mode, int sliceaxis, int slicelevel, const double& zoomin, const double& zoomax, const bool coninv,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const unsigned int& bufsize, unsigned char* buffer);

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID,
                                                    const int& mode, int sliceaxis, int slicelevel, const double& zoomin, const double& zoomax, const bool coninv,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer);

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
                                                    const int& mode, int sliceaxis, int slicelevel, const double& zoomin, const double& zoomax, const bool coninv,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer);

// enumerate modes to avoid confusion
enum {
	Magnitude = 1,
	InclField = 2,
	ExclField = 3,
	Isolines  = 4,
	Mobility  = 5,
	CurveFit  = 6
};

int main(int argc, char* argv[])
{
	// command line error check
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " [--help] [--slice] [--zoom] [--contour] [--mag|--field|--exclude] infile [outfile]\n";
		std::exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help") {
		std::cout << argv[0] << ": convert MMSP grid data to PNG image format.\n";
		std::cout << "Usage:    " << argv[0] << " [--help] [--slice=X,n] [--zoom] [--contour=n,a,b] [--mag|--field=n|--exclude=n] [--invert] infile [outfile]\n\n";
		std::cout << "Examples: " << argv[0] << " --help\n"
		          << "             displays this message.\n";
		std::cout << "          " << argv[0] << " infile\n"
		          << "             writes magnitude of infile values to infile.png, bounded on [0,1].\n"
		          << "             If infile is 3-D, holds X-coord constant at the midpoint.\n";
		std::cout << "          " << argv[0] << " --mob infile\n"
		          << "             writes mobility associated with magnitude of infile values to infile.png\n";
		std::cout << "          " << argv[0] << " --zoom infile\n"
		          << "             writes magnitude of infile values to infile.png, bounded by data min and max exactly.\n"
		          << "             Use in conjunction with any other flag for data with very small values.\n";
		std::cout << "          " << argv[0] << " --zoom=a,b infile\n"
		          << "             writes magnitude of infile values to infile.png, bounded on [a,b].\n"
		          << "             Use in conjunction with any other flag for data with very small values.\n";
		std::cout << "          " << argv[0] << " --slice=X,32 infile outfile\n"
		          << "             writes magnitude of infile values with constant X-coord of 32 to outfile, bounded on [0,1].\n"
		          << "             Use for 3-D data in conjunction with any other flag. Y and Z also work.\n"
		          << "             Defaults to the mid-plane parallel to the specified axis.\n";
		std::cout << "          " << argv[0] << " infile.dat\n"
		          << "             writes magnitude of infile.dat values to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --mag infile.dat \n"
		          << "             writes magnitude of infile.dat values to infile.png, bounded on [0,1].\n"
		          << "             Use with only one field (--field=n or --exclude all others) to visualize its absolute value.\n";
		std::cout << "          " << argv[0] << " --field=0 infile.dat\n"
		          << "             writes field 0 of infile.dat values, only, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --field=1,2,3 infile.dat\n"
		          << "             writes magnitude of fields 1,2,3 of infile.dat values, only, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --exclude=0 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --exclude=0 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except field 0, to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --exclude=1,2,3 infile.dat\n"
		          << "             writes magnitude of infile.dat values, except fields 1,2,3,to infile.png, bounded on [0,1].\n";
		std::cout << "          " << argv[0] << " --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             makes black pixels where field 0 equals 0.25,0.50,0.75. Use with any other flag.\n";
		std::cout << "          " << argv[0] << " --contour=-1,0.25,0.50,0.75 infile.dat\n"
		          << "             makes black pixels where magnitude equals 0.25,0.50,0.75, ignoring the field index.\n";
		std::cout << "          " << argv[0] << " --coninv --contour=0,0.25,0.50,0.75 infile.dat\n"
		          << "             makes white pixels where field 0 equals 0.25,0.50,0.75. Use with any other flag.\n";
		std::cout << "          " << argv[0] << " --invert infile.dat\n"
		          << "             writes magnitude of infile.dat values to infile.png, bounded on [255,0].\n"
		          << "             Use in conjunction with any other flag: inversion is the last operation before output.\n";
		std::cout << "Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n";
		std::exit(0);
	}

	int datindex = 0; // in typical usage, input filename comes immediately after executable
	int pngindex = 0; // in typical usage, output filename comes immediately after input
	int mode = 0;
	bool invert=false; // filter to invert the whole image before writing
	bool coninv=false; // filter to set contour pixels white instead of black
	double contol = 0.06; // contour tolerance (5-7% works OK)
	int sliceaxis = -1;
	int slicelevel = -1;
	double zoomin=0.0, zoomax=1.0;
	int lvlfield = -1; // store the field on which contour are calculated
	std::set<int> fieldset; // for --field or --exclude
	std::set<double> levelset; // for contours

	// Parse command line flags
	for (int i=1; i<argc; i++) {
		std::string flag(argv[i]);
		if (flag == "--mag") {
			mode=std::max(int(Magnitude),mode);
		} else if (flag == "--invert") {
			invert=true;
		} else if (flag.substr(0,6) == "--zoom") {
			// Interpolate grayscale between data minimum and maximum, instead of [0,1]
			if (flag.length()==6 || flag[6]!='=') {
				zoomin = std::sqrt(std::numeric_limits<double>::max());
				zoomax = -zoomin;
			} else {
				std::string val = flag.substr(7,255);
				char* cstr = new char [val.length()+1];
				std::strcpy(cstr, val.c_str());
				char* ptr = strtok(cstr,",");
				zoomin = atof(ptr);
				ptr = strtok(NULL,",");
				zoomax = atof(ptr);
				ptr = strtok(NULL,",");
				assert(ptr==NULL);
				delete [] cstr;
			}
		} else if (flag.substr(0,7) == "--slice") {
			// Create image from 2D slice of 3D dataset
			if (flag.length()==7 || flag[7]!='=') {
				std::cerr<<"Error: expected --slice=X[|Y|Z][,n].\n"<<std::endl;
				std::exit(-1);
			}
			switch(flag[8]) {
				case 'x': sliceaxis=0;
				          break;
				case 'y': sliceaxis=1;
				          break;
				case 'z': sliceaxis=2;
				          break;
				case 'X': sliceaxis=0;
				          break;
				case 'Y': sliceaxis=1;
				          break;
				case 'Z': sliceaxis=2;
				          break;
			}
			if (sliceaxis==-1) {
				std::cerr<<"Error: axis "<<flag[8]<<" not recognized.\n"<<std::endl;
				std::exit(-1);
			}
			if (flag.length()>10)
				slicelevel = atoi(flag.substr(10,20).c_str());
			else
				std::cout<<"Slicing through the midpoint of "<<flag[8]<<" axis."<<std::endl;
		} else if (flag.substr(0,8) == "--coninv") {
			// Overlay white contours instead of black
			coninv = true;
		} else if (flag.substr(0,8) == "--contol") {
			// Set tolerance or contour line width
			if (flag.length()==8 || flag[8]!='=') {
				std::cerr<<"Error: expected --contol=a\n"<<std::endl;
				std::exit(-1);
			}
			std::string val = flag.substr(10,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			contol=(atof(ptr));
			ptr = strtok(NULL,",");
			assert(ptr==NULL);
			delete [] cstr;
		} else if (flag.substr(0,9) == "--contour") {
			// Overlay contours of specified field at specified value(s)
			if (flag.length()==9 || flag[9]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=field,a[,b,...]\n"<<std::endl;
				std::exit(-1);
			}
			std::string val = flag.substr(10,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			lvlfield=(atoi(ptr));
			ptr = strtok(NULL,",");
			while (ptr != NULL) {
				levelset.insert(atof(ptr));
				ptr = strtok(NULL,",");
			}
			delete [] cstr;
			if (lvlfield<-1) {
				std::cerr<<"Error: field "<<lvlfield<<" invalid.\n"<<std::endl;
				std::exit(-1);
			}
		} else if (flag.substr(0,7) == "--field") {
			// Render only the specified field(s)
			if (flag[7]!='=') {
				std::cerr<<"Error: expected "<<flag<<"=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(int(InclField),mode);
 			std::string val = flag.substr(8,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			while (ptr != NULL) {
				fieldset.insert(atoi(ptr));
				ptr = strtok(NULL,",");
			}
			delete [] cstr;
		} else if (flag.substr(0,9) == "--exclude") {
			// Render everything but the specified field(s)
			if (flag.length()==9 || flag[9]!='=') {
				std::cerr<<"Error: expected --exclude=n[,p,q,...]\n"<<std::endl;
				std::exit(-1);
			}
			mode=std::max(int(ExclField),mode);
			std::string val = flag.substr(10,255);
			char* cstr = new char [val.length()+1];
			std::strcpy(cstr, val.c_str());
			char* ptr = strtok(cstr,",");
			while (ptr != NULL) {
				fieldset.insert(atoi(ptr));
				ptr = strtok(NULL,",");
			}
			delete [] cstr;
		} else if (flag.substr(0,5) == "--fit") {
			// Fit a tricrystal grain boundary profile to the analytical expression
			mode=std::max(int(CurveFit),mode);
		} else if (flag.substr(0,5) == "--mob") {
			// Plot mobility instead of microstructure
			mode=std::max(int(Mobility),mode);
		} else if (flag.find(".png") != std::string::npos) {
			// if a PNG filename is provided, use it as the image destination
			pngindex = i;
		} else {
			// MMSP data input file
			datindex = i;
		}
	}

// file open error check
	if (datindex==0) {
		std::cerr << "File input error: no input specified."<<".\n";
		exit(-1);
	}
	std::ifstream input(argv[datindex]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[datindex] << ".\n";
		exit(-1);
	}

// generate output file name
	std::stringstream pngname;
	if (pngindex==0) {
		std::string datname(argv[datindex]);
		int extpos = datname.find_last_of(".");
		if (datname.find_first_of("0123456789",extpos) != std::string::npos)
			pngname << datname << ".png";
		else
			pngname << datname.substr(0, extpos) << ".png";
	} else
		pngname << argv[pngindex];

// read data type
	std::string type;
	getline(input, type, '\n');

// grid type error check
	if (type.substr(0, 4) != "grid") {
		std::cerr << "File input error: file does not contain grid data." << std::endl;
		exit(-1);
	}

// parse data type
	bool bool_type = (type.find("bool") != std::string::npos);
	bool char_type = (type.find("char") != std::string::npos);
	bool unsigned_char_type = (type.find("unsigned char") != std::string::npos);
	bool int_type = (type.find("int") != std::string::npos);
	bool unsigned_int_type = (type.find("unsigned int") != std::string::npos);
	bool long_type = (type.find("long") != std::string::npos);
	bool unsigned_long_type = (type.find("unsigned long") != std::string::npos);
	bool short_type = (type.find("short") != std::string::npos);
	bool unsigned_short_type = (type.find("unsigned short") != std::string::npos);
	bool float_type = (type.find("float") != std::string::npos);
	bool double_type = (type.find("double") != std::string::npos);
	bool long_double_type = (type.find("long double") != std::string::npos);

	bool scalar_type = (type.find("scalar") != std::string::npos);
	bool vector_type = (type.find("vector") != std::string::npos);
	bool sparse_type = (type.find("sparse") != std::string::npos);

	if (not bool_type    and
	        not char_type    and  not unsigned_char_type   and
	        not int_type     and  not unsigned_int_type    and
	        not long_type    and  not unsigned_long_type   and
	        not short_type   and  not unsigned_short_type  and
	        not float_type   and
	        not double_type  and  not long_double_type) {
		std::cerr << "File input error: unknown grid data type." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input >> dim;

	// read number of fields
	int fields;
	input >> fields;

	// read grid sizes
	int x0[3] = {0, 0, 0};
	int x1[3] = {0, 0, 0};
	for (int d=0; d<dim; d++)
		input >> x0[d] >> x1[d];

	// read cell spacing
	float dx[3] = {1.0, 1.0, 1.0};
	for (int d=0; d<dim; d++)
		input >> dx[d];

	// ignore trailing endlines
	input.ignore(10, '\n');

	// if slice was not specified, handle it
	if (dim==3 && sliceaxis<0) {
		sliceaxis=2;
		if (slicelevel<x0[sliceaxis])
			slicelevel = (x1[sliceaxis] - x0[sliceaxis])/2;
	}

	int image_size[2] = {1,1};
	int i=0;
	for (int d=0; d<dim; d++) {
		if (d==sliceaxis)
			continue;
		image_size[i] = x1[d]-x0[d];
		i++;
	}
	unsigned int theSize=image_size[0] * image_size[1];
	unsigned char* buffer = new unsigned char[theSize];

	// write grid data
	if (scalar_type or (not vector_type and not sparse_type))
	{ // must be scalar or built-in
		/*
		if (bool_type) {
			if (dim == 1) {
				MMSP::grid<1,bool> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,bool> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,bool> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (unsigned_char_type) {
			if (dim == 1) {
				MMSP::grid<1,unsigned char> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,unsigned char> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,unsigned char> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (char_type) {
			if (dim == 1) {
				MMSP::grid<1,char> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,char> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,char> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (unsigned_int_type) {
			if (dim == 1) {
				MMSP::grid<1,unsigned int> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,unsigned int> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,unsigned int> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (int_type) {
			if (dim == 1) {
				MMSP::grid<1,int> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,int> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,int> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (unsigned_long_type) {
			if (dim == 1) {
				MMSP::grid<1,unsigned long> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,unsigned long> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,unsigned long> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (long_type) {
			if (dim == 1) {
				MMSP::grid<1,long> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,long> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,long> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (unsigned_short_type) {
			if (dim == 1) {
				MMSP::grid<1,unsigned short> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,unsigned short> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,unsigned short> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (short_type) {
			if (dim == 1) {
				MMSP::grid<1,short> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,short> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,short> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (float_type) {
			if (dim == 1) {
				MMSP::grid<1,float> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,float> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,float> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (long_double_type) {
			if (dim == 1) {
				MMSP::grid<1,long double> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,long double> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,long double> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		else if (double_type) {
			if (dim == 1) {
				MMSP::grid<1,double> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,double> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,double> GRID(argv[datindex]);
				convert_scalars(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, theSize, buffer);
			}
		}
		*/
	}

	else if (vector_type)
	{
		/*
		if (bool_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<bool> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<bool> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<bool> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_char_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<unsigned char> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<unsigned char> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<unsigned char> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (char_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<char> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<char> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<char> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_int_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<unsigned int> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<unsigned int> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<unsigned int> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (int_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<int> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<int> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<int> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_long_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<unsigned long> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<unsigned long> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<unsigned long> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (long_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<long> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<long> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<long> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_short_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<unsigned short> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<unsigned short> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<unsigned short> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (short_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<short> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<short> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<short> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (float_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<float> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<float> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<float> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (long_double_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<long double> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<long double> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<long double> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (double_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::vector<double> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::vector<double> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::vector<double> > GRID(argv[datindex]);
				convert_vectors(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		*/
	}
	else if (sparse_type)
	{
		if (bool_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<bool> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<bool> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<bool> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_char_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<unsigned char> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (char_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<char> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<char> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<char> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_int_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<unsigned int> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (int_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<int> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<int> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<int> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_long_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<unsigned long> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (long_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<long> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<long> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<long> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (unsigned_short_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<unsigned short> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (short_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<short> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<short> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<short> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (float_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<float> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<float> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<float> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (long_double_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<long double> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<long double> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<long double> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
		else if (double_type) {
			if (dim == 1) {
				MMSP::grid<1,MMSP::sparse<double> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 2) {
				MMSP::grid<2,MMSP::sparse<double> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			} else if (dim == 3) {
				MMSP::grid<3,MMSP::sparse<double> > GRID(argv[datindex]);
				convert_sparses(GRID, mode, sliceaxis, slicelevel, zoomin, zoomax, coninv, levelset, lvlfield, contol, fieldset, theSize, buffer);
			}
		}
	}
	if (invert)
		for (unsigned int n=0; n<theSize; n++)
			buffer[n] = 255-buffer[n];

	int result = writePNG(image_size[0], image_size[1], 1,buffer, pngname.str().c_str());

	// clean up
	delete [] buffer;

	return result;
}

int writePNG(const int w, const int h, const int bpp, unsigned char* imData, const char* filename)
{
	// using libpng
	// After "A simple libpng example program,"
	// http://zarb.org/~gc/html/libpng.html
	// and the libpng manual, http://www.libpng.org/pub/png

	png_byte color_type = PNG_COLOR_TYPE_GRAY;
	// valid choices: PNG_COLOR_TYPE_GRAY       (bit depths 1,2,4, 8, 16)
	//                PNG_COLOR_TYPE_GRAY_ALPHA (bit depths 8, 16)
	//                PNG_COLOR_TYPE_PALETTE    (bit depths 1,2,4, 8)
	//                PNG_COLOR_TYPE_RGB        (bit_depths 8, 16)
	//                PNG_COLOR_TYPE_RGB_ALPHA  (bit_depths 8, 16)
	//                PNG_COLOR_MASK_PALETTE
	//                PNG_COLOR_MASK_COLOR
	//                PNG_COLOR_MASK_ALPHA

	png_byte bit_depth = 8; // valid choices: 1,2,4, 8, 16
	png_structp png_ptr;
	png_infop info_ptr;

	png_bytepp row_pointers = new png_bytep[h];
	for (int j=0; j<h; j++)
		row_pointers[j] = &imData[j*w];

	// Setup PNG file
	FILE *fp = fopen(filename, "wb");
	if (!fp) {
		std::cerr<<"Error making image: check permissions on "<<filename<<std::endl;
		return (-1);
	}
	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) {
		std::cerr<<"Error making image: png_create_write_struct failed."<<std::endl;
		return (-1);
	}
	info_ptr = png_create_info_struct(png_ptr);
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to init_io."<<std::endl;
		return (-1);
	}
	png_init_io(png_ptr, fp);

	// Write PNG header
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to write header."<<std::endl;
		return (-1);
	}
	png_set_IHDR(png_ptr, info_ptr, w, h,
	                 bit_depth, color_type, PNG_INTERLACE_NONE,
	                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

	png_write_info(png_ptr, info_ptr);

	// Write image
	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to write data."<<std::endl;
		return (-1);
	}
	png_write_image(png_ptr, row_pointers);

	if (setjmp(png_jmpbuf(png_ptr))) {
		std::cerr<<"Error making image: unable to finish writing."<<std::endl;
		return (-1);
	}
	png_write_end(png_ptr, NULL);

	// Clean up
	delete [] row_pointers;

	fclose(fp);

	return 0;
}

template <int dim, typename T> void convert_scalars(const MMSP::grid<dim,T>& GRID,
                                                    const int& mode, int sliceaxis, int slicelevel, const double& zoomin, const double& zoomax, const bool coninv,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const unsigned int& bufsize, unsigned char* buffer)
{
	double min=zoomin;
	double max=zoomax;

	for (int n=0; n<MMSP::nodes(GRID); n++) {
		T val = GRID(n);
		if (mode==Magnitude) // mag
			val = std::abs(val);
		if (val>max)
			max=val;
		else if (val<min)
			min=val;
	}

	std::cout<<"Rescaling on ["<<min<<','<<max<<"]."<<std::endl;

	if (dim==1) {
		unsigned int n=0;
		MMSP::vector<int> x(1,0);
		for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
			T val = GRID(x);
			if (mode==Magnitude) //mag
				val = std::abs(val);
			assert(n<bufsize);
			buffer[n] = 255*((val-min)/(max-min));
			if (levelset.size()>0) //contour
				for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
					if (std::fabs(val-*it)/std::fabs(*it)<contol)
						buffer[n] = 255-buffer[n];
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		MMSP::vector<int> x(2,0);
		for (x[1] = MMSP::g1(GRID,1)-1; x[1] >= MMSP::g0(GRID,1); x[1]--)
			for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
				T val = GRID(x);
				if (mode==Magnitude) //mag
					val = std::abs(val);
				assert(n<bufsize);
				buffer[n] = 255*((val-min)/(max-min));
				if (levelset.size()>0) //contour
					for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
						if (std::fabs(val-*it)/std::fabs(*it)<contol)
							buffer[n] = 255-buffer[n];
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		MMSP::vector<int> x(3,0);
		for (x[2] = MMSP::g0(GRID,2); x[2] < MMSP::g1(GRID,2); x[2]++)
			for (x[1] = MMSP::g1(GRID,1)-1; x[1] >= MMSP::g0(GRID,1); x[1]--)
				for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					T val = GRID(x);
					if (mode==Magnitude) //mag
						val = std::abs(val);
					assert(n<bufsize);
					buffer[n] = 255*((val-min)/(max-min));
					if (levelset.size()>0) //contour
						for (std::set<double>::iterator it=levelset.begin(); it!=levelset.end(); it++)
							if (std::fabs(val-*it)/std::fabs(*it)<contol)
								buffer[n] = coninv ? 255 : 0;
					n++;
				}
	}
}

template <int dim, typename T> void convert_vectors(const MMSP::grid<dim,MMSP::vector<T> >& GRID,
                                                    const int& mode, int sliceaxis, int slicelevel, const double& zoomin, const double& zoomax, const bool coninv,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer)
{
	double min=zoomin;
	double max=zoomax;
	int included = (mode==InclField) ? fieldset.size() : MMSP::fields(GRID)-fieldset.size();

	for (int n=0; n<MMSP::nodes(GRID); n++) {
		double sum=0.0;
		if (mode==0) {        //  no option specified
			for (int i=0; i<MMSP::fields(GRID); i++)
				if (included>1)
					sum += pow(GRID(n)[i],2.0);
				else
					sum = GRID(n)[i];
		} else if (mode==Magnitude) { //  --mag
			for (int i=0; i<MMSP::fields(GRID); i++)
				sum += pow(GRID(n)[i],2.0);
		} else if (mode==InclField) { //  --field
			if (included>1)
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
					sum += pow(GRID(n)[*it],2.0);
			else
				sum = GRID(n)[*fieldset.begin()];
		} else if (mode==ExclField) { //  --exclude
			for (int i=0; i<MMSP::fields(GRID); i++) {
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end()) {
					if (included>1)
						sum += pow(GRID(n)[i],2.0);
					else
						sum = GRID(n)[i];
				}
			}
		}
		if (mode==Magnitude || included!=1)
			sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}

	std::cout<<"Rescaling on ["<<min<<','<<max<<"]."<<std::endl;

	if (dim==1) {
		unsigned int n=0;
		MMSP::vector<int> x(1,0);
		for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
			double sum=0.0;
			if (mode==0) { //         default selection
				for (int i=0; i<MMSP::fields(GRID); i++) {
					if (included>1)
						sum += pow(GRID(x)[i],2.0);
					else
						sum = GRID(x)[i];
				}
			} else if (mode==Magnitude) { //  --mag
				for (int i=0; i<MMSP::fields(GRID); i++)
					sum += pow(GRID(x)[i],2.0);
			} else if (mode==InclField) { //  --field
				if (fieldset.size()==1)
					sum = GRID(x)[*fieldset.begin()];
				else
					for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
						sum += pow(GRID(x)[*it],2.0);
			} else if (mode==ExclField) { //  --exclude
				for (int i=0; i<MMSP::fields(GRID); i++) {
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end()) {
						if (included>1)
							sum += pow(GRID(x)[i],2.0);
						else
							sum = GRID(x)[i];
					}
				}
			}
			if (mode==Magnitude || included!=1)
				sum = std::sqrt(sum);
			assert(n<bufsize);
			buffer[n] = 255*((sum-min)/(max-min));
			if (levelset.size()>0) {// --contour
				double value = (lvlfield>-1) ? GRID(x)[lvlfield] : sum;
				for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
					if (std::fabs(value - *itl)/std::fabs(*itl)<contol)
						buffer[n] = coninv ? 255 : 0;
			}
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		MMSP::vector<int> x(2,0);
		for (x[1] = MMSP::g1(GRID,1)-1; x[1] >= MMSP::g0(GRID,1); x[1]--)
			for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
				double sum=0.0;
				if (mode==0) { //         default selection
					for (int i=0; i<MMSP::fields(GRID); i++) {
						if (included>1)
							sum += pow(GRID(x)[i],2.0);
						else
							sum = GRID(x)[i];
					}
				} else if (mode==Magnitude) { //  --mag
					for (int i=0; i<MMSP::fields(GRID); i++)
						sum += pow(GRID(x)[i],2.0);
				} else if (mode==InclField) { //  --field
					if (fieldset.size()==1) {
						sum = GRID(x)[*fieldset.begin()];
					} else {
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
							sum += pow(GRID(x)[*it],2.0);
					}
				} else if (mode==ExclField) { //  --exclude
					for (int i=0; i<MMSP::fields(GRID); i++) {
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end()) {
							if (included>1)
								sum += pow(GRID(x)[i],2.0);
							else
								sum = GRID(x)[i];
						}
					}
				}
				if (mode==Magnitude || included!=1)
					sum = std::sqrt(sum);
				assert(n<bufsize);
				buffer[n] = 255*((sum-min)/(max-min));
				if (levelset.size()>0) {// --contour
					double value = (lvlfield>-1) ? GRID(x)[lvlfield] : sum;
					for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
						if (std::fabs(value - *itl)/std::fabs(*itl)<contol)
							buffer[n] = coninv ? 255 : 0;
				}
				n++;
			}
	} else if (dim==3) {
		unsigned int n=0;
		MMSP::vector<int> x(3,0);
		for (x[2] = MMSP::g0(GRID,2); x[2] < MMSP::g1(GRID,2); x[2]++)
			for (x[1] = MMSP::g1(GRID,1)-1; x[1] >= MMSP::g0(GRID,1); x[1]--)
				for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double sum=0.0;
					if (mode==0) { //         default selection
						for (int i=0; i<MMSP::fields(GRID); i++) {
							if (included>1)
								sum += pow(GRID(x)[i],2.0);
							else
								sum = GRID(x)[i];
						}
					} else if (mode==Magnitude) { //  --mag
						for (int i=0; i<MMSP::fields(GRID); i++)
							sum += pow(GRID(x)[i],2.0);
					} else if (mode==InclField) { //  --field
						if (fieldset.size()==1)
							sum = GRID(x)[*fieldset.begin()];
						else
							for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++)
								sum += pow(GRID(x)[*it],2.0);
					} else if (mode==ExclField) { //  --exclude
						for (int i=0; i<MMSP::fields(GRID); i++) {
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end()) {
								if (included>1)
									sum += pow(GRID(x)[i],2.0);
								else
									sum = GRID(x)[i];
							}
						}
					}
					if (mode==Magnitude || included!=1)
						sum = std::sqrt(sum);
					assert(n<bufsize);
					buffer[n] = 255*((sum-min)/(max-min));
					if (levelset.size()>0) {// --contour
						double value = (lvlfield>-1) ? GRID(x)[lvlfield] : sum;
						for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
							if (std::fabs(value - *itl)/std::fabs(*itl)<contol)
								buffer[n] = coninv ? 255 : 0;
					}
					n++;
				}
	}
}

template <int dim, typename T> void convert_sparses(const MMSP::grid<dim,MMSP::sparse<T> >& GRID,
                                                    const int& mode, int sliceaxis, int slicelevel, const double& zoomin, const double& zoomax, const bool coninv,
                                                    const std::set<double>& levelset, const int& lvlfield, const double& contol, const std::set<int>& fieldset,
                                                    const unsigned int& bufsize, unsigned char* buffer)
{
	double min=zoomin;
	double max=zoomax;
	int included = fieldset.size();
	const T mu_lo = 0.01;
	const T mu_hi = 1.00;
	const T mu_x = 0.6422; // mean of mobility Gaussian
	const T mu_s = 0.0175; // deviation of mobility Gaussian

	double minPhi = 10.0;
	MMSP::vector<int> tj(dim,0);

	for (int n=0; n<MMSP::nodes(GRID); n++) {
		double sum=0.0;
		if (mode==InclField) { //  --field
			for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
				if (included>1)
					sum += pow(GRID(n)[*it],2.0);
				else
					sum = GRID(n)[*it];
			}
		} else if (mode==ExclField) { //  --exclude
			for (int h=0; h<GRID(n).length(); h++) {
				int i=GRID(n).index(i);
				std::set<int>::iterator it=fieldset.find(i);
				if (it == fieldset.end())
					sum += pow(GRID(n).value(h),2.0);
			}
		} else if (mode==CurveFit) { //  --fit
			for (int h=0; h<GRID(n).length(); h++)
				sum += pow(GRID(n).value(h) - 0.5,2.0);
		} else if (mode==Mobility) {
			sum = mobility(mu_lo, mu_hi, mu_x, mu_s, GRID(n).getMagPhi());
		} else { //         --mag or default selection
			for (int h=0; h<GRID(n).length(); h++)
				sum += pow(GRID(n).value(h),2.0);
		}
		if (GRID(n).getMagPhi()<minPhi) {
			tj = MMSP::position(GRID, n);
			minPhi = GRID(tj).getMagPhi();
		}

		if (included!=1 && mode!=Mobility)
			sum = std::sqrt(sum);
		if (sum>max)
			max=sum;
		else if (sum<min)
			min=sum;
	}

	std::cout<<"Rescaling on ["<<min<<','<<max<<"]."<<std::endl;

	if (dim==1) {
		unsigned int n=0;
		MMSP::vector<int> x(1,0);
		for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
			double sum=0.0;
			if (mode==InclField) { //  --field
				for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
					if (included>1)
						sum += pow(GRID(x)[*it],2.0);
					else
						sum = GRID(x)[*it];
				}
			} else if (mode==ExclField) { //  --exclude
				for (int h=0; h<GRID(x).length(); h++) {
					int i=GRID(x).index(i);
					std::set<int>::iterator it=fieldset.find(i);
					if (it == fieldset.end())
						sum += pow(GRID(x).value(h),2.0);
				}
			} else if (mode==Mobility) {
				sum = mobility(mu_lo, mu_hi, mu_x, mu_s, GRID(x).getMagPhi());
			} else { //          --mag or default
				for (int h=0; h<GRID(x).length(); h++)
					sum += pow(GRID(x).value(h),2.0);
			}

			if ((mode!=InclField && mode!=Mobility) || included!=1)
				sum = std::sqrt(sum);
			assert(n<bufsize);
			buffer[n] = 255*((sum-min)/(max-min));
			if (levelset.size()>0) {// --contour
				double value = (lvlfield>-1) ? GRID(x)[lvlfield] : sum;
				for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
					if (std::fabs(value - *itl)/std::fabs(*itl)<contol)
						buffer[n] = coninv ? 255 : 0;
			}
			n++;
		}
	} else if (dim==2) {
		unsigned int n=0;
		MMSP::vector<int> x(2,0);
		for (x[1] = MMSP::g1(GRID,1)-1; x[1] >= MMSP::g0(GRID,1); x[1]--) {
			for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
				double sum=0.0;
				if (mode==InclField) { //  --field
					for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
						if (included>1)
							sum += pow(GRID(x)[*it],2.0);
						else
							sum = GRID(x)[*it];
					}
				} else if (mode==ExclField) { //  --exclude
					for (int h=0; h<GRID(x).length(); h++) {
						int i=GRID(x).index(i);
						std::set<int>::iterator it=fieldset.find(i);
						if (it == fieldset.end())
							sum += pow(GRID(x).value(h),2.0);
					}
				} else if (mode==Mobility) {
					sum = mobility(mu_lo, mu_hi, mu_x, mu_s, GRID(x).getMagPhi());
				} else { //          --mag or default
					for (int h=0; h<GRID(x).length(); h++)
						sum += pow(GRID(x).value(h),2.0);
				}
				if ((mode!=InclField && mode!=Mobility) || included!=1)
					sum = std::sqrt(sum);
				assert(n<bufsize);
				buffer[n] = 255*((sum-min)/(max-min));
				if (levelset.size()>0) {// --contour
					double value = (lvlfield>-1) ? GRID(x)[lvlfield] : sum;
					for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
						if (std::fabs(value - *itl)/std::fabs(*itl)<contol)
							buffer[n] = coninv ? 255 : 0;
				}
				n++;
			}
		}

		if (mode==CurveFit) {
			// Compute grain boundary profile
			x = tj;
			printf("Vertex sits at (%i, %i)\n", tj[0], tj[1]);

			int nexp = std::min(100, MMSP::x1(GRID)-x[0]); // number of points
			double* xexp = new double[nexp]; // raw data
			double* yexp = new double[nexp]; // raw data

			for (int i=0; i<nexp; i++) {
				MMSP::vector<int> s(x);
				double bottom = 2.0;
				while (s[1]<MMSP::y1(GRID)-10) {
					if (GRID(s).getMagPhi() < bottom) {
						bottom = GRID(s).getMagPhi();
						xexp[i] = MMSP::dx(GRID)*(s[0] - tj[0]);
						yexp[i] = MMSP::dy(GRID)*(s[1] - tj[1]);
					}
					s[1]++;
				}
				x[0]++;
			}

			curvefitter tricrystal(nexp, xexp, yexp);

			double a, theta;
			tricrystal.fit(a, theta);

			for (int i=0; i<nexp; i++) {
				int xi = tj[0] + xexp[i]/MMSP::dx(GRID);
				int yi = tj[1] + trijunction_profile(a, theta, xexp[i])/MMSP::dy(GRID);
				int width = MMSP::x1(GRID)-MMSP::x0(GRID);
				int height = MMSP::y1(GRID)-MMSP::y0(GRID);
				int index = (height-yi)*width + xi;
				buffer[index] = 0;
			}

			delete [] xexp;
			delete [] yexp;
		}
	} else if (dim==3) {
		unsigned int n=0;
		MMSP::vector<int> x(3,0);
		for (x[2] = MMSP::g0(GRID,2); x[2] < MMSP::g1(GRID,2); x[2]++)
			for (x[1] = MMSP::g1(GRID,1)-1; x[1] >= MMSP::g0(GRID,1); x[1]--)
				for (x[0] = MMSP::g0(GRID,0); x[0] < MMSP::g1(GRID,0); x[0]++) {
					if (x[sliceaxis]!=slicelevel) // clumsy, but effective
						continue;
					double sum=0.0;
					if (mode==InclField) { //  --field
						for (std::set<int>::iterator it=fieldset.begin(); it!=fieldset.end(); it++) {
							if (included>1)
								sum += pow(GRID(x)[*it],2.0);
							else
								sum = GRID(x)[*it];
						}
					} else if (mode==ExclField) { //  --exclude
						for (int h=0; h<GRID(x).length(); h++) {
							int i=GRID(x).index(i);
							std::set<int>::iterator it=fieldset.find(i);
							if (it == fieldset.end())
								sum += pow(GRID(x).value(h),2.0);
						}
					} else if (mode==CurveFit) { //  --fit
						for (int h=0; h<GRID(x).length(); h++)
							sum += pow(GRID(x).value(h) - 0.5,2.0);
					} else if (mode==Mobility) {
						sum = mobility(mu_lo, mu_hi, mu_x, mu_s, GRID(x).getMagPhi());
					} else { //          --mag
						for (int h=0; h<GRID(x).length(); h++)
							sum += pow(GRID(x).value(h),2.0);
					}

					if ((mode!=InclField && mode!=Mobility) || included!=1)
						sum = std::sqrt(sum);
					assert(n<bufsize);
					buffer[n] = 255*((sum-min)/(max-min));
					if (levelset.size()>0) {// --contour
						double value = (lvlfield>-1) ? GRID(x)[lvlfield] : sum;
						for (std::set<double>::iterator itl=levelset.begin(); itl!=levelset.end(); itl++)
							if (std::fabs(value - *itl)/std::fabs(*itl)<contol)
								buffer[n] = coninv ? 255 : 0;
					}
					n++;
				}
	}
}

// Implementation of curvefitter class member functions
int trijprof_f(const gsl_vector * k, void *data, gsl_vector* f)
{
	int n = ((struct data *)data)->n;
	double *x = ((struct data *)data)->x;
	double *y = ((struct data *)data)->y;

	double a = gsl_vector_get(k, 0);
	double theta = gsl_vector_get(k, 1);

	for (int i=0; i<n; i++) {
		// y = -a*(-theta + pi/2)/(2*theta) + a*acos(exp(-2*theta*x/a)*sin(theta))/(2*theta)
		double Fi = trijunction_profile(a, theta, x[i]);
		gsl_vector_set(f, i, Fi - y[i]);
	}
	return GSL_SUCCESS;
}

int trijprof_df(const gsl_vector * k, void *data, gsl_matrix* J)
{
	int n = ((struct data *)data)->n;
	double *x = ((struct data *)data)->x;

	double a = gsl_vector_get(k, 0);
	double theta = gsl_vector_get(k, 1);

	for (int i=0; i<n; i++) {
		// dy/da = -(-theta + pi/2)/(2*theta) + acos(exp(-2*theta*x/a)*sin(theta))/(2*theta) - x*exp(-2*theta*x/a)*sin(theta)/(a*sqrt(1 - exp(-4*theta*x/a)*sin(theta)**2))
		double dFda = -(-theta + M_PI/2.0)/(2.0*theta)
		            + acos(exp(-2.0*theta*x[i]/a)*sin(theta))/(2.0*theta)
		            - x[i]*exp(-2.0*theta*x[i]/a)*sin(theta)/(a*sqrt(1.0 - exp(-4.0*theta*x[i]/a)*pow(sin(theta),2)));

		// dy/dθ = a/(2*theta) - a*(exp(-2*theta*x/a)*cos(theta) - 2*x*exp(-2*theta*x/a)*sin(theta)/a)/(2*theta*sqrt(1 - exp(-4*theta*x/a)*sin(theta)**2)) + a*(-theta + pi/2)/(2*theta**2) - a*acos(exp(-2*theta*x/a)*sin(theta))/(2*theta**2)
		double dFdt = a/(2.0*theta) - a*(exp(-2.0*theta*x[i]/a)*cos(theta) - 2.0*x[i]*exp(-2.0*theta*x[i]/a)*sin(theta)/a)
		                              /(2.0*theta*sqrt(1.0 - exp(-4.0*theta*x[i]/a)*pow(sin(theta),2)))
		            + a*(-theta + M_PI/2.0)/(2.0*theta*theta)
		            - a*acos(exp(-2.0*theta*x[i]/a)*sin(theta))/(2.0*theta*theta);

		gsl_matrix_set (J, i, 0, dFda);
		gsl_matrix_set (J, i, 1, dFdt);
	}
	return GSL_SUCCESS;
}

curvefitter::curvefitter(int N, double* X, double* Y) :
	n(N), p(2), maxiter(20), T(gsl_multifit_fdfsolver_lmsder)
{
		d.n = N;
		d.x = X;
		d.y = Y;

		J = gsl_matrix_alloc(n, p);
		covar = gsl_matrix_alloc(p, p);
		k = gsl_vector_alloc(p);
		gsl_vector_set(k, 0, 64.0); // grain size guess
		gsl_vector_set(k, 1, M_PI/3.0); // dihedral angle guess

		f.f = &trijprof_f;
		f.df = &trijprof_df;
		f.n = n;
		f.p = p;
		f.params = &d;

		s = gsl_multifit_fdfsolver_alloc(T, n, p);
}

void curvefitter::fit(double& a, double& theta)
{
	const double xtol = 1.0e-8;
	const double gtol = 1.0e-8;
	const double ftol = 1.0e-8;

	// initialize solver with starting point and weights
	gsl_multifit_fdfsolver_set (s, &f, k);

	// compute initial residual norm
	res_f = gsl_multifit_fdfsolver_residual(s);
	chi0 = gsl_blas_dnrm2(res_f);

	// solve
	status = gsl_multifit_fdfsolver_driver(s, maxiter, xtol, gtol, ftol, &info);
	gsl_multifit_fdfsolver_jac(s, J);
	gsl_multifit_covar (J, 0.0, covar);

	// compute final residual norm
	chi = gsl_blas_dnrm2(res_f);

	#define FIT(i) gsl_vector_get(s->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

	fprintf(stdout, "summary from %s fit: %s after %zu iterations (%s)\n",
	        gsl_multifit_fdfsolver_name(s), gsl_strerror(status), gsl_multifit_fdfsolver_niter(s),
	        (info == 1) ? "threshold step size" : "threshold gradient");
	{
		double dof = n - p;
		double c = GSL_MAX_DBL(1, chi / sqrt(dof));

		fprintf(stdout, "chisq/dof = %g\n",  pow(chi, 2.0) / dof);

		fprintf (stdout, "  a   = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
		fprintf (stdout, "theta = %.5f +/- %.5f\n", 180.0/M_PI*FIT(1), 180.0/M_PI*c*ERR(1));
	}

	a = FIT(0);
	theta = FIT(1);
}

curvefitter::~curvefitter()
{
		gsl_multifit_fdfsolver_free(s);
		gsl_matrix_free(covar);
		gsl_matrix_free(J);
		gsl_vector_free(k);
}
