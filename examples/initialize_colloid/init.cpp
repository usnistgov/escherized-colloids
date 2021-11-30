/**
 * Create a colloid by choosing an isohedral tile, then loading
 * a motif into it. These colloids have their edges decorated using
 * the "4+1" rule.
 *
 * @author Nathan A. Mahynski
 *
 * Note from Tactile:
 * "Note that the program might randomly generate tiles that
 * self-intersect. That's not a bug in the library, it's just a
 * bad choice of tiling vertex parameters and edge shapes."
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "tiling.hpp"
#include "motif.hpp"
#include "colloid.hpp"

using namespace std;
using namespace csk;

static void show_usage(char *name)
{
    std::cerr << "Usage: " << name << " <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-n,--number\t\tISOHEDRAL TILE\tSpecify the IH Tile number\n"
	      << "\t-f,--filename\t\tFILENAME\tSpecify the file to write to\n"
	      << "\t-m,--motif\t\tMOTIF\tJSON file with motif in it\n"
	      << "\t-px,--px\tPARAMETERX\tSpecify tile parameter x"
              << std::endl;
}

int main( int argc, char **argv )
{
	if ( argc < 3 ) {
		show_usage(argv[0]);
		return 1;
	}
	
	int ih_number = -1, idx = -1;
	map<int, float> p;
	string filename = "NONE", pname = "-p", motif = "NONE";

	for ( int i=0; i < argc ; ++i ) {
		string arg = argv[i];
		if ( (arg == "-h") || (arg == "--help") ) {
			show_usage(argv[0]);
			return 0;
		} else if ( (arg == "-n") || (arg == "--number") ) {
			if ( i+1 < argc ) {
				ih_number = atoi(argv[i+1]);
				for ( size_t j = 0; j < sizeof(tiling_types)/sizeof(tiling_types[0]); ++j) {
					if ( int(tiling_types[j]) == ih_number ) {
						idx = j;
						break;
					}
				}
				if ( idx < 0 ) {
					cerr << ih_number << " is an invalid IH tile number" << endl;
					return 1;
				}
			} else {
				cerr << arg << " option requires one argument" << endl;
				return 1;
			}
		} else if ( (arg == "-m") || (arg == "--motif") ) {
			if ( i+1 < argc ) {
				motif = string(argv[i+1]);
			} else {
				cerr << arg << " option requires one argument" << endl;
				return 1;
			}
		} else if ( (arg == "-f") || (arg == "--filename") ) {
			if ( i+1 < argc ) {
				filename = string(argv[i+1]);
			} else {
				cerr << arg << " option requires one argument" << endl;
                                return 1;
			}
		} else if ( "-p" == arg.substr(0, 2) ) {
			int n = atoi( arg.substr(2, arg.length()).c_str() );
			if ( (n < 0) || (n > 5) ) { // At most there are 6 tile parameters
				cerr << "invalid parameter index " << arg << endl;
				return 1;
			} else {
				if ( i+1 < argc ) {
					p[n] = atof(argv[i+1]);
				} else {
					cerr << arg << " option requires one argument" << endl;
	                                return 1;
				}
			}
		} 
	}

	
	// Construct a tiling of the given type.
        IsohedralTiling t( tiling_types[idx] );

	// Create an array to hold a copy of the tiling vertex parameters.
        double params[ t.numParameters() ], default_params[ t.numParameters() ];
	t.getParameters( default_params );

	// Iterate over parameters to make sure we have them all, given the tile choice
	for ( int i=0; i < int(t.numParameters()); ++i ) {
		if ( p.find(i) == p.end() ) {
			cerr << "p" << i << " not found; using a default value" << endl;
			params[i] = default_params[i];
		} else {
			params[i] = p[i];
		}
	}	

	// Load the motif
	Motif m;
	m.load(motif);

	// Build the colloid
	try {
		Colloid c(m, t, 0.3);
		c.dump("colloid.json");
		c.dumpXYZ("colloid.xyz");
	} catch ( const exception& e ) {
		cout << e.what();
		return 1;
	}

	return 0;
}

