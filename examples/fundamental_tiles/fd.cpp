/**
 * Create an example of all fundamental tilings.
 *
 * @author Nathan A. Mahynski
 *
 * Note from Tactile:
 * "Note that the program might randomly generate tiles that
 * self-intersect. That's not a bug in the library, it's just a
 * bad choice of tiling vertex parameters and edge shapes."
 */

#include <iostream>

#include "tiling.hpp"
#include "src/colloid.hpp"
#include "src/utils.hpp"

using namespace std;
using namespace csk;

int main( int argc, char **argv )
{
	for (unsigned int j = 0; j < sizeof(tiling_types) / sizeof(tiling_types[0]); ++j) {
		// Supported types in Tactile
		int ih_type = static_cast<int>(tiling_types[j]);
		for (unsigned int i = 0; i < sizeof(FD_TYPES) / sizeof(FD_TYPES[0]); ++i) {
			if (ih_type == FD_TYPES[i]) { // Check if that type is fundamental
				try {
					// Construct a default tiling of the given type.
					IsohedralTiling t( ih_type );

					// Load the motif
					Motif m;
					m.load("../../motif_library/c1_random.json");
					
					// Create tile
					vector<double> u0(t.numEdgeShapes(), 0.25);
					vector<double> df(t.numEdgeShapes(), 0.0); // No curvature
					Colloid c(m, t, u0, df, 0.1, false);
					stringstream ss, tt;
					ss << "colloid_" << int(ih_type) << ".xyz";
					c.dumpXYZ(ss.str(), true);
					tt << "colloid_" << int(ih_type) << ".json";
					c.dump(tt.str());
				} catch (const customException &e) {
					std::cout << "Error for tile type " << int(ih_type) << " : " << e.getMessage() << std::endl;
					

					return -1;
				}
			}
		}
	}

	return 0;
}

