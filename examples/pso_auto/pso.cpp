/**
 * Particle Swarm Optimization (PSO) is a stochastic swarm intelligence algorithm for global optimization of potentially ill-behaved nonlinear functions.
 */

#define OPTIM_ENABLE_ARMA_WRAPPERS
#include "optim.hpp"

#include <iostream>
#include <vector>
#include <cmath>

#include "tiling.hpp"
#include "src/motif.hpp"
#include "src/colloid.hpp"

struct fn_data {
	Colloid* c;
	int verbosity;
	double A_target;
	double penalty;
	double df_min;
};

double area_error2(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
	/**
	 * Compute the square deviation of the area of a tile from a target value.
	 * 
	 * @param vals_inp Parameter vector, see `Colloid`.
	 * @param grad_out Empty vector, PSO does not use gradients.
	 * @param opt_data Array of [Colloid].
	 */

	fn_data* data = reinterpret_cast<fn_data*>(opt_data);

	vector<double> p = data->c->getParameters();
	for( size_t i = 0; i < p.size(); ++i ) {
		p[i] = vals_inp(i);
	}
	try {
		data->c->setParameters(p, data->df_min);
	} catch ( customException& e) {
		// Report the error and return a large value to allow the code to try something
		// else without crashing. 
                // Exceptions can be thrown when a tile is self-intersecting, for example, so 
                // we highly penalize this configuration and move on.
                if (data->verbosity > 0) {
			std::cerr << e.getMessage() << std::endl;
		}
		return data->penalty;
	}

	double raw_tile_area = data->c->tileArea();
	double f_in = data->c->fractionMotifInside(20);
	double f_out = 1.0 - f_in;

	/**
	 * The logic behind this functional form is as follows.  We want to target a given tile area
	 * so (area - target)**2 is the base function to minimize.  However, we need to penalize the 
	 * sytem if the motif cannot be fully enclosed.  The more the motif is outside, the bigger
	 * the penalty should be so a term that looks like ~f_out*penalty will monotonically increase
	 * as the tile excludes more of the motif.  We should also consider that small changes to the
	 * tile shape may not change the discrete number of points inside or outside the tile, yet
	 * we want to bias the system to move toward expanding the tile if some parts of the motif
	 * do not fit.  This is not always perfect, e.g., if the tile is plenty big but the motif is
	 * awkwardly placed in a corner.  Nonetheless, an additional bias ~ 1/area should
	 * help the system behave as we would like. This is mostly important when using gradient-
	 * based optimization approaches, which PSO doesn't use. Still, it is nice to have this so
	 * performance can be compared to optimization routines which do use derivatives. Note that
	 * f_out is outside the parentheses because if it is 0, no additional penalty should be 
	 * applied. Note that tile area is usually on the order of unity for motifs in the included
	 * library, but penalty should be chosen to be commensurate with it.
	 */
	double loss_function = pow(f_out*(data->penalty + 1.0/raw_tile_area) + raw_tile_area - data->A_target, 2);

	return loss_function;
}

static void show_usage(char *name)
{
    std::cerr << "Usage: " << name << " <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-n,--number\t\tISOHEDRAL TILE\tSpecify the IH Tile number\n"
	      << "\t-m,--motif\t\tMOTIF\tJSON file with motif in it\n"
	      << "\t-px,--px\tPARAMETERX\tSpecify tile parameter x"
              << std::endl;
}

int main(int argc, char **argv)
{
	// Command line arguments
	if ( argc < 3 ) {
		show_usage(argv[0]);
		return -1;
	}
	
	int ih_number = -1, idx = -1;
	map<int, float> p;
	string pname = "-p", motif = "NONE";

	for (int i=0; i < argc ; ++i) {
		string arg = argv[i];
		if ( (arg == "-h") || (arg == "--help") ) {
			show_usage(argv[0]);
			return 0;
		} else if ( (arg == "-n") || (arg == "--number") ) {
			if ( i+1 < argc ) {
				ih_number = atoi(argv[i+1]);
				for ( size_t j = 0; j < sizeof(csk::tiling_types)/sizeof(csk::tiling_types[0]); ++j) {
					if ( int(csk::tiling_types[j]) == ih_number ) {
						idx = j;
						break;
					}
				}
				if ( idx < 0 ) {
					std::cerr << ih_number << " is an invalid IH tile number" << endl;
					return -1;
				}
			} else {
				std::cerr << arg << " option requires one argument" << endl;
				return -1;
			}
		} else if ( (arg == "-m") || (arg == "--motif") ) {
			if ( i+1 < argc ) {
				motif = string(argv[i+1]);
			} else {
				std::cerr << arg << " option requires one argument" << endl;
				return -1;
			}
		} else if ( "-p" == arg.substr(0, 2) ) {
			int n = atoi( arg.substr(2, arg.length()).c_str() );
			if ( (n < 0) || (n > 5) ) { // At most there are 6 tile parameters
				std::cerr << "invalid parameter index " << arg << endl;
				return -1;
			} else {
				if ( i+1 < argc ) {
					p[n] = atof(argv[i+1]);
				} else {
					std::cerr << arg << " option requires one argument" << endl;
	                                return -1;
				}
			}
		} 
	}

	// Create the tile
	IsohedralTiling t(ih_number); // Use default tile parameters - usually good

	// Create an array to hold a copy of the tiling vertex parameters.
        double params[t.numParameters()], default_params[t.numParameters()];
	t.getParameters(default_params);

	// Iterate over parameters to make sure we have them all, given the tile choice
	for (unsigned int i=0; i < int(t.numParameters()); ++i) {
		if (p.find(i) == p.end()) {
			std::cerr << "p" << i << " not found; using a default value" << endl;
			params[i] = default_params[i];
		} else {
			params[i] = p[i];
		}
	}	
	t.setParameters(params);

	// Data
	fn_data data;
	data.verbosity = 0; // Set to > 0 to print error messages / information
	data.penalty = 1000.0;
	data.A_target = 0.0; // Minimize the area ("escherization" problem)
	data.df_min = 0.1; // Enforce a minimum curvature

	// Create the colloid
	arma::vec x_1;
	data.c = new Colloid();
	try {
		Motif m;
		m.load(motif); // Select a motif
		data.c->setMotif(m);
		data.c->setTile(t);
		vector<double> u0(t.numEdgeShapes(), 0.1);
                vector<double> df(t.numEdgeShapes(), 0.1); 
		data.c->setU0(u0);
		data.c->setDform(df);
		data.c->setDU(0.1);
		data.c->setTileScale(1.0);
		data.c->init();
		x_1 = data.c->getParameters(); // Initial guess is result after initialization
	} catch (const exception& e) {
		std::cerr << "unable to initialize colloid" << std::endl;
		return -1;
	}

	arma::vec lb = arma::zeros(x_1.size(), 1); // Lower bounds
	lb(0) = 0.0; // Motif scaled_com_x
	lb(1) = 0.0; // Motif scaled_com_y
	lb(2) = 0.0; // Motif theta
	for (unsigned int i=0; i < t.numParameters(); ++i) {
		lb(3+i) = 0.0; // See https://isohedral.ca/software/tactile/ to estimate visually
	}

	// du/2 <= u0 <= 1-5*du-du/2.
	for (unsigned int i=0; i < t.numEdgeShapes(); ++i) {
		lb(3+t.numParameters()+i) = 0.1; // edge u0
	}

	for (unsigned int i=0; i < t.numEdgeShapes(); ++i) {
                lb(3+t.numParameters()+t.numEdgeShapes()+i) = -0.5; // edge df
        }

	// Tile scale - use a factor on the initial scale found
	lb(3+t.numParameters()+2*t.numEdgeShapes()) = 0.25*x_1[x_1.size()-1]; 

	arma::vec ub = arma::zeros(x_1.size(), 1); // Upper bounds
	ub(0) = 1.0; // Motif scaled_com_x
	ub(1) = 1.0; // Motif scaled_com_y
	ub(2) = 2.0*arma::datum::pi; // Motif theta

	for (unsigned int i=0; i < t.numParameters(); ++i) {
                ub(3+i) = 2.0; // See https://isohedral.ca/software/tactile/ to estimate visually
        }

	for (unsigned int i=0; i < t.numEdgeShapes(); ++i) {
                ub(3+t.numParameters()+i) = 0.4; // 1.0 - 5*0.1 - 0.1/2 = 0.45
        }

	for (unsigned int i=0; i < t.numEdgeShapes(); ++i) {
                ub(3+t.numParameters()+t.numEdgeShapes()+i) = 0.5; // edge df
        }

	// Tile scale - use a factor on the initial scale found
	ub(3+t.numParameters()+2*t.numEdgeShapes()) = 2.0*x_1[x_1.size()-1]; 

	optim::algo_settings_t settings_1;
	
	settings_1.vals_bound = true; 
	settings_1.lower_bounds = lb;
	settings_1.upper_bounds = ub;
 	settings_1.print_level = 2;

	settings_1.pso_settings.center_particle = true; //false;
	settings_1.pso_settings.n_pop = 500; // population size of each generation.
	settings_1.pso_settings.n_gen = 30; // number of vectors to generate (iterations).
 
	bool success_1 = optim::pso(x_1, area_error2, &data, settings_1);
 
	if (success_1) {
		std::cout << "pso: Area minimization completed successfully." << std::endl;
		std::vector<double> f(x_1.size(), 0.0);
		for (unsigned int i=0; i < f.size(); ++i) {
			f[i] = x_1[i];
		}
		try {
			// "Revisions" to respect symmetry might change the parameters
			// so set() and get() once to return the "true" values being
			// used.
			data.c->setParameters(f);
			data.c->dumpXYZ("success.xyz", true);
			data.c->dump("success.json");
		} catch (const customException &e) {
			std::cerr << e.getMessage() << std::endl;
			return -1;
		}
	} else {
		std::cout << "pso: Area minimization completed unsuccessfully." << std::endl;
	}

        x_1 = data.c->getParameters();	
	arma::cout << "pso: solution to Area minimization test:\n" << x_1 << arma::endl;
	
	delete data.c;
}
