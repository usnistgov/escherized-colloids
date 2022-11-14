/**
 * Particle Swarm Optimization (PSO) is a stochastic swarm intelligence algorithm for global optimization of potentially ill-behaved nonlinear functions.
 */

#define OPTIM_ENABLE_ARMA_WRAPPERS
#include "optim.hpp"

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "tiling.hpp"
#include "src/motif.hpp"
#include "src/colloid.hpp"

struct fn_data {
	Colloid* c;
	int verbosity;
	double A_target;
	double penalty;
	double df_min;
	double eps;
	double n;
	double cutoff_ratio;
	double skin;
	double side_penalty;
	double side_ratio_limit;
};

static void show_usage(char *name)
{
    std::cerr << "Usage: " << name << " <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
	      << "\t-l,--label\t\tLABEL\tPrefix to label output files with\n"
              << "\t-n,--number\t\tISOHEDRAL TILE\tSpecify the IH Tile number\n"
	      << "\t-m,--motif\t\tMOTIF\tJSON file with motif in it\n"
	      << "\t-u,--du\t\tDU\tGap between boundary points\n"
	      << "\t-e,--eps\t\tEPSILON\tThe interaction energy scale between motif and boundary points\n"
	      << "\t-px,--px\tPARAMETERX\tSpecify tile parameter x\n"
	      << "\t-a,--min_angle\tMINIMUM_ANGLE\tMinimum angle allowable at tile corners in DEGREES\n"
              << "\t-g,--min_gap\tMINIMUM_GAP\tMinimum distance allowable between centers of tile edges."
	      << std::endl;
}

const double pairwise(const double r, const double eps, const double n, const double cutoff_ratio, const double skin) 
{
	/**
	 * Compute the energy of interaction between two points.
	 * We shift this 1/r^n potential so it is always finite and equal to -eps at r=skin.
	 * The cutoff is determined so that U(r_cut)/(-eps) = 1/cutoff_ratio.
	 */
	assert(eps >= 0.0);
	assert(n >= 1.0);
	assert(cutoff_ratio >= 1.0);
	assert(skin >= 0.0);	
	const double rc = pow(cutoff_ratio, 1.0/n) - skin + 1.0;
	if ((r >= rc) || (r < skin)) {
		return 0.0;
	} else {
		return -eps/pow(r + 1.0 - skin, n);
	}
}

const double energy(Colloid& c, const double eps, const double n, const double cutoff_ratio, const double skin) 
{
	/**
	 * Compute the pairwise energy between the boundary points and the motif.
	 */
	const vector<vector<double>> mc = c.getMotif().getCoords();
	const vector<vector<double>> bc = c.getBoundaryCoords();

	double u = 0.0;
	for (unsigned int i = 0; i < mc.size(); ++i) {
		for (unsigned int j = 0; j < bc.size(); ++j) {
			const double r = sqrt(distance2(mc[i], bc[j])); 
			u = std::min(u, pairwise(r, eps, n, cutoff_ratio, skin));
		}
	}

	return u;
}

const double side_length_penalty(Colloid& c, const double penalty, const double max_min_limit = 10.0) {
	/**
	 * Penalize the configuration by the ratio between the longest and shortest distance
	 * between control points on the tile.
	 */
	assert(penalty >= 0.0);
	assert(max_min_limit >= 1.0);

	const vector<vector<double>> cp = c.getTileControlPoints();
	double min_d2 = PIP_INF, max_d2 = 0.0;
	for (size_t i=0; i < cp.size(); i++) {
		for (size_t j=i+1; j < cp.size(); ++j) {
			const double d2 = distance2(cp[i], cp[j]);
			if (d2 < min_d2) {
				min_d2 = d2;
			}
			if (d2 > max_d2) {
				max_d2 = d2;
			}
		}
	}
	
	const double delta = sqrt(max_d2/min_d2)/max_min_limit - 1.0;
	if (delta > 0.0) {
		return penalty*delta;
	}
	return 0.0;
}

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

	// First check if the tile shape is even "legal". (no intersection, bad angles, all edges > min_gap, no "overrun", etc.)
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

	// If the tile shape is "legal" then we need to optimize the motif inside of it.
	double raw_tile_area = data->c->tileArea();
	double f_in = data->c->fractionMotifInside(20, data->skin);
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
	 *
	 * We also add an energetic "benefit" to try to move the boundary points to be "close" to
	 * the motif's points for practical purposes. The 1/r^n potential used means that energy is
	 * strictly minimized for an individual bead when it coincides with the tile edge (+skin); 
	 * for a many-body system the minimum could be to place the motif so it is centered on the 
	 * edge with some points inside and others outside of the tile, which is undesirable. We add
	 * the (negative) energy only when the motif is entirely "inside" the tile. 
	 *
	 * The side length penalty can be used to try to keep the edges of the tile from going to
	 * zero. This is applied only when the motif is fully "inside" the tile. 
	 */
	const double motif_fit_penalty = f_out*(data->penalty + 1.0/raw_tile_area);
	const double effective_area = motif_fit_penalty + raw_tile_area;
	double loss_function = pow(effective_area - data->A_target, 2);
        if (f_in == 1.0) {
		loss_function += energy(*data->c, data->eps, data->n, data->cutoff_ratio, data->skin);
		loss_function += side_length_penalty(*data->c, data->side_penalty, data->side_ratio_limit);
	}

	return loss_function;
}

int main(int argc, char **argv)
{
	// Command line arguments
	if ( argc < 3 ) {
		show_usage(argv[0]);
		return 1;
	}
	
	int ih_number = -1, idx = -1;
	map<int, float> p;
	string pname = "-p", motif = "NONE", label = "";
	double du = 0.1, eps = 0.0, min_angle = 0.0, min_gap = 0.0;

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
					return 1;
				}
			} else {
				std::cerr << arg << " option requires one argument" << endl;
				return 1;
			}
		} else if ( (arg == "-m") || (arg == "--motif") ) {
			if ( i+1 < argc ) {
				motif = string(argv[i+1]);
			} else {
				std::cerr << arg << " option requires one argument" << endl;
				return 1;
			}
		} else if ( (arg == "-l") || (arg == "--label") ) {
			if ( i+1 < argc ) {
				label = string(argv[i+1]);
			} else {
				std::cerr << arg << " option requires one argument" << endl;
				return 1;
			}
		} else if ( (arg == "-u") || (arg == "--du") ) {
			if ( i+1 < argc ) {
                                du = atof(argv[i+1]);
                        } else {
                                std::cerr << arg << " option requires one argument" << endl;
                                return 1;
                        }
		} else if ( (arg == "-a") || (arg == "--min_angle") ) {
                        if ( i+1 < argc ) {
                                min_angle = atof(argv[i+1])*M_PI/180.0;
				assert(min_angle >= 0.0);
				assert(min_angle <= M_PI);
                        } else {
                                std::cerr << arg << " option requires one argument" << endl;
                                return 1;
                        }
		} else if ( (arg == "-g") || (arg == "--min_gap") ) {
                        if ( i+1 < argc ) {
                                min_gap = atof(argv[i+1]);
                                assert(min_gap >= 0.0);
                        } else {
                                std::cerr << arg << " option requires one argument" << endl;
                                return 1;
                        }
		} else if ( (arg == "-e") || (arg == "--eps") ) {
                        if ( i+1 < argc ) {
                                eps = atof(argv[i+1]);
                        } else {
                                std::cerr << arg << " option requires one argument" << endl;
                                return 1;
                        }
		} else if ( "-p" == arg.substr(0, 2) ) {
			int n = atoi( arg.substr(2, arg.length()).c_str() );
			if ( (n < 0) || (n > 5) ) { // At most there are 6 tile parameters
				std::cerr << "invalid parameter index " << arg << endl;
				return 1;
			} else {
				if ( i+1 < argc ) {
					p[n] = atof(argv[i+1]);
				} else {
					std::cerr << arg << " option requires one argument" << endl;
	                                return 1;
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
	data.penalty = 100.0;
	data.A_target = 0.0; // Minimize the area ("escherization" problem)
	data.df_min = 0.1; // Enforce a minimum curvature
	data.eps = eps; // Interaction energy - 0 implies no interaction between tile and boundary
	data.n = 24;
	data.cutoff_ratio = 100.0;
	data.side_penalty = 10.0*eps; 
	data.side_ratio_limit = 5.0;

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
		data.c->setDU(du);
		data.c->setTileScale(1.0);
		data.c->setMinAngle(min_angle);
		data.c->setMinGap(min_gap);
		data.c->init();
		data.skin = 1.13*m.minDistance()/2.0; // Fudge factor for skin
		x_1 = data.c->getParameters(); // Initial guess is result after initialization
	} catch (const customException& e) {
		std::cerr << "unable to initialize colloid : " << e.getMessage() << std::endl;
		return 1;
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
		lb(3+t.numParameters()+i) = du/2.0; // edge u0
	}

	for (unsigned int i=0; i < t.numEdgeShapes(); ++i) {
                lb(3+t.numParameters()+t.numEdgeShapes()+i) = -0.5; // edge df
        }

	// Tile scale - use a factor on the initial scale found
	lb(3+t.numParameters()+2*t.numEdgeShapes()) = 0.5*x_1[x_1.size()-1]; 

	arma::vec ub = arma::zeros(x_1.size(), 1); // Upper bounds
	ub(0) = 1.0; // Motif scaled_com_x
	ub(1) = 1.0; // Motif scaled_com_y
	ub(2) = 2.0*arma::datum::pi; // Motif theta

	for (unsigned int i=0; i < t.numParameters(); ++i) {
                ub(3+i) = 2.0; // See https://isohedral.ca/software/tactile/ to estimate visually
        }

	for (unsigned int i=0; i < t.numEdgeShapes(); ++i) {
                ub(3+t.numParameters()+i) = 1.0 - 5*du - du/2; // edge u0
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

	settings_1.pso_settings.center_particle = false;
	settings_1.pso_settings.n_pop = 10000; // population size of each generation.
	settings_1.pso_settings.n_gen = 100; // number of vectors to generate (iterations).
 
	bool success_1 = optim::pso(x_1, area_error2, &data, settings_1);
	bool success_2 = (data.c->fractionMotifInside(20, data.skin) > 0.99999999);

	if (success_1 && success_2) {
		std::cout << "pso: Area minimization completed successfully." << std::endl;
		std::vector<double> f(x_1.size(), 0.0);
		for (unsigned int i=0; i < f.size(); ++i) {
			f[i] = x_1[i];
		}
		try {
			// "Revisions" to respect symmetry might change the parameters
			// so set() and get() once to return the "true" values being
			// used.
			data.c->setParameters(f); // Relative position of motif unchanged
			stringstream ss;
			ss << label << "_success.xyz";
			data.c->dumpXYZ(ss.str(), true);
			ss.str("");
			ss << label << "_success.json";
			data.c->dump(ss.str());

			vector<vector<double>> c_, b_;
                	vector<string> t_;
	                data.c->unitCell(&c_, &t_, &b_, 5, 5, 1.0e-8, false, true);	
			ss.str("");
			ss << label << "_5x5_wb_unit_cell.xyz";
			dumpXYZ(c_, t_, ss.str());
			data.c->unitCell(&c_, &t_, &b_, 5, 5, 1.0e-8, false, false);
			ss.str("");
			ss << label << "_5x5_unit_cell.xyz";
        	        dumpXYZ(c_, t_, ss.str());
		} catch (const customException &e) {
			std::cerr << e.getMessage() << std::endl;
			return 1;
		}
	} else {
		std::cout << "pso: Area minimization completed unsuccessfully : ";
		if (!success_2) {
			std::cout << "motif was not fully enclosed." << std::endl;
		} else {
			std::cout << "pso failed to converge." << std::endl;
		}
	}

        x_1 = data.c->getParameters();	
	arma::cout << "pso: solution to Area minimization test:\n" << x_1 << arma::endl;
	
	delete data.c;

	return 0;
}
