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
	double A_target;
	double penalty;
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
	data->c->setParameters(p);

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
	double objective_function = pow(f_out*(data->penalty + 1.0/raw_tile_area) + raw_tile_area - data->A_target, 2);

	return objective_function;
}

int main()
{
	// Data
	fn_data data;
	data.penalty = 1000.0;
	data.A_target = 0.0; // Minimize the area ("escherization" problem)

	// Create a colloid
	arma::vec x_1;
	data.c = new Colloid();
	try {
		Motif m;
		m.load("../../motif_library/d1_vitruvian.json"); // Select a motif
		IsohedralTiling t( 7 ); // Use default tile parameters - usually good
		data.c->setMotif(m);
		data.c->setTile(t);
		data.c->init();
		x_1 = data.c->getParameters(); // Initial guess is result after initialization
	} catch ( const exception& e ) {
		std::cerr << "unable to initialize colloid" << std::endl;
		return -1;
	}

	arma::vec lb = arma::zeros(x_1.size(), 1); // Lower bounds
	lb(0) = 0.0; // Motif scaled_com_x
	lb(1) = 0.0; // Motif scaled_com_y
	lb(2) = 0.0; // Motif theta
	lb(3) = 0.0; // IH07 v_0 - see https://isohedral.ca/software/tactile/ to estimate visually
	lb(4) = 0.0; // IH07 v_1 - see https://isohedral.ca/software/tactile/ to estimate visually
	lb(5) = 0.1; // edge u0
	lb(6) = 0.5; // Tile scale

	arma::vec ub = arma::zeros(x_1.size(), 1); // Lower bounds
	ub(0) = 1.0; // Motif scaled_com_x
	ub(1) = 1.0; // Motif scaled_com_y
	ub(2) = 2.0*arma::datum::pi; // Motif theta
	ub(3) = 2.0; // IH07 v_0 - see https://isohedral.ca/software/tactile/ to estimate visually
	ub(4) = 2.0; // IH07 v_1 - see https://isohedral.ca/software/tactile/ to estimate visually
	ub(5) = 0.5; // edge u0
	ub(6) = 3.0; // Tile scale

	optim::algo_settings_t settings_1;
	
	settings_1.vals_bound = true; 
	settings_1.lower_bounds = lb;
	settings_1.upper_bounds = ub;
 	settings_1.print_level = 0;

	settings_1.pso_settings.center_particle = false;
	settings_1.pso_settings.n_pop = 50; // population size of each generation.
	settings_1.pso_settings.n_gen = 40; // number of vectors to generate (iterations).
 
	bool success_1 = optim::pso(x_1, area_error2, &data, settings_1);
 
	if (success_1) {
		std::cout << "pso: Area minimization completed successfully." << std::endl;
		data.c->dumpXYZ("success.xyz", true);
		data.c->dump("success.json");
	} else {
		std::cout << "pso: Area minimization completed unsuccessfully." << std::endl;
	}
 
	arma::cout << "pso: solution to Area minimization test:\n" << x_1 << arma::endl;

	delete data.c;
}
