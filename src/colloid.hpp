/**
 * A colloid is composed of an isohedral tile boundary and internal motif.
 * 
 * This file contains functions to build and manipulate colloids. Tile 
 * information is inherited from the Tactile library, while a Motif
 * class is described in motif.hpp. This should be self-contained so 
 * that, for example, external optimizers can manipulate the colloid's
 * parameters controlling its shape, etc. and compute properties from
 * that.
 * 
 * @author Nathan A. Mahynski
 */

#ifndef COLLOID_H_
#define COLLOID_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <bits/stdc++.h>

#include "tiling.hpp"
#include "motif.hpp"
#include "utils.hpp"
#include "json.hpp"

using namespace csk;
using namespace std;
using namespace glm;
using json = nlohmann::json;

class Colloid
{
public:
        Colloid();
	Colloid( Motif m, IsohedralTiling t , double tile_scale, double tile_u0 );
        ~Colloid();

	vector<double> unscale_coords_( const vector<double> &scaled_coords );
	vector<double> scale_coords_( const vector<double> &unscaled_coords );
        void setParameters( const vector<double> &params );
        const vector<double> getParameters();

        void setMotif( Motif m );
	const Motif getMotif();

        void setTile( IsohedralTiling t );
	const IsohedralTiling getTile();
	void setTileScale( const double s ) { tile_scale_ = s; } // Assign the tile_scale_.
	double tileArea();

        bool isMotifInside( const int N );
	double fractionMotifInside( const int N );

	void init() { buildBoundary_(); initMotif_( 5.0, 0.2, 1000, 20 ); built_ = true; } // Initialize the colloid.

	void load( const string filename );
	void dump( const string filename );
	void dumpXYZ( const string filename );	

	vector<double> boundaryCOM();

private:
	void buildBoundary_();
	void initMotif_( double max_scale_factor, double min_scale_factor, int n_scale_incr, int N );
	void perimeter_( double u0, double du, int n, double scale, 
		vector<int>* boundary_ids, 
		vector<vector<double>>* boundary_coords, 
		vector<vector<double>>* tile_control_points );
	vector<vector<dvec2>> perimeter_edges_( double u0, double du, int n, double scale );

	bool tile_assigned_; // Has the tile been assigned yet?
	bool motif_assigned_; // Has the motif been assigned yet?
	bool built_; // Has the colloid been constructed at least once?

	double sphere_deform_; // Normalized amount a sphere "deforms" the edge.
	double edge_du_; // Parameterized (Bezier) gap between boundary points.
	double edge_u0_; // Parameterized (Bezier) starting point for boundary points.
	double tile_scale_; // The default Tactile tile is isotropically scaled by this factor.

	vector<int> boundary_ids_; // Chemical identities of boundary points.
	vector<double> params_; // Unrolled parameter vector.

	vector<vector<double>> boundary_coords_; // Coordinates of points on tile's boundary.
	vector<vector<double>> tile_control_points_; // Control points on Bezier curves defining edges.

        IsohedralTiling tile_; // Isohedral tile from Tactile library.
        Motif m_; // The colloid's motif.	
};

#endif
