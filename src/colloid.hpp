/*!
 * author: Nathan A. Mahynski
 *
 * Functionalized colloid.
 *
 * Described as a combination of a tile and motif.
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

using namespace csk;
using namespace std;
using namespace glm;

class Colloid
{
public:
        Colloid();
        ~Colloid();

        // Unrolled parameter array for optimization, etc.
        void setParameters( const vector<double> &params );
        const vector<double> getParameters();

        // Assign the motif
        void setMotif( Motif m );
	const Motif getMotif() { return m_; }

        // Assign the tile
        void setTile( IsohedralTiling t );
	const IsohedralTiling getTile() { return tile_; }

        bool isMotifInside( const int N ); // Check if the motif is inside
        double tileArea(); // Compute the area of the tile

	void build() { buildBoundary_(); buildMotif_(); }

private:
	void buildBoundary_();
	void buildMotif_();
	vector<vector<dvec2>> perimeter_( double u0, double du, int n );

	bool tile_assigned_, motif_assigned_;
        IsohedralTiling tile_;
        Motif m_;

        vector<double> params_; // Unrolled parameter vector

	vector<int> boundary_ids_; // Chemical identities of boundary points
	vector<vector<double>> boundary_coords_; // Points on tile's boundary
	double sphere_deform_, edge_du_, edge_u0_;
	vector<vector<double>> tile_control_points_;
};

#endif
