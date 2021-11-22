/*!
 * author: Nathan A. Mahynski
 */

#include "colloid.hpp"

bool pip(const vector<dvec2> &polygon, const vector<double>& point) {
	// Based on https://stackoverflow.com/questions/11716268/point-in-polygon-algorithm
	bool c = false;
	for (size_t i=0, j=polygon.size()-1; i < polygon.size(); j=i++) {
		if( 
			((polygon[i].y >= point[1]) != (polygon[j].y >= point[1])) &&
		  	(point[0] <= (polygon[j].x - polygon[i].x)*(point[1] - polygon[i].y)/(polygon[j].y-polygon[i].y) + polygon[i].x)
		) {
			c = !c;
		}
	}
	return c;
}

dvec2 bezier(dvec2 p0, dvec2 p1, dvec2 p2, dvec2 p3, double u) {
	/*!
	 * Compute parameterized location along Bezier curve.
	 */
	
	// From https://www.geeksforgeeks.org/cubic-bezier-curve-implementation-in-c/
	dvec2 coord;	
	coord.x = pow(1.0-u, 3)*p0.x + 3.0*u*pow(1.0-u, 2)*p1.x + 3.0*(1.0-u)*pow(u, 2)*p2.x + pow(u, 3)*p3.x;
	coord.y = pow(1.0-u, 3)*p0.y + 3.0*u*pow(1.0-u, 2)*p1.y + 3.0*(1.0-u)*pow(u, 2)*p2.y + pow(u, 3)*p3.y;

	return coord;
}

Colloid::Colloid(): tile_(IsohedralTiling(0)) {
	tile_assigned_ = false;
	motif_assigned_ = false;
	sphere_deform_ = 0.25; // The amount a sphere "deforms" the edge
	edge_du_ = 0.1; // Parameterized (Bezier) gap between edge points
}

Colloid::~Colloid() {
}

void Colloid::setParameters( const vector<double> &params ) {
	// This will have to be more complex - when updating, make sure the tile and motif are
	// consistent with each other (same COM, for example).


	if (!motif_assigned_ || !tile_assigned_) {
		// Make sure the tile and motif have been assigned first.
		throw customException("must assign colloid and tile before parameters");
	}

	// params = com_x, com_y, theta of the motif
	const vector<double> motif_params = { params[0], params[1], params[2] };
	m_.setParameters( motif_params );

	// next params are for the IH tile
	double tile_params[tile_.numParameters()];
	for ( int i = 3; i < tile_.numParameters()+3; ++i ) {
		tile_params[i-3] = params[i];
	}
	tile_.setParameters( tile_params );

	// finally, params tell us where to start placing points on the edge curves
	edge_u0_ = params[3+tile_.numParameters()];

	// getParams updates the parameters internally
	getParameters();
}

const vector<double> Colloid::getParameters() {
	if (!motif_assigned_ || !tile_assigned_) {
		// Make sure the tile and motif have been assigned first.
		throw customException("must assign colloid and tile before parameters");
	}

	params_.clear();

	vector<double> dummy;
	dummy = m_.getParameters();
	for ( size_t i=0; i < dummy.size(); ++i ) {
		params_.push_back(dummy[i]);
	}
	
	double tile_dummy[tile_.numParameters()];
	tile_.getParameters(tile_dummy);
	for ( size_t i=0; i < tile_.numParameters(); ++i ) {
		params_.push_back(tile_dummy[i]);
	}

	params_.push_back(edge_u0_);

	return params_;
}

void Colloid::setMotif( Motif m ) {
	motif_assigned_ = true;
	m_ = m; // Create a copy
}

void Colloid::setTile( IsohedralTiling t ) {
	tile_assigned_ = true;
	tile_ = t; // Create a copy
}

bool Colloid::isMotifInside(const int N=20) // Put this many points on each edge to create a polygonal approximation of the tile's (smooth) boundary
{
	const double du = (1.0 - 0.0)/N;
	vector<vector<dvec2>> p = perimeter_(0.0, du, N-1);
	vector<dvec2> polygon;

	// Unroll points from each tile edge into single polygon vector
	polygon.clear();
	for( size_t i=0; i < p.size(); ++i ) {
		for( size_t j=0; j < p[i].size(); ++j ) {
			polygon.push_back(p[i][j]);
		}
	}

	// Check each point in motif 
	const vector<vector<double>> c = m_.getCoords();
	for( size_t i=0; i < c.size(); ++i ) {
		if (!pip(polygon, c[i])) {
			return false;
		}
	}

	return true;
}

double Colloid::tileArea()
{
	// Based on https://www.wikihow.com/Calculate-the-Area-of-a-Polygon
	// The corner control points are the vertices of the generating polygon.
	const int n = tile_control_points_.size();
	double sum1 = 0.0, sum2 = 0.0;
	for( int i=0; i < n; ++i ) {
		if ( i == n-1 ) {
			sum1 += tile_control_points_[i][0]*tile_control_points_[0][1];
			sum2 += tile_control_points_[i][1]*tile_control_points_[0][0];
		} else {
			sum1 += tile_control_points_[i][0]*tile_control_points_[i+1][1];
			sum2 += tile_control_points_[i][1]*tile_control_points_[i+1][0];
		}
	}
	return (sum1 - sum2)/2.;
}

void Colloid::buildBoundary_()
{
	vector<vector<dvec2>> edges = perimeter_(edge_u0_, edge_du_, 1+4+1);

	// Use a vector to hold the control points of the final tile outline.
	vector<dvec2> shape;

	// Iterate over the edges of a single tile, asking the tiling to
	// tell you about the geometric information needed to transform 
	// the edge shapes into position.  Note that this iteration is over
	// whole tiling edges.  It's also to iterator over partial edges
	// (i.e., halves of U and S edges) using t.parts() instead of t.shape().
	vector<int> identity;
	int total_edges = 0;
	for( auto i : tile_.shape() ) {
		// 1 + 4 colors + 1 stop codon
		vector<int> pattern = {
			0,
			i->getId()*4+1,
			i->getId()*4+2,
			i->getId()*4+3,
			i->getId()*4+4,
			0}; // 0 for stop codons
		
		// Interacting points on U and S edges are centro-symmetrically labelled
		if (i->getShape() == 1 || i->getShape() == 2) {
			pattern[3] = pattern[2];
			pattern[4] = pattern[1];
		}

		// Get the relevant edge shape created above using i->getId().
		const vector<dvec2> ed = edges[ i->getId() ];

		// Also get the transform that maps to the line joining consecutive
		// tiling vertices.
		const glm::dmat3& T = i->getTransform();

		// If i->isReversed() is true, we need to run the parameterization
		// of the path backwards.
		if( i->isReversed() ) {
			for( size_t idx = 1; idx < ed.size(); ++idx ) {
				shape.push_back( T * dvec3( ed[ed.size()-1-idx], 1.0 ) );
                        }
			for( size_t idx = 0; idx < pattern.size(); ++idx ) {
                                identity.push_back( pattern[pattern.size()-1-idx] );
                        }
			identity.push_back(-1); // -1 for control points
		} else {
			for( size_t idx = 1; idx < ed.size(); ++idx ) {
				shape.push_back( T * dvec3( ed[idx], 1.0 ) );
			}
			for( size_t idx = 0; idx < pattern.size(); ++idx ) {
                                identity.push_back( pattern[idx] );
                        }
			identity.push_back(-1); // -1 for control points
		}

		total_edges += 1;
	}

	// Count any gaps in the identities from S and U edges
	map<int, int> adjust;
	int ctr = 0;
	for (int i = 0; i <= *max_element(identity.begin(), identity.end()); ++i) {
		if (find(identity.begin(), identity.end(), i) == identity.end()) {
			++ctr; 
		} else {
			adjust[i] = ctr;
		}
	}

	boundary_ids_.clear();
	boundary_coords_.clear();
	tile_control_points_.clear();

	int skip = 0; 
	for ( size_t j = 0; j < shape.size() ; ++j) {
		vector<double> c = {shape[j].x, shape[j].y};
		if (identity[j] >= 0) { // Ignore control points
			if ( (identity[j] == 0) && (skip%2 == 0) ) { // Skip every other stop codon in fixed orientation
				boundary_coords_.push_back(c);
				boundary_ids_.push_back(identity[j]-adjust[identity[j]]);
				++skip;
			} else if (identity[j] == 0) {
				++skip;
			} else {
				boundary_coords_.push_back(c);
				boundary_ids_.push_back(identity[j]-adjust[identity[j]]);
			}
		} else {
			// Save control points to compute area, etc.
			tile_control_points_.push_back(c);
		}
	}
}

void Colloid::buildMotif_() 
{
	// Put motif inside of the boundary.
	return;
}

vector<vector<dvec2>> Colloid::perimeter_(double u0, double du, int n) {
	// Create a vector to hold some edge shapes.  The tiling tells you
	// how many distinct edge shapes you need, but doesn't know anything
	// about how those shapes might be represented.  It simply assumes
	// that each one will be a curve from (0,0) to (1,0).  The tiling
	// provides tools to let you map those curves into position around
	// the outline of a tile.  All the curves below have exactly four
	// control points.
	vector<vector<dvec2>> edges;

	// Generate edge shapes.
	for( U8 idx = 0; idx < tile_.numEdgeShapes(); ++idx ) {
		vector<dvec2> ej;

		// Define Bezier Curve that is sort of like a sphere impacting
		// This could definitely be changed, but will introduce more
		// free parameters.
		ej.push_back( dvec2( 0, 0 ) );
		ej.push_back( dvec2( 1/3., sphere_deform_ ));
		ej.push_back( dvec2( 2/3., sphere_deform_ ));
		ej.push_back( dvec2( 1, 0 ) );

		// Now, depending on the edge shape class, enforce symmetry 
		// constraints on edges.
		switch( tile_.getEdgeShape( idx ) ) {
		case J: 
			break;
		case U:
			ej[2].x = 1.0 - ej[1].x;
			ej[2].y = ej[1].y;
			break;
		case S:
			ej[2].x = 1.0 - ej[1].x;
			ej[2].y = -ej[1].y;
			break;
		case I:
			ej[1].y = 0.0;
			ej[2].y = 0.0;
			break;
		}

		// From Bezier curve determine location of 1+4+1 points along the
		// deformed edge. Here, i'll assume the stop codon is in contact
		// with the other 4. Later we will decide which stop codon to drop.
		vector<dvec2> coords;
		coords.push_back( dvec2( 0, 0 ) );
		for( int k = 0; k < n; ++k ) {
			coords.push_back(bezier(ej[0], ej[1], ej[2], ej[3], u0+k*du));	
		}
		coords.push_back( dvec2( 1, 0 ) );

		edges.push_back(coords);
	}

	return edges;
}

