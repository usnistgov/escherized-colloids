/*!
 * author: Nathan A. Mahynski
 *
 * Motif inside the tile.
 *
 * Described as a series of labelled points, akin to an XYZ file
 * or "molecule" without any bonds, etc.
 */

#ifndef MOTIF_H_
#define MOTIF_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include "utils.hpp"

using namespace std;

class Motif
{
public:
        Motif();
        ~Motif();

        void loadXYZ( const string filename ); // Load from XYZ file
        void dumpXYZ( const string filename ); // Print to XYZ file

        void setCoords( const vector<vector<double>> &coords );
	const vector<vector<double>> getCoords() { return coords_; }

        void setTypes( const vector<string> &types );
	const vector<string> getTypes() { return types_; }

	void setParameters ( const vector<double> &params );
        const vector<double> getParameters(); // Return COM and angle

private:
	void computeCOM_();
	void rotate_( const double theta );
        void translate_( const vector<double> &dx );

        vector<vector<double>> coords_;
        vector<string> types_;
        vector<double> com_;
	double theta_;
};

#endif
