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
	
	Motif(Motif& other); 
	void copy(Motif& other);

        void loadXYZ( const string filename ); // Load from XYZ file
        void dumpXYZ( const string filename ); // Print to XYZ file

        void setCoords( const vector<vector<double>> &coords, double theta );
	const vector<vector<double>> getCoords() const { return coords_; }

        void setTypes( const vector<string> &types );
	const vector<string> getTypes() const { return types_; }

	void setParameters ( const vector<double> &params );
        const vector<double> getParameters(); // Return COM and angle

	vector<double> getCOM() { computeCOM_(); return com_; }
	void rotate( const double theta );
        void translate( const vector<double> &dx );

private:
	void computeCOM_();

        vector<vector<double>> coords_;
        vector<string> types_;
        vector<double> com_;
	double theta_;
};

#endif
