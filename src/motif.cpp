/**
 * A motif for the colloid is what is "inside" the tile.
 * 
 * This file contains functions to build and manipulate motifs.
 * 
 * @author Nathan A. Mahynski
 */

#include "motif.hpp"

Motif::Motif() {
}

Motif::~Motif() {
}

Motif::Motif(Motif& other) {
	copy(other);
}

void Motif::copy(Motif& other) {
	vector<double> p = other.getParameters();
	setCoords(other.getCoords(), p[2]);
	setTypes(other.getTypes());
	//setParameters(other.getParameters()); // Handled by the above routines
}

void Motif::loadXYZ( const string filename ) {
	// Load from XYZ file	
	ifstream xyz(filename);
	string line, t;
	int ctr = 0, N = 0;
	vector<vector<double>> coords;
	vector<string> types;
	
	try {
		while (getline(xyz, line))
		{
			istringstream iss(line);
			vector<double> c = { 0, 0, 0 };
			if (ctr == 0) {
				iss >> N; // First line
			} else if (ctr == 1) {
				;
			} else {
				iss >> t >> c[0] >> c[1] >> c[2];
				types.push_back(t);
				coords.push_back(c);
			}
			ctr += 1;
		}
	} catch (const exception& e) {
		throw(customException("corrupt xyz file"));
	}

	if (ctr-2 != N) {
		throw(customException("line error in xyz file"));	
	}

	setCoords(coords, 0);
	setTypes(types);

	return;      
}

void Motif::dumpXYZ( const string filename ) { 
	// Print to XYZ file
	if (coords_.empty()) {
		throw(customException("coordinates have not been assigned"));
	}
	if (types_.empty()) {
		throw(customException("types have not been assigned"));
	}

	try {
		ofstream xyz(filename);

		xyz << coords_.size() << endl;
		xyz << endl;
		
		for( size_t i = 0; i < coords_.size(); ++i ) {
			xyz << types_[i] << "\t" << coords_[i][0] << "\t" << coords_[i][1] << "\t" << coords_[i][2] << endl;
		}
	} catch (const exception& e) {
		throw(customException("unable to write to xyz file"));
	}

	return;
}

void Motif::setCoords( const vector<vector<double>> &coords, double theta=0.0 ) {
	coords_.clear();
	for( size_t i = 0; i < coords.size(); ++i ) {
		coords_.push_back(coords[i]);
	}

	// Compute COM so it is initialized
	computeCOM_();

	// Initial configuration defines initial orientiation
	theta_ = theta;

	return;
}
        
void Motif::setTypes( const vector<string> &types ) {
	types_.clear();
	for ( size_t i = 0; i < types.size(); ++i ) {
		types_.push_back(types[i]);	
	}

	return;
}

void Motif::setParameters( const vector<double> &params ) { 
	if (coords_.empty()) {
		throw(customException("cannot assign motif parameters because it has not been assigned"));
	}

	//computeCOM_(); 
	const vector<double> shift = { -com_[0], -com_[1] }, new_com = { params[0], params[1]};

	// Shift to zero and rotate to theta
	translate(shift);
	rotate(params[2]);
	translate(new_com);

	// Recompute COM - should be same as params 
	//computeCOM_();
	com_[0] = params[0];
	com_[1] = params[1];
}

const vector<double> Motif::getParameters() { 
	// Return COM and angle
	computeCOM_();
	vector<double> params = { com_[0], com_[1], theta_ };

	return params;
}

void Motif::computeCOM_() {
	vector<double> com = {0.0, 0.0};
	const size_t N = coords_.size();
	for( size_t i = 0; i < N; ++i ) {
		for( size_t j = 0; j < 2; ++j ) {
			com[j] += coords_[i][j]/N;
		}
		
	}

	com_ = com;
}

void Motif::rotate( const double theta ) {
	// Assume theta (cc) is absolute relative to initialization
	const double dtheta = theta - theta_;
	double x = 0, y = 0;

	for( size_t i = 0; i < coords_.size(); ++i ) {
		x = coords_[i][0]*cos(dtheta) - coords_[i][1]*sin(dtheta);
		y = coords_[i][0]*sin(dtheta) + coords_[i][1]*cos(dtheta); 
		coords_[i][0] = x;
		coords_[i][1] = y;
	}

	theta_ = theta;

	return;
}
 
void Motif::translate( const vector<double> &dx ) {
	// Translate coordinates
	for( size_t i = 0; i < coords_.size(); ++i ) {
		for( size_t j = 0; j < 2; ++j ) {
			coords_[i][j] += dx[j];
		}
		
	}

	return;
}

