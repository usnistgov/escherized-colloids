/**
 * @author Nathan A. Mahynski
 */

#include "motif.hpp"

Motif::Motif() {
}

Motif::~Motif() {
}

Motif::Motif( Motif& other ) {
	/**
	 * Instantiate a motif by copying from another.
	 *
	 * @param other Motif to copy.
	 */

	copy(other);
}

void Motif::copy( Motif& other ) {
	/**
	 * Copy the internal properties of another motif to this one.
	 *
	 * @params other Motif to copy.
	 */

	vector<double> p = other.getParameters();
	setCoords( other.getCoords(), p[ 2 ] );
	setTypes( other.getTypes() );
	//setParameters(other.getParameters()); // Handled by the above routines
}

void Motif::loadXYZ( const string filename ) {
	/**
	 * Load the motif coordinates from an xyz file.
	 * 
	 * @param filename Name of file to read from.
	 *
	 * @throws customException if file is corrupt or in an unexpected format.
	 */

	ifstream xyz( filename );
	string line, t;
	int ctr = 0, N = 0;
	double dummy;
	vector<vector<double>> coords;
	vector<string> types;
	
	try {
		while( getline( xyz, line ) )
		{
			istringstream iss( line );
			vector<double> c = { 0, 0 };
			if( ctr == 0 ) {
				iss >> N; // First line
			} else if (ctr == 1) {
				;
			} else {
				iss >> t >> c[0] >> c[1] >> dummy;
				types.push_back( t );
				coords.push_back( c );
			}
			ctr += 1;
		}
	} catch( const exception& e ) {
		throw( customException( "corrupt xyz file" ) );
	}

	if( ctr-2 != N ) {
		throw( customException( "line error in xyz file" ) );	
	}

	setCoords( coords, 0 );
	setTypes( types );

	return;      
}

void Motif::dumpXYZ( const string filename ) { 
	/**
	 * Dump the motif coordinates to an xyz file.
	 * 
	 * @param filename Name of file to dump to.
	 *
	 * @throws customException if unable to write or motif is missing information. 
	 */

	if( coords_.empty() ) {
		throw( customException( "coordinates have not been assigned" ) );
	}
	if( types_.empty() ) {
		throw( customException( "types have not been assigned" ) );
	}

	try {
		ofstream xyz( filename );

		xyz << coords_.size() << endl;
		xyz << endl;
		
		for( size_t i = 0; i < coords_.size(); ++i ) {
			xyz << types_[ i ] << "\t" << coords_[ i ][ 0 ] << "\t" << coords_[ i ][ 1 ] << "\t" << 0 << endl;
		}
	} catch( const exception& e ) {
		throw( customException( "unable to write to xyz file " ) );
	}

	return;
}

void Motif::setCoords( const vector<vector<double>> &coords, double theta=0.0 ) {
	/**
	 * Assign coordinates and orientation to the motif.
	 *
	 * The orientation provide is assumed to be an angle `theta` relative to
	 * some fixed reference frame.  When the motif is rotated, it is rotated
	 * to some absolute value of theta (effectively from (0, 2*pi]).
	 *
	 * @param coords Matrix of (x,y) coordinates.
	 * @param theta Absolute orientation to assume.
	 */

	coords_.clear();
	for( size_t i = 0; i < coords.size(); ++i ) {
		coords_.push_back( coords[ i ] );
	}

	// Compute COM so it is initialized
	computeCOM_();

	// Initial configuration defines initial orientiation
	theta_ = theta;

	return;
}
        
void Motif::setTypes( const vector<string> &types ) {
	/**
	 * Assign the chemical identities of the points in the motif.
	 * 
	 * @param types Vector of identities as strings.
	 */

	types_.clear();
	for( size_t i = 0; i < types.size(); ++i ) {
		types_.push_back( types[ i ] );	
	}

	return;
}

void Motif::setParameters( const vector<double> &params ) { 
	/**
	 * Assign motif parameters.
	 * 
	 * This is an unrolled list of parameters indicating the
	 * [com_x, com_y, theta] of the motif. The motif will be adjusted
	 * (translated and rotated) to have these values.
	 *
	 * @param params Vector of parameters described above. 
	 */

	if( coords_.empty() ) {
		throw( customException( "cannot assign motif parameters because it has not been assigned" ) );
	}

	//computeCOM_(); 
	const vector<double> shift = { -com_[0], -com_[1] }, new_com = { params[0], params[1] };

	// Shift to zero and rotate to theta
	translate( shift );
	rotate( params[ 2 ] );
	translate( new_com );

	// Recompute COM - should be same as params 
	//computeCOM_();
	com_[ 0 ] = params[ 0 ];
	com_[ 1 ] = params[ 1 ];
}

const vector<double> Motif::getParameters() { 
	/**
	 * Retrieve the parameters describing the motif, i.e., [com_x, com_y, theta].
	 * 
	 * @returns Vector of parameters described above.
	 */

	// Return COM and angle
	computeCOM_();
	vector<double> params = { com_[ 0 ], com_[ 1 ], theta_ };

	return params;
}

void Motif::computeCOM_() {
	/**
	 * Compute the center of mass of the motif and store it internally.
	 */

	vector<double> com = { 0.0, 0.0 };
	const size_t N = coords_.size();
	for( size_t i = 0; i < N; ++i ) {
		for( size_t j = 0; j < 2; ++j ) {
			com[ j ] += coords_[ i ][ j ]/N;
		}
		
	}

	com_ = com;
}

void Motif::rotate( const double theta ) {
	/**
	 * Rotate the motif to an absolute angle.
	 * 
	 * The angle, `theta`, is in an absolute reference frame.
	 * The motif is rotated by the difference between this `theta` and its
	 * current orientation angle. Positive angles are counterclockwise.
	 *
	 * @param theta Absolute angle to rotate to.
	 */

	const double dtheta = theta - theta_;
	double x = 0, y = 0;

	for( size_t i = 0; i < coords_.size(); ++i ) {
		x = coords_[ i ][ 0 ]*cos( dtheta ) - coords_[ i ][ 1 ]*sin( dtheta );
		y = coords_[ i ][ 0 ]*sin( dtheta ) + coords_[ i ][ 1 ]*cos( dtheta ); 
		coords_[ i ][ 0 ] = x;
		coords_[ i ][ 1 ] = y;
	}

	theta_ = theta;

	return;
}
 
void Motif::translate( const vector<double> &dx ) {
	/**
	 * Translate the motif by some amount.
	 * 
	 * Unlike rotation, `dx` indicates how much to translate relative to the current 
	 * position. The center of mass provides an absolute location.
	 *
	 * @param dx Vector of (dx, dy) to translate by.
	 */

	for( size_t i = 0; i < coords_.size(); ++i ) {
		for( size_t j = 0; j < 2; ++j ) {
			coords_[ i ][ j ] += dx[ j ];
		}
		
	}

	return;
}

void Motif::load( const string filename ) {
	/**
	 * Load a motif from a JSON file.
	 * 
	 * @param filename Name of file to read from.
	 *
	 * @throws customException if anything goes wrong.
	 */

	ifstream in( filename );
	json j;
	in >> j;

	try {
		double theta = j[ "theta" ].get<double>();
		vector<vector<double>> coords = j[ "coords" ].get<vector<vector<double>>>();
		vector<string> types = j[ "types" ].get<vector<string>>();
		setCoords( coords, theta );
		setTypes( types );	
		symmetry_ = j[ "symmetry" ].get<string>();
	} catch ( const exception& e ) {
		throw( customException( "unable to load Motif" ) );
	}

}

