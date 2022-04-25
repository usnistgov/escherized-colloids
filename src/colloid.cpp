/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 */

#include "src/colloid.hpp"

double mod(double x, double v) {
  /**
   * Modulus operation for negative values and floating point numbers.
   */
  double y = x;
  while (y < 0) {
    y += v;
  }
  while (y > v) {
    y -= v;
  }

  return y;
}

dvec2 bezier(dvec2 p0, dvec2 p1, dvec2 p2, dvec2 p3, double u) {
  /**
   * Compute parameterized location along a Bezier curve.
   *
   * In particular, this code is based on the discussion found at
   * https://www.geeksforgeeks.org/cubic-bezier-curve-implementation-in-c/
   *
   * @param p0 First Bezier control point.
   * @param p1 Second Bezier control point.
   * @param p2 Third Bezier control point.
   * @param p3 Fourth Bezier control point.
   * @param u Position along the curve; 0 is the start of the curve, 1 is the
   * end.
   */

  dvec2 coord;
  coord.x = pow(1.0 - u, 3) * p0.x + 3.0 * u * pow(1.0 - u, 2) * p1.x +
            3.0 * (1.0 - u) * pow(u, 2) * p2.x + pow(u, 3) * p3.x;
  coord.y = pow(1.0 - u, 3) * p0.y + 3.0 * u * pow(1.0 - u, 2) * p1.y +
            3.0 * (1.0 - u) * pow(u, 2) * p2.y + pow(u, 3) * p3.y;

  return coord;
}

vector<vector<double>> bbox(const vector<vector<double>>& points) {
  /**
   * Compute the bounding box around a series of points. By convention,
   * the (min_x, min_y) coordinate is returned first.
   *
   * @param points All (x,y) coordinates to consider.
   *
   * @returns Bounding box coordinates.
   */

  vector<double> max = {points[0][0], points[0][1]},
                 min = {points[0][0], points[0][1]};

  for (size_t i = 0; i < points.size(); ++i) {
    for (size_t j = 0; j < 2; ++j) {
      if (points[i][j] < min[j]) {
        min[j] = points[i][j];
      }
      if (points[i][j] > max[j]) {
        max[j] = points[i][j];
      }
    }
  }

  vector<vector<double>> box = {
      {min[0], min[1]}, {max[0], min[1]}, {min[0], max[1]}, {max[0], max[1]}};

  return box;
}

Colloid::Colloid() : tile_(IsohedralTiling(1)) {
  /**
   * Instantiate a new colloid.
   *
   * The Tactile library requires an isohedral tile number be given at
   * instantiation so the colloid is initialized, by default, to have tile IH1.
   * This can be changed later. Also by default, sphere_deform_ = 0.25, edge_du_
   * = 0.1.
   *
   * You should later set relevant parameters and call `init()` to initialize
   * the colloid.
   */

  defaults_();
}

Colloid::Colloid(Motif m, IsohedralTiling t, vector<double> tile_u0, vector<double> edge_df, const double du, const bool debug)
    : tile_(IsohedralTiling(1)) {
  /**
   * Instantiate a new colloid if all parameters are known (preferred method).
   *
   * The Tactile library requires an isohedral tile number be given at
   * instantiation so the colloid is initialized, by default, to have tile IH1.
   * This can be changed later.
   */

  // If you create a Colloid object manually, follow these steps to correctly initialize it.
  defaults_();
  setU0(tile_u0);
  setDform(edge_df);
  setDU(du);

  setMotif(m);
  setTile(t);
  init(debug); // This will find a good value for tile_scale_
}

const vector<double> Colloid::getU0() {
  if (edge_u0_.size() == 0) {
    throw customException("u0 has not been set");
  }
  return edge_u0_;
}

const vector<double> Colloid::getDform() {
  if (edge_df_.size() == 0) {
    throw customException("df has not been set");
  }
  return edge_df_;
}

const double Colloid::getDU() {
  if (edge_du_ < 0) {
    throw customException("du has not been set");
  }
  return edge_du_;
}

void Colloid::init(bool debug) {
  buildBoundary_();
  initMotif_(10.0, 0.1, 1000, 20, debug);
  built_ = true;
} // Initialize the colloid.

void Colloid::defaults_() {
  /**
   * Assign default paramter values.
   */

  tile_assigned_ = false;
  motif_assigned_ = false;
  built_ = false;
  setTileScale(1.0);
  setU0({}); // size of 0 serves to indicate it has not been set yet
  setDform({}); // size of 0 serves to indicate it has not been set yet
  setDU(-1.0); // <0 serves as flag to indicate it has not been set yet
}

Colloid::~Colloid() {}

vector<double> Colloid::boundaryCOM() {
  /**
   * Compute the center of mass (COM) of the boundary points.
   *
   * The boundary points include the functional points and the
   * "stop codons".  No control points or other parts of the
   * edges used by Tactile are used here.
   *
   * @returns com (x,y) coordinates of the COM.
   *
   * @throws customException if boundary coordinates have not been initialized.
   */

  if (!boundary_coords_.empty()) {
    vector<double> com = {0.0, 0.0};
    const size_t N = boundary_coords_.size();
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < 2; ++j) {
        com[j] += boundary_coords_[i][j] / N;
      }
    }

    return com;
  } else {
    throw(customException("boundary has not been defined yet"));
  }
}

vector<double> Colloid::unscale_coords_(const vector<double>& scaled_coords) {
  /**
   * Unscale the motif's coordinates.
   */

  vector<vector<double>> box = bbox(tile_control_points_);
  const double wx = box[1][0] - box[0][0], wy = box[2][1] - box[1][1];
  const vector<double> us = {box[0][0] + scaled_coords[0] * wx,
                             box[0][1] + scaled_coords[1] * wy};

  return us;
}

vector<double> Colloid::scale_coords_(const vector<double>& unscaled_coords) {
  /**
   * Scale the motif's coordinates.
   */

  vector<vector<double>> box = bbox(tile_control_points_);
  const double wx = box[1][0] - box[0][0], wy = box[2][1] - box[1][1];
  const vector<double> scaled = {(unscaled_coords[0] - box[0][0]) / wx,
                                 (unscaled_coords[1] - box[0][1]) / wy};

  return scaled;
}

void Colloid::setParameters(const vector<double>& params, const double df_min) {
  /**
   * Assign all parameters defining the colloid.
   *
   * The parameters are an unrolled (double) vector useful for optimization
   * schemes. They are: [{motif scaled_com_x, motif scaled_com_y}, motif_theta,
   * {v0, v1, etc. for tile}, {edge_u0}, {edge_df}, tile_scale]. The scaled com's are values
   * [0, 1] which describe the location in terms of the bounding box around the
   * tile. Motif COM coordinates are only used for tiles which are fundamental
   * domains.
   *
   * The motif will track its coordinates in absolute units, but the colloid
   * uses reduced units since it knows about the tile and means that these
   * variables can be given reasonable bounds to an optimizer in advance.
   *
   * The following checks are performed and will throw an exception if unmet:
   * 1. "Patch" of decorations falls on the Bezier curve between [0, 1).
   * 2. |df| >= df_min. This is to enforce that edges SHOULD have curvature.
   * 3. tile_scale > 0.
   * 4. Tile should not be self-intersecting.
   *
   * Note: A colloid can be instantiated with such parameters (which can be helpful
   * for initialization, e.g., with straight edges so df = 0), but setParameters()
   * is used during optimizations so this is where we enforce such considerations.
   *
   * @param params Parameter vector described above.
   * @param df_min Minimum magnitude of df allowed for an edge. 
   *
   * @throws customException if tile or motif has not been assigned yet or
   * any of the parameters are invalid (e.g., |df| < df_min, or tile_scale < 0).
   */

  if (!motif_assigned_ ||
      !tile_assigned_) { // Makes sure the tile and colloid have been assigned
    throw customException(
        "must assign tile and motif before setting new parameters");
  }

  // Make sure we have the right number of parameters in the vector
  if (params.size() != static_cast<unsigned int>(3 + tile_.numParameters() + 2*tile_.numEdgeShapes() + 1)) {
    throw customException(
        "incorrect number of parameters provided");
  }

  // Tile
  double tile_params[tile_.numParameters()];
  for (int i = 3; i < tile_.numParameters() + 3; ++i) {
    tile_params[i - 3] = params[i];
  }
  tile_.setParameters(tile_params);

  // A few validity checks:

  // 1a. Make sure u0 + du*(N-1) < 1-du/2, i.e., when you place the decorations they
  // do not go past the "end" of the bezier curve.  For S and U edges, I have
  // doubled du (since each Bezier curve is a "part", or half, of a tile edge)
  // to more consistently space points with respect to I and J edges.

  // 1b. For consistency, we also check that u0 > du/2. "du" may be thought of as the
  // diameter of the circular decoration.

  // The net valid range is du/2 <= u0 <= 1-5*du-du/2; e.g., 0.05 < u0 < 0.45 if du = 0.1.

  // In practice we place 1+4+1 total points "starting" from u0 along Bezier curve and delete
  // one of the extra stop codons later once a fixed orientation for them is determined.
  // If the "furthest" one is kept then it is at u0 + (1+4+1 - 1)*du which we need to test,
  // since at this point we do not know if we are going to keep it or not.
  vector<double> u0_params(params.begin()+3+tile_.numParameters(), params.begin()+3+tile_.numParameters()+tile_.numEdgeShapes());
  const int n_decor = 1+4+1;
  for (U8 idx = 0; idx < tile_.numEdgeShapes(); ++idx) {
    int n_edge = 0;
    stringstream ss;
    switch (tile_.getEdgeShape(idx)) {
      case J: case I:
        n_edge = n_decor;
        if ( !(u0_params[idx] + getDU()*(n_edge-1) <= 1.0 - getDU()/2.0) ) {
          ss << "decorations on (J/I) edge index " << int(idx) << " will go past the end, u0 = " << u0_params[idx] << ", du = " << getDU();
          throw customException(ss.str());
        }
        if ( !(u0_params[idx] >= getDU()/2.0) ) {
          ss << "decorations on (J/I) edge index " << int(idx) << " will begin before the start, u0 = " << u0_params[idx] << ", du = " << getDU();
          throw customException(ss.str());
        }
        break;
      case U: case S: // In perimeter_edges_() we doubled dU since edge is half the size
        // n_edge <= n_decor/2 due to integer rounding (intentional)
        n_edge = (n_decor-2)/2+1; // This is the n used in perimeter_edges_()
        // In these cases we "start" from the midpoint at "1" not "0" but the direction does
        // not matter, just the total "length" of the patch.
        if ( !(u0_params[idx] + 2*getDU()*(n_edge-1) <= 1.0 - getDU()) ) {
          ss << "decorations on (S/U) edge index " << int(idx) << " will go past the end, u0 = " << u0_params[idx] << ", du = " << getDU();
          throw customException(ss.str());
        }
        if ( !(u0_params[idx] >= getDU()) ) {
          ss << "decorations on (S/U) edge index " << int(idx) << " will begin before the start, u0 = " << u0_params[idx] << ", du = " << getDU();
          throw customException(ss.str());
        }
        break;
    }
  }
  setU0(u0_params);

  // 2. Check |df| > df_min
  // To practically allow the optimizer to change signs, just set a gap near 0 which is "off limits"
  vector<double> df_params(params.begin()+3+tile_.numParameters()+tile_.numEdgeShapes(), params.begin()+3+tile_.numParameters()+2*tile_.numEdgeShapes());
  for (unsigned int i=0; i < df_params.size(); ++i) {
    if (fabs(df_params[i]) < fabs(df_min)) {
      throw customException("|df| is too small");
    }
  }
  setDform(df_params);

  // 3. Check scale is > 0
  if (params[params.size()-1] <= 0.0) {
    throw customException("tile scale must be positive");
  }
  setTileScale(params[params.size()-1]);

  // 4. Check if the tile is self-intersecting
  const int intersections = countIntersections(10);
  if (intersections > 0) {
    stringstream ss;
    ss << "tile has " << intersections << " self-intersections";
    throw customException(ss.str());
  }

  // Build boundary (need updated tile_control_points_) before computing scaled
  // coordinates
  buildBoundary_();

  // Motif - convert scaled to absolute coordinates
  // Setting the tile parameters, etc. changes the motif's absolute location.
  // (1) Already set new scale, (2) unscale to get absolute coordinates, then
  // (3) assign the motif position and orientation, but (4) modify them as
  // needed to conform to the desired symmetry.  Note that this means the
  // params passed to setParameters() are NOT necessarily what the true
  // parameters are; however, the input is uniquely mapped to some new value(s)
  // so this is reproducible.  Just be aware that the input from say an optimizer
  // and the resulting colloid might be different because of this change.
  // Since getParameters() collects information from the motif directly this
  // should work just fine.
  const vector<double> scaled_coords = {params[0], params[1]};
  const vector<double> us = unscale_coords_(scaled_coords);
  vector<double> motif_params = {us[0], us[1], params[2]};
  if (!isTileFundamental()) {
    constrain_(&motif_params);
  }
  m_.setParameters(motif_params);
}

void Colloid::constrain_(vector<double>* motif_params) {
  /**
   * Some motif parameters are initially proposed, then updated by constraints
   * imposed by the nature of the tile. The vector will be modified so the
   * motif can be updated with those values elsewhere. The tile and motif are
   * not modified here.
   *
   * buildBoundary_() should be called first to establish the control points,
   * since these are the vertices of the polygon (in unscaled units) used to
   * compute locations and orientations of motifs.
   *
   * @param[in, out] motif_params Motif parameters.
   */

  const int ih_number = static_cast<int>(tile_.getTilingType());

  vector<double> orig_coords = {motif_params->at(0), motif_params->at(1)},
                 results, p0(2, 0.0), p1(2, 0.0);
  int induced = 0;
  string prefix = "";

  const double jitter = 1.0e-12; // Jitter is used to help numerically protect PIP routine

  // Conventions: p0 -> p1 forms a minimum counterclockwise angle wrt x-axis.

  if (ih_number == 64) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[3][0] - tile_control_points_[2][0]) / 2.0 +
           tile_control_points_[2][0] + jitter;
    p0[1] = (tile_control_points_[3][1] + tile_control_points_[2][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[1][0] - tile_control_points_[0][0]) / 2.0 +
            tile_control_points_[0][0] + jitter;
    p1[1] = (tile_control_points_[1][1] + tile_control_points_[0][1]) / 2.0 + jitter;
  } else if (ih_number == 12) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[4][0] - tile_control_points_[5][0]) / 2.0 +
           tile_control_points_[5][0] + jitter;
    p0[1] = (tile_control_points_[4][1] + tile_control_points_[5][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[2][0] - tile_control_points_[1][0]) / 2.0 +
            tile_control_points_[1][0] + jitter;
    p1[1] = (tile_control_points_[2][1] + tile_control_points_[1][1]) / 2.0 + jitter;
  } else if (ih_number == 14) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) horizontal
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[2][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[0][1] - tile_control_points_[2][1]) / 2.0 +
           tile_control_points_[2][1] + jitter;

    p1[0] = (tile_control_points_[3][0] + tile_control_points_[5][0]) / 2.0 + jitter;
    p1[1] = (tile_control_points_[5][1] - tile_control_points_[3][1]) / 2.0 +
            tile_control_points_[3][1] + jitter;
  } else if (ih_number == 68) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) horizontal
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = tile_control_points_[0][0] + jitter;
    p0[1] = tile_control_points_[0][1] + jitter;

    p1[0] = tile_control_points_[2][0] + jitter;
    p1[1] = tile_control_points_[2][1] + jitter;
  } else if (ih_number == 13) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[4][0] - tile_control_points_[5][0]) / 2.0 +
           tile_control_points_[5][0] + jitter;
    p0[1] = (tile_control_points_[4][1] + tile_control_points_[5][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[2][0] - tile_control_points_[1][0]) / 2.0 +
            tile_control_points_[1][0] + jitter;
    p1[1] = (tile_control_points_[2][1] + tile_control_points_[1][1]) / 2.0 + jitter;
  } else if (ih_number == 15) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) horizontal
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[3][0] + tile_control_points_[5][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[3][1] - tile_control_points_[5][1]) / 2.0 +
           tile_control_points_[5][1] + jitter;

    p1[0] = (tile_control_points_[2][0] + tile_control_points_[0][0]) / 2.0 + jitter;
    p1[1] = (tile_control_points_[2][1] - tile_control_points_[0][1]) / 2.0 +
            tile_control_points_[0][1] + jitter;
  } else if (ih_number == 66) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[3][0] - tile_control_points_[2][0]) / 2.0 +
           tile_control_points_[2][0] + jitter;
    p0[1] = (tile_control_points_[3][1] + tile_control_points_[2][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[1][0] - tile_control_points_[0][0]) / 2.0 +
            tile_control_points_[0][0] + jitter;
    p1[1] = (tile_control_points_[1][1] + tile_control_points_[0][1]) / 2.0 + jitter;
  } else if (ih_number == 69) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) horizontal
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = tile_control_points_[1][0] + jitter;
    p0[1] = tile_control_points_[1][1] + jitter;

    p1[0] = tile_control_points_[3][0] + jitter;
    p1[1] = tile_control_points_[3][1] + jitter;
  } else if (ih_number == 26) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) horizontal
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[1][0] + tile_control_points_[2][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[1][1] - tile_control_points_[2][1]) / 2.0 +
           tile_control_points_[2][1] + jitter;

    p1[0] = tile_control_points_[4][0] + jitter;
    p1[1] = tile_control_points_[4][1] + jitter;
  } else if (ih_number == 67) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[2][0] - tile_control_points_[1][0]) / 2.0 +
           tile_control_points_[1][0] + jitter;
    p0[1] = (tile_control_points_[2][1] + tile_control_points_[1][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[0][0] - tile_control_points_[3][0]) / 2.0 +
            tile_control_points_[3][0] + jitter;
    p1[1] = (tile_control_points_[0][1] + tile_control_points_[3][1]) / 2.0 + jitter;
  } else if (ih_number == 91) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[1][0] - tile_control_points_[0][0]) / 2.0 +
           tile_control_points_[0][0] + jitter;
    p0[1] = (tile_control_points_[1][1] + tile_control_points_[0][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[2][0] + jitter;
    p1[1] = tile_control_points_[2][1] + jitter;
  } else if (ih_number == 16) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = tile_control_points_[4][0] + jitter;
    p0[1] = tile_control_points_[4][1] + jitter;

    p1[0] = tile_control_points_[1][0] + jitter;
    p1[1] = tile_control_points_[1][1] + jitter;
  } else if (ih_number == 36) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = tile_control_points_[1][0] + jitter;
    p0[1] = tile_control_points_[1][1] + jitter;

    p1[0] = tile_control_points_[3][0] + jitter;
    p1[1] = tile_control_points_[3][1] + jitter;
  } else if (ih_number == 29) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[4][0] - tile_control_points_[3][0]) / 2.0 +
           tile_control_points_[3][0] + jitter;
    p0[1] = (tile_control_points_[4][1] + tile_control_points_[3][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[1][0] + jitter;
    p1[1] = tile_control_points_[1][1] + jitter;
  } else if (ih_number == 71) {
    // Mirror line along diagonal of tile. Tactile has it as a (fixed) 135 deg
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = tile_control_points_[2][0] + jitter;
    p0[1] = tile_control_points_[2][1] + jitter;

    p1[0] = tile_control_points_[0][0] + jitter;
    p1[1] = tile_control_points_[0][1] + jitter;
  } else if (ih_number == 82) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[0][0] - tile_control_points_[2][0]) / 2.0 +
           tile_control_points_[2][0] + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[2][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[1][0] + jitter;
    p1[1] = tile_control_points_[1][1] + jitter;
  } else if (ih_number == 32) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) horizontal
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = tile_control_points_[0][0] + jitter;
    p0[1] = tile_control_points_[0][1] + jitter;

    p1[0] = tile_control_points_[2][0] + jitter;
    p1[1] = tile_control_points_[2][1] + jitter;
  } else if (ih_number == 40) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[0][0] - tile_control_points_[2][0]) / 2.0 +
           tile_control_points_[2][0] + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[2][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[1][0] + jitter;
    p1[1] = tile_control_points_[1][1] + jitter;
  } else if (ih_number == 72) {
    // Mirror lines cross in the middle of tile. Tactile has it as a (fixed) vertical
    // and horizontal lines. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 2;  // S(P|M) = d2

    p0[0] = (tile_control_points_[0][0] - tile_control_points_[1][0]) / 2.0 +
           tile_control_points_[1][0] + jitter;
    p0[1] = (tile_control_points_[0][1] - tile_control_points_[3][1]) / 2.0 +
           tile_control_points_[3][1] + jitter;

    p1[0] = (tile_control_points_[0][0] + tile_control_points_[3][0]) / 2.0 + jitter;
    p1[1] = p0[1];
  } else if (ih_number == 17) {
    // Mirror lines cross in the middle of tile. Tactile has it as a (fixed) vertical
    // and horizontal lines. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 2;  // S(P|M) = d2

    p0[0] = (tile_control_points_[2][0] - tile_control_points_[1][0]) / 2.0 +
           tile_control_points_[1][0] + jitter;
    p0[1] = (tile_control_points_[5][1] - tile_control_points_[1][1]) / 2.0 +
           tile_control_points_[1][1] + jitter;

    p1[0] = tile_control_points_[0][0] + jitter;
    p1[1] = tile_control_points_[0][1] + jitter;
  } else if (ih_number == 74) {
    // Mirror lines cross in the middle of tile. Tactile has it as a (fixed) vertical
    // and horizontal lines. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 2;  // S(P|M) = d2

    p0[0] = (tile_control_points_[3][0] + tile_control_points_[1][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[2][1] + tile_control_points_[0][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[2][0] + jitter;
    p1[1] = tile_control_points_[2][1] + jitter;
  } else if (ih_number == 73) {
    // Mirror lines cross in the middle of tile. Tactile has it as a (fixed) vertical
    // and horizontal lines. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 2;  // S(P|M) = d2

    p0[0] = (tile_control_points_[2][0] + tile_control_points_[1][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[1][1] + tile_control_points_[0][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[2][0] + tile_control_points_[3][0]) / 2.0 + jitter;
    p1[1] = (tile_control_points_[2][1] + tile_control_points_[3][1]) / 2.0 + jitter;
  } else if (ih_number == 37) {
    // Mirror lines cross in the middle of tile. Tactile has it as a (fixed) vertical
    // and horizontal lines. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 2;  // S(P|M) = d2

    p0[0] = (tile_control_points_[3][0] + tile_control_points_[1][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[3][1] + tile_control_points_[1][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[2][0] + jitter;
    p1[1] = tile_control_points_[2][1] + jitter;
  } else if (ih_number == 18) {
    // Mirror lines cross in the middle of tile. Tactile has them as fixed lines.
    // Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 3;  // S(P|M) = d3

    p0[0] = (tile_control_points_[1][0] + tile_control_points_[4][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[1][1] + tile_control_points_[4][1]) / 2.0 + jitter;

    p1[0] = (tile_control_points_[5][0] + tile_control_points_[4][0]) / 2.0 + jitter;
    p1[1] = (tile_control_points_[5][1] + tile_control_points_[4][1]) / 2.0 + jitter;
  } else if (ih_number == 93) {
    // Mirror lines cross in the middle of tile. Tactile has them as fixed lines.
    // Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 3;  // S(P|M) = d3

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0]) / 3.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1]) / 3.0 + jitter;

    p1[0] = tile_control_points_[1][0] + jitter;
    p1[1] = tile_control_points_[1][1] + jitter;
  } else if (ih_number == 76) {
    // Mirror lines cross in the middle of tile. Tactile has them as fixed lines.
    // Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 4;  // S(P|M) = d4

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    p1[0] = (tile_control_points_[2][0] + tile_control_points_[3][0])/2. + jitter;
    p1[1] = (tile_control_points_[2][1] + tile_control_points_[3][1])/2. + jitter;
  } else if (ih_number == 20) {
    // Mirror lines cross in the middle of tile. Tactile has them as fixed lines.
    // Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 6;  // S(P|M) = d6

    p0[0] = (tile_control_points_[1][0] + tile_control_points_[4][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[1][1] + tile_control_points_[4][1]) / 2.0 + jitter;

    p1[0] = tile_control_points_[4][0] + jitter;
    p1[1] = tile_control_points_[4][1] + jitter;
  } else if (ih_number == 8) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[5][0] + tile_control_points_[4][0] + tile_control_points_[2][0] + tile_control_points_[1][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[5][1] + tile_control_points_[4][1] + tile_control_points_[2][1] + tile_control_points_[1][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 57) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 9) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[4][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[5][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[4][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[5][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 59) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[2][0]) / 2.0 + jitter;
    p0[1] = (tile_control_points_[1][1] + tile_control_points_[3][1]) / 2.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 58) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 61) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 34) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 2;  // S(P|M) = c2

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 10) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 3;  // S(P|M) = c3

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[5][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[5][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 90) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 3;  // S(P|M) = c3

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0]) / 3.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1]) / 3.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 62) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 4;  // S(P|M) = c4

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[1][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[1][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else if (ih_number == 11) {
    // Rotation center in the middle of tile.
    prefix = "c";
    induced = 6;  // S(P|M) = c6

    p0[0] = (tile_control_points_[0][0] + tile_control_points_[5][0] + tile_control_points_[2][0] + tile_control_points_[3][0]) / 4.0 + jitter;
    p0[1] = (tile_control_points_[0][1] + tile_control_points_[5][1] + tile_control_points_[2][1] + tile_control_points_[3][1]) / 4.0 + jitter;

    // p1 is ignored
    p1[0] = p0[0];
    p1[1] = p0[1];
  } else {
    throw(customException("unrecognized tile type"));
  }

  // Perform the revision
  results = revise_motif_params_(p0, p1, orig_coords, motif_params->at(2), prefix, induced);

  // 3. Update coordinates and orientation
  motif_params->at(0) = results[0];
  motif_params->at(1) = results[1];
  motif_params->at(2) = results[2];
}

const vector<double> Colloid::revise_motif_params_(const vector<double>& p0,
                                      const vector<double>& p1,
                                      const vector<double>& orig_coords,
                                      const double current_theta,
                                      const string prefix, const int induced) {
  /**
   * Strategy:
   * Induced = c(n>1), place motif at rotation center, no forced rotation
   * necessary Induced = d(n>1), place motif at mirror intersection Induced =
   * d(n=1), place motif along the mirror line
   *
   * Proposal:
   * For d(n>1) could enforce that p0 should be provided such that this should
   * be the motif COM. Similarly, this convention could be helpful for c(n)
   * where p0 should provide the motif COM.
   *
   * p0 and p1 are on a mirror line; if d1 they are on the only line; if d(n>1)
   * p0 should be the mirror intersection point and p1 is a point on some mirror
   * that defines a line; if c(n > 1) p0 represents the rotation center and 
   * p1 is ignored.
   *
   * If the tile is incommensurate with the motif symmetry an exception is thrown;
   * For example, if IH64 (where S(P|M) = d1) and the tile has c(n) symmetry. 
   */
  vector<double> projected_coords(2, 0);
  double absolute_theta = 0.0;

  const int n = compatibility_check(m_, prefix, induced);

  if (prefix.compare("d") == 0) {
    // 1. Put motif COM on mirror lines
    if (induced == 1) {
      // If only 1 mirror line we have a DoF in terms of where on that line.
      projected_coords = project_to_line(p0, p1, orig_coords);
    } else {
      // With multiple mirrors, their intersection defines the COM of the motif.
      // By convention, we have chosen p0 to be that intersection.
      projected_coords = p0;
    }

    // 2. Round orientation to nearest allowable absolute theta value
    // Simply allow rotations that make mirror lines coincide (mod 2 leaves
    // motif unchanged)
    const double curr_theta = thetaBounds(current_theta); // 0 <= t < 2PI
    double min_diff = pow(2.0 * M_PI, 2), diff = 0.0,
           allowed_angle = 0.0;
    for (int i = 0; i <= 2 * n; // when d_inf occurs searches each integer degree
         ++i) {
      allowed_angle = i * M_PI / n + tile_mirror_alignment(p0, p1);

      double min_t = allowed_angle - curr_theta;
      if (fabs(min_t - 2*M_PI) < fabs(min_t)) {
        // fabs(allowed_angle - (curr_theta+2*M_PI) < fabs(min_t)
        min_t -= 2*M_PI;
      }

      diff = pow(min_t, 2);
      if (diff < min_diff) {
        min_diff = diff;
        absolute_theta = thetaBounds(allowed_angle);
      }
    }
  } else {
    // 1. Assign motif COM to rotation center
    projected_coords = p0;

    // 2. No rotation is required
    absolute_theta = thetaBounds(current_theta);
  }

  // 3. Update and return
  vector<double> results(3, 0.0);
  results[0] = projected_coords[0];
  results[1] = projected_coords[1];
  results[2] = absolute_theta;

  return results;
}

const vector<double> Colloid::getParameters() {
  /**
   * Retrieve all parameters defining the colloid. Also recomputes it
   * internally.
   *
   * The parameters are an unrolled (double) vector useful for optimization
   * schemes. They are: [ {motif scaled_com_x, motif scaled_com_y}, motif_theta,
   * {v0, v1, etc. for tile}, {edge_u0}, {edge_df}, tile_scale]. See `setParameters()` for an
   * explanation of the scaled coordinates. Motif COM coordinates are only
   * included for tiles which are fundamental domains.
   *
   * @returns Parameter vector.
   *
   * @throws customException if tile or motif has not been assigned yet.
   */

  if (!motif_assigned_ || !tile_assigned_) {
    throw customException(
        "must assign assign tile and motif before getting parameters");
  }

  params_.clear();

  vector<double> dummy, rescaled;
  dummy = m_.getParameters();
  for (size_t i = 0; i < 2; ++i) {
    params_.push_back(dummy[i]); // Unscaled COM coordinates
  }
  rescaled = scale_coords_(params_);
  params_[0] = rescaled[0]; // scaled_com_x
  params_[1] = rescaled[1]; // scaled_com_y
  params_.push_back(dummy[2]); // theta

  double tile_dummy[tile_.numParameters()];

  tile_.getParameters(tile_dummy);
  for (size_t i = 0; i < tile_.numParameters(); ++i) {
    params_.push_back(tile_dummy[i]);
  }
  for (size_t i = 0; i < edge_u0_.size(); ++i) {
    params_.push_back(edge_u0_[i]);
  }
  for (size_t i = 0; i < edge_df_.size(); ++i) {
    params_.push_back(edge_df_[i]);
  }
  params_.push_back(getTileScale());

  return params_;
}

void Colloid::setMotif(Motif m) {
  /**
   * Assign the motif by copying an existing one.
   *
   * The motif's initial orientation should be consistent with any symmetry
   * constraints.  It's initial orientation is always used as a reference.
   *
   * @param m Motif for the colloid to use.
   */

  motif_assigned_ = true;
  m_.copy(m);  // Create a copy - motif may not have been initialized yet
}

Motif Colloid::getMotif() {
  /**
   * Retieve the colloid's motif.
   *
   * @returns The colloid's motif.
   *
   * @throws customException if the motif has not been assigned yet.
   */

  if (!motif_assigned_) {
    throw(customException("motif has not been assigned yet"));
  } else {
    return m_;
  }
}

void Colloid::setTile(IsohedralTiling t) {
  /**
   * Assign the tile.
   *
   * @param t Isohedral tile from Tactile library.
   */

  tile_assigned_ = true;
  tile_ = t;  // Create a copy - already initialized at instantiation
}

const IsohedralTiling Colloid::getTile() {
  /**
   * Retieve the colloid's tile.
   *
   * @returns The colloid's tile.
   *
   * @throws customException if the tile has not been assigned yet.
   */

  if (!tile_assigned_) {
    throw(customException("tile has not been assigned yet"));
  } else {
    return tile_;
  }
}

const int Colloid::countIntersections(const int N) {
  /**
   * Approximate the tile's boundary as a polygon with a discrete
   * number of points on each edge and count the number of interecting
   * line segments.
   *
   * @param Number of points to discretize each edge into.
   *
   * @returns The total number of intersecting pairs of line segments.
   */

  // Build the perimeter
  vector<vector<dvec2>> polygon = buildTilePolygon(N);

  // Compare each line segment to all other line segments
  dvec2 p1, q1, p2, q2;
  int total = 0;
  for (unsigned int edge_idx=0; edge_idx < polygon.size(); ++edge_idx) {
    for (unsigned int i=0; i < polygon[edge_idx].size()-1; ++i) {
      // First segment
      p1 = polygon[edge_idx][i];
      q1 = polygon[edge_idx][i+1];

      for (unsigned int comp_idx=edge_idx; comp_idx < polygon.size(); ++comp_idx) {
        // Compare to other segments on this edge and other edges
        unsigned int start = 0; // On different edge, look at all segments
        if (comp_idx == edge_idx) {
          start = i+1; // On same edge, look only at "future" segments (+1).
        }
        for (unsigned int j=start; j < polygon[comp_idx].size()-1; ++j) {
          p2 = polygon[comp_idx][j];
          q2 = polygon[comp_idx][j+1];

          int v = commonVertices(p1, q1, p2, q2);
          if (v == 0) { // If 1 vertex is shared, do not count as intersection
            if (doIntersect(p1, q1, p2, q2)) {
              total++;
            }
          } else if (v == 2) {
            // If 2 vertices are shared, this is the same segment and is an error
            throw customException("cannot check intersection, segments are identical");
          }
        }
      }
    }
  }

  return total;
}

bool Colloid::isMotifInside(const int N = 10) {
  /**
   * Test if the motif is inside the tile boundary.
   *
   * The tile is discretized into points and treated as a polygon.
   * All points on the motif are then tested if they are inside that
   * polygon; if they all are, then the motif is considered to be
   * inside the tile.
   *
   * @param N Number of points to discretize each edge into.
   *
   * @returns Boolean if the motif is completely inside the tile.
   *
   * @throws customException if tile or motif has not been assigned yet.
   */

  if (!motif_assigned_ || !tile_assigned_) {
    throw customException("must assign assign tile and motif first");
  }

  vector<vector<dvec2>> polygon = buildTilePolygon(N);

  // Check each point in motif
  const vector<vector<double>> c = m_.getCoords();
  for (size_t i = 0; i < c.size(); ++i) {
    dvec2 point = {c[i][0], c[i][1]};
    if (!isInside(polygon, point)) {
      return false;
    }
  }

  return true;
}

const vector<vector<dvec2>> Colloid::buildTilePolygon(const int N) {
  /**
   * Build a set of edges that create a polygonal approximation of the tile.
   * This places a total of N points evenly (in Bezier parameterization) along each edge.
   */

  // 0 -- p_1 -- p_2 -- ... p_(n-2) -- 1
  // Place N points between ends of the edge - space evenly
  assert(N > 2); // N refers to total and there are always 2 endpoints
  const double du = (1.0 - 0.0) / (N-1); // Bezier ends are at 0 and 1, place these
  const vector<double> u0(getTile().numEdgeShapes(), du); // Start placing from du

  vector<int> boundary_ids;
  vector<vector<double>> boundary_coords;
  vector<vector<double>> tile_control_points;
  vector<vector<dvec2>> polygon;

  // Place N-2 points between the ends
  perimeter_(u0, getDform(), du, N-2, getTileScale(), &boundary_ids, &boundary_coords,
             &tile_control_points, &polygon);

  return polygon;
}

double Colloid::fractionMotifInside(const int N = 20) {
  /**
   * Compute what fraction of the motif's points are inside the tile boundary.
   *
   * The tile is discretized into points and treated as a polygon.
   *
   * @param N Number of points to discretize each edge into.
   *
   * @returns double Number of points inside the tile.
   *
   * @throws customException if tile or motif has not been assigned yet.
   */

  if (!motif_assigned_ || !tile_assigned_) {
    throw customException("must assign assign tile and motif first");
  }

  double outside = 0.0;
  vector<vector<dvec2>> polygon = buildTilePolygon(N);

  // Check each point in motif
  const vector<vector<double>> c = m_.getCoords();
  for (size_t i = 0; i < c.size(); ++i) {
    dvec2 point = {c[i][0], c[i][1]};
    if (!isInside(polygon, point)) {
      outside += 1.0;
    }
  }

  return 1.0 - outside / c.size();
}

bool Colloid::isTileFundamental() {
  /**
   * Check if the tile is one of the 46 fundamental domain tiles or not.
   *
   * @returns Boolean of if the tile is a fundamental domain.
   *
   * @throws customException if tile has not been assigned yet.
   */
  if (!tile_assigned_) {
    throw(customException("tile has not been assigned yet"));
  }

  bool fundamental = false;
  const int tt = static_cast<int>(tile_.getTilingType());
  for (unsigned int i = 0; i < sizeof(FD_TYPES) / sizeof(FD_TYPES[0]); ++i) {
    if (tt == FD_TYPES[i]) {
      fundamental = true;
    }
  }

  return fundamental;
}

double Colloid::tileArea() {
  /**
   * Compute the area of the tile (colloid).
   *
   * Since all edge deformations should "cancel out" by symmetry, we use the
   * control points of the edges to define the polygonal patch that defines
   * the area of the colloid. In particular, this code is based on the
   * discussion at https://www.wikihow.com/Calculate-the-Area-of-a-Polygon.
   *
   * @returns Area of the tile.
   *
   * @throws customException if tile has not been assigned yet.
   */

  if (!tile_assigned_) {
    throw(customException("tile has not been assigned yet"));
  }

  const int n = tile_control_points_.size();
  double sum1 = 0.0, sum2 = 0.0;
  for (int i = 0; i < n; ++i) {
    if (i == n - 1) {
      sum1 += tile_control_points_[i][0] * tile_control_points_[0][1];
      sum2 += tile_control_points_[i][1] * tile_control_points_[0][0];
    } else {
      sum1 += tile_control_points_[i][0] * tile_control_points_[i + 1][1];
      sum2 += tile_control_points_[i][1] * tile_control_points_[i + 1][0];
    }
  }
  return (sum1 - sum2) / 2.;
}

void Colloid::buildBoundary_() {
  /**
   * Build the colloid's boundary points.
   *
   * The boundary is made using the "4+1" rule. so there are 4 interacting
   * points and 1 "stop codon" on each edge.
   */

  if (!tile_assigned_) { // Makes sure the tile has been assigned
    throw customException(
        "must assign tile before building boundary");
  }

  // Place one extra stop codon for now and choose which to use later
  vector<vector<dvec2>> polygon;
  perimeter_(getU0(), getDform(), getDU(), 1 + 4 + 1, getTileScale(), &boundary_ids_,
             &boundary_coords_, &tile_control_points_, &polygon);
}

void Colloid::perimeter_(vector<double> u0, vector<double> df, double du, int n, double scale,
                         vector<int>* boundary_ids,
                         vector<vector<double>>* boundary_coords,
                         vector<vector<double>>* tile_control_points,
                         vector<vector<dvec2>>* polygon) {
  /**
   * Compute points along the tile's perimeter.
   *
   * A Bezier curve defines each edge on the tile, following the Tactile
   * library. Each edge assumed to be "impacted" by a sphere to create a
   * curvature of a single sign along its edge. The identities of each
   * symmetrically unique point along the boundary are identified in ascending
   * order, without gap. 0 is used for stop codons, positive integers for
   * interacting points, -1 is used for control points at the ends of Bezier 
   * curves.
   *
   * @param[in] u0 Starting point for boundary points along Bezier curve. Should
   * be in (0,1) for each "unique edge" as defined by Tactile.
   * @param[in] df Deformation magnitude (positive or negative) along Bezier curve. 
   * Should be provided for each "unique edge" as defined by Tactile.
   * @param[in] du Gap along Bezier curve between boundary points. Should be
   * < (1-u0)/n.
   * @param[in] n Total number of points to place along each edge (incl. stop codons).
   * @param[in] scale The default Tactile tile is isotropically scaled by this
   * factor.
   * @param[out] boundary_ids Chemical identities (integers > 0) of boundary
   * points. These are listed in ascending order from 1 up, without gap.
   * @param[out] boundary_coords Coordinates of points on tile's boundary.
   * @param[out] tile_control_points Terminal control points from the ends of Bezier  
   * curves defining the edges. 
   * @param[out] polygon For each edge ("part") consecutive points (including ends) 
   * are given. This includes the points placed + 2 (one for each endpoint / terminal 
   * control point). Edges do not necessarily have the same orientation.
   */

  // Put points on perimeter (this returns stop codons AND control points).
  // Performs checks that u0 and df consistent with tile.
  vector<vector<dvec2>> edges = perimeter_edges_(u0, df, du, n, scale);

  // Iterate over the edges of a single tile, asking the tiling to
  // tell you about the geometric information needed to transform
  // the edge shapes into position.  Note that this iteration is over
  // tiling "parts" which goes over entire edges for I and J, and each "half"
  // for S and U.
  vector<int> identity, cp_idx;
  vector<dvec2> shape;
  polygon->clear();
  for (auto i : tile_.parts()) {
    polygon->resize(polygon->size()+1);
    vector<int> pattern; // e.g., [-1, 0, {1, 2, 3, ...}, 0, -1] where -1 = CP, 0 = SC

    // Get the relevant edge shape created above using i->getId().
    const vector<dvec2> ed = edges[i->getId()];

    // Also get the transform that maps to the line joining consecutive
    // tiling vertices.
    const glm::dmat3& T = i->getTransform();

    // Interacting points on U and S edges are centro-symmetrically labelled
    if (i->getShape() == 1 || i->getShape() == 2) {
      // Only put a stop codon at the end and made half an edge
      const int half_size = (n-2)/2+1; // Round down if needed
      pattern.resize(half_size, 0);

      // Stop codons are type 0
      pattern[half_size-1] = 0;

      for (int j=1; j < half_size; ++j) {
        pattern[j-1] = i->getId() * (n-2) + j; // n-2 non-STOP codons
      }
    } else if (is_undirected_mirror(tile_.getTilingType(), i->getId())) {
      // Undirected I edge
      assert(i->getShape() == 3);

      const int half_size = (n-2)/2+1; // Round down if needed
      pattern.resize(2*half_size, 0);

      // Stop codons are type 0
      pattern[0] = 0;
      pattern[pattern.size()-1] = 0;

      // Symmetric labelling around midpoint
      for (int j=1; j < half_size; ++j) {
        pattern[j] = i->getId() * (n-2) + j; // n-2 non-STOP codons
        pattern[pattern.size()-1-j] = pattern[j];
      }
    } else {
      // Directed I edge or J edge
      pattern.resize(n, 0);

      // Stop codons are type 0
      pattern[0] = 0;
      pattern[n-1] = 0;

      for (int j=1; j < n-1; ++j) {
        pattern[j] = i->getId() * (n-2) + j; // n-2 non-STOP codons
      }
    }

    // If i->isReversed() is true, we need to run the parameterization
    // of the path backwards.
    if (i->isReversed()) {
      identity.push_back(-1);  // -1 for control points
      cp_idx.push_back(identity.size()-1);
      for (size_t idx = 0; idx < ed.size(); ++idx) { // 1
        shape.push_back(T * dvec3(ed[ed.size() - 1 - idx], 1.0 * scale));
        polygon->back().push_back(shape.back());
      }
      for (size_t idx = 0; idx < pattern.size(); ++idx) {
        identity.push_back(pattern[pattern.size() - 1 - idx]);
      }
      identity.push_back(-1);  // -1 for control points
      cp_idx.push_back(identity.size()-1);
    } else {
      identity.push_back(-1);  // -1 for control points
      cp_idx.push_back(identity.size()-1);
      for (size_t idx = 0; idx < ed.size(); ++idx) { // 1
        shape.push_back(T * dvec3(ed[idx], 1.0 * scale));
        polygon->back().push_back(shape.back());
      }
      for (size_t idx = 0; idx < pattern.size(); ++idx) {
        identity.push_back(pattern[idx]);
      }
      identity.push_back(-1);  // -1 for control points
      cp_idx.push_back(identity.size()-1);
    }
  }

  // Tactile does not seem to provide an orientation-consistent method of
  // figuring out the ends of the Bezier curves naturally; the edges are
  // traversed in this order, though, so we can use this to infer the
  // order of the points. If a better way is found this can be replaced.
  // This order is needed to compute the tile area correctly.  Other things
  // like isMotifInside() to not need this ordering.

  // Get control points at the "ends" of edges in counterclockwise order
  int iter = 0;
  vector<vector<int>> cp_edges;
  for (auto i : tile_.parts()) {  // Iterator goes in counterclockwise edge order from Tactile
    vector<int> edges;
    if (i->getShape() == 1 || i->getShape() == 2) { // S or U edge
      if (i->isSecondPart()) {
        vector<int> idx = {cp_idx[iter-2], cp_idx[iter-1], cp_idx[iter], cp_idx[iter+1]};
        edges = unique_(shape, idx); // "Delete" the overlapping set of 2 to get edges
        assert(edges.size() == 2);
        cp_edges.push_back(edges); // Save the outer 2 points
      }
    } else { // J or I edge
      edges = {cp_idx[iter], cp_idx[iter+1]};
      cp_edges.push_back(edges); // Save the outer 2 points
    }
    iter += 2;
  }

  // Traverse in counterclockwise order
  dvec2 starting_point(0.0, 0.0);
  tile_control_points->clear();
  for (unsigned int i=0; i < cp_edges.size(); ++i) {
    if (i > 0) { // From the second edge onward...
      vector<int> last_edge = cp_edges[i-1], curr_edge = cp_edges[i];

      unsigned int ov_ = 0, new_ = 0;
      new_idx_(shape, last_edge, curr_edge, ov_, new_);
      tile_control_points->push_back({shape[curr_edge[new_]].x, shape[curr_edge[new_]].y}); // next CP in CC order
      if (i == 1) {
        starting_point = shape[curr_edge[ov_]]; // One that is duplicated one this iteration is the starting point you end up returning to
      }
    }
  }
  tile_control_points->push_back({starting_point.x, starting_point.y});

  // Count any gaps in the identities from S and U edges so identities are contiguous
  map<int, int> adjust;
  int ctr = 0;
  for (int i = 0; i <= *max_element(identity.begin(), identity.end()); ++i) {
    if (find(identity.begin(), identity.end(), i) == identity.end()) {
      ++ctr;
    } else {
      adjust[i] = ctr;
    }
  }

  boundary_ids->clear();
  boundary_coords->clear();

  int skip = 0;
  for (size_t j = 0; j < shape.size(); ++j) {
    vector<double> c = {shape[j].x, shape[j].y};
    if (identity[j] >= 0) {  // Ignore control points
      if ((identity[j] == 0) &&
          (skip % 2 ==
           0)) {  // Skip every other stop codon in fixed orientation
        boundary_coords->push_back(c);
        boundary_ids->push_back(identity[j] - adjust[identity[j]]);
        ++skip;
      } else if (identity[j] == 0) {
        ++skip; // Skip every other stop codon in fixed orientation
      } else {
        boundary_coords->push_back(c);
        boundary_ids->push_back(identity[j] - adjust[identity[j]]);
      }
    }
  }
}

vector<int> unique_(const vector<dvec2>& shape, const vector<int>& idx, const double eps) {
  /**
   * Find indices of points which are unique.
   */
  vector<int> u;
  for (unsigned int i=0; i < idx.size(); ++i) {
    u.push_back(idx[i]);
    for (unsigned int j=0; j < idx.size(); ++j) {
      if (i != j) {
        double d2 = pow(shape[idx[i]].x - shape[idx[j]].x, 2) + pow(shape[idx[i]].y - shape[idx[j]].y, 2);
        if (d2 < eps*eps) {
          u.pop_back();
          break;
        }
      }
    }
  }

  return u;
}

void new_idx_(const vector<dvec2>& shape, vector<int>& last_edge, vector<int>& curr_edge, unsigned int& ov_, unsigned int& new_, const double eps) {
  /**
   * Compare end points of 2 consecutive edges (in counterclockwise order) to determine the
   * "leading" point, and the point the 2 edges share.
   */
  assert(last_edge.size() == 2);
  assert(curr_edge.size() == 2);

  for (unsigned int i=0; i < 2; ++i) {
    for (unsigned int j=0; j < 2; ++j) {
      const double d2 = pow(shape[curr_edge[i]].x - shape[last_edge[j]].x, 2) + pow(shape[curr_edge[i]].y - shape[last_edge[j]].y, 2);
      if (d2 < eps) {
        ov_ = i; // i overlaps one from last_edge
        new_ = (i+1)%2; // The other one is the new one
        return;
      }
    }
  }
  throw customException("unable to find overlapping point");
}

void Colloid::initMotif_(double max_scale_factor = 10.0,
                         double min_scale_factor = 0.1, int n_scale_incr = 1000,
                         int N = 20, bool debug = false) {
  /**
   * Initialize the motif by isotropically expanding the tile to "just" enclose it.
   *
   * Motif orientation is unaffected for FD tiles, however, if there are constraints
   * from non-FD tiles the motif may be changed to match this.
   *
   * This attempts to place the motif entirely inside the tile's boundary after
   * it has been initialized. The default tile settings from the Tactile library
   * are generally convenient and convex so it is recommended that these be kept
   * until this initialization is run. If a motif does not fit, then the tile is
   * expanded until it "just fits"; similarly, if it initially fits, the tile is
   * shrunk until it cannot be shrunk any more.
   *
   * Note that overlapping tile edges are not accounted for here - only scaling
   * is performed so it can be wise to start from df = 0 (straight edges) with
   * the Tactile default shape parameters.
   *
   * @param max_scale_factor Maximum factor to scale the default Tactile tile
   * size to trying to fit the motif inside of it.
   * @param min_scale_factor Minimum factor to scale the default Tactile tile
   * size to trying to fit the motif inside of it.
   * @param n_scale_incr Number of intermediate scale factors to try.
   * @param N Number of points to place along edges when discretizing into a
   * polygon.
   * @param debug Whether or not to print out XYZ file if the initialization
   * fails.  This always writes to _debug_.xyz.
   *
   * @throws customException if initialization fails for any reason.
   */

  if (!tile_assigned_ || !motif_assigned_) { // Makes sure the tile & motif have been assigned
    throw customException(
        "must assign tile and motif before initializing the motif");
  }

  if (isTileFundamental()) {
    // If no restrictions, place motif COM on tile COM and shrink/expand to fit.
    // Tactile default parameters provide convenient, convex, starting points
    // at least one control point is at (0,0).  This way, we can scale points
    // easily by multiplying by a number.

    // 1. Place motif so its COM is on tile COM
    // tile COM is estimated from the functional points on the perimeter
    vector<double> motif_com = m_.getCOM(), tile_com = boundaryCOM(), dx;
    for (size_t i = 0; i < motif_com.size(); ++i) {
      dx.push_back(tile_com[i] - motif_com[i]);
    }
    m_.translate(dx);

    // 2. Expand/contract the tile until motif "just" fits
    double min_scale = getTileScale() * min_scale_factor;
    double max_scale = getTileScale() * max_scale_factor;
    double orig_scale = getTileScale();
    bool found = false;
    if (isMotifInside(N)) {
      // Shrink the tile to fit
      double last_scale = getTileScale();
      for (int i = 0; i <= n_scale_incr; ++i) {
        double new_scale =
            (orig_scale - (orig_scale - min_scale) / n_scale_incr * i);
        setTileScale(new_scale);

        vector<double> params = m_.getParameters();
        m_.setParameters({new_scale*tile_com[0], new_scale*tile_com[1], params[2]});

        if (!isMotifInside(N)) {
          // Move motif back to last position
          m_.setParameters({last_scale*tile_com[0], last_scale*tile_com[1], params[2]});
          setTileScale(last_scale);
          found = true;
          break;
        }
        last_scale = getTileScale();
      }
      if (!found) {
        if (debug) {
          buildBoundary_();
          dumpXYZ("_debug_.xyz", true);
        }
        throw(customException("unable to shrink tile around motif"));
      }
    } else {
      // Expand the tile to fit
      for (int i = 0; i <= n_scale_incr; ++i) {
        double new_scale = (orig_scale + (max_scale - orig_scale) / n_scale_incr * i);
        setTileScale(new_scale);

        vector<double> params = m_.getParameters();
        m_.setParameters({new_scale*tile_com[0], new_scale*tile_com[1], params[2]});

        if (isMotifInside(N)) {
          // Motif and tile now in acceptable positions
          found = true;
          break;
        }
      }
      if (!found) {
        if (debug) {
          buildBoundary_();
          dumpXYZ("_debug_.xyz", true);
        }
        throw(customException("unable to expand tile around motif"));
      }
    }
  } else {
    // We use setParameters() because internally revisions occur to keep parameters
    // within constraints determined by symmetry.

    // 1. Place motif so its COM is on tile COM
    vector<double> orig_colloid_params = getParameters();
    vector<double> curr_params = orig_colloid_params, last_params;

    const vector<double> tile_com = scale_coords_(boundaryCOM());
    curr_params[0] = tile_com[0];
    curr_params[1] = tile_com[1];
    setParameters(curr_params); // setParameters() works in "scaled coordinates"
    last_params = curr_params;

    // 2. Expand/contract the tile until motif "just" fits
    double min_scale = getTileScale() * min_scale_factor;
    double max_scale = getTileScale() * max_scale_factor;
    double orig_scale = getTileScale();
    bool found = false;
    if (isMotifInside(N)) {
      // Shrink the tile to fit
      for (int i = 0; i <= n_scale_incr; ++i) {
        double new_scale =
            (orig_scale - (orig_scale - min_scale) / n_scale_incr * i);
        curr_params[curr_params.size()-1] = new_scale;
        setParameters(curr_params); // No need to changed scaled coords of Motif!

        if (!isMotifInside(N)) {
          // Move motif back to last position
          setParameters(last_params);
          found = true;
          break;
        }
        last_params = curr_params;
      }
      if (!found) {
        if (debug) {
          buildBoundary_();
          dumpXYZ("_debug_.xyz", true);
        }
        throw(customException("unable to shrink tile around motif"));
      }
    } else {
      // Expand the tile to fit
      for (int i = 0; i <= n_scale_incr; ++i) {
        double new_scale = (orig_scale + (max_scale - orig_scale) / n_scale_incr * i);
        curr_params[curr_params.size()-1] = new_scale;
        setParameters(curr_params); // No need to changed scaled coords of Motif!

        if (isMotifInside(N)) {
          // Motif and tile now in acceptable positions
          found = true;
          break;
        }
      }
      if (!found) {
        if (debug) {
          buildBoundary_();
          dumpXYZ("_debug_.xyz", true);
        }
        throw(customException("unable to expand tile around motif"));
      }
    }
  }

  buildBoundary_();  // Re-build based on final tile_scale_
  return;
}

vector<vector<dvec2>> Colloid::perimeter_edges_(vector<double> u0, vector<double> df, double du, int n,
                                                double scale) {
  /**
   * Create the Bezier curves that will serve as tile edges.
   *
   * These curves canonically go from (0,0) to (1,0) with some curvature between
   * them, and are uniformly scaled.  They will be rotated into place later on
   * by the Tactile library.
   *
   * @param u0 Starting point for boundary points along each Bezier curve. Should be
   * in (0,1). S and U edges start from their center and essentially ignore this.
   * Undirected I edges do the same and also ignore this parameter.
   * @param df Deformation (positive or negative) along each Bezier curve. 
   * @param du Gap along Bezier curve between boundary points. Should be < 1.
   * @param n Total number of points to place along each a completed edge. Note 
   * that U and S edges are cut in "half" and n refers to the total edge, so we
   * generate only one half of these points here.
   * @param scale The default Tactile tile is isotropically scaled by this
   * factor.
   *
   * @returns Discretized points along each edge, including Bezier control
   * points at the ends. For S and U: [CP_left, c1, c2, ..., c_((n-2)/2+1)=STOP, 
   * CP_right], while for I and J: [CP_left, c1=STOP, c2, ..., c_n=STOP, 
   * CP_right]. Points move from CP_left to CP_right.
   */

  if (!tile_assigned_) { // Makes sure the tile has been assigned
    throw customException(
        "must assign tile before building perimeter");
  }

  // Create a vector to hold some edge shapes.  The tiling tells you
  // how many distinct edge shapes you need, but doesn't know anything
  // about how those shapes might be represented.  It simply assumes
  // that each one will be a curve from (0,0) to (1,0).  The tiling
  // provides tools to let you map those curves into position around
  // the outline of a tile.  All the curves below have exactly four
  // control points.
  vector<vector<dvec2>> edges;

  if (u0.size() != getTile().numEdgeShapes()) {
    throw customException("incorrect number of u0 parameters provided");
  }
  if (df.size() != getTile().numEdgeShapes()) {
    throw customException("incorrect number of df parameters provided");
  }
  assert(n > 0);

  double u0_ = 0, du_ = 0;
  int n_ = 0;

  // Generate edge shapes.
  for (U8 idx = 0; idx < tile_.numEdgeShapes(); ++idx) {
    vector<dvec2> ej;

    // Canonical Tactile coordinates
    dvec2 cp_left(0.0, 0.0), cp_right(1.0 * scale, 0.0), dummy;

    double sphere_deform = df[idx];

    // Define Bezier Curve that is sort of like a sphere impacting
    // This could definitely be changed, but will introduce more
    // free parameters.
    ej.push_back(dvec2(0.0, 0.0));
    ej.push_back(dvec2(1 / 3. * scale, sphere_deform * scale));
    ej.push_back(dvec2(2 / 3. * scale, sphere_deform * scale));
    ej.push_back(dvec2(1.0 * scale, 0.0));

    // Now, depending on the edge shape class, enforce symmetry
    // constraints on edges.
    const bool undirected_edge = is_undirected_mirror(tile_.getTilingType(), idx);
    switch (tile_.getEdgeShape(idx)) {
      case J:
        u0_ = u0[idx];
        du_ = du;
        n_ = n;
        break;
      case U: // double dU since edge is half the size
        du_ = -2.0*du; // Walk "away" from midpoint
        u0_ = 1.0+du_/2.0; // Start U from x=1 = midpoint and walk away
        assert(n >= 2);
        n_ = (n-2)/2+1; // Only put 1 stop codon at the end

        // Reverse CP order - needed for line segments / self-intersection tests
        dummy = cp_left;
        cp_left = cp_right;
        cp_right = dummy;
        break;
      case S: // double dU since edge is half the size
        du_ = -2.0*du; // Walk "away" from midpoint
        u0_ = 1.0+du_/2.0; // Start S from x=1 = midpoint and walk away
        assert(n >= 2);
        n_ = (n-2)/2+1; // Only put 1 stop codon at the end

        // Reverse CP order - needed for line segments / self-intersection tests
        dummy = cp_left;
        cp_left = cp_right;
        cp_right = dummy;
        break;
      case I:
        if (undirected_edge) {
          // Some edges on non-FD tiles are symmetric about their midpoint
          // and should be treated like "straight U" edges.
          ej[1].y = 0.0;
          ej[2].y = 0.0;
          du_ = du;
          assert(n >= 2);
          const int n_side = (n-2)/2+1; // On one half
          u0_ = 0.5-(n_side-0.5)*du_; // Start from left hand side of midpoint x=0.5
          n_ = 2*n_side;
        } else { // "Normal" mirror edge
          ej[1].y = 0.0;
          ej[2].y = 0.0;
          u0_ = u0[idx];
          du_ = du;
          n_ = n;
        }
        break;
    }

    // From Bezier curve determine location of points along the
    // deformed edge. Here, i'll assume the stop codons are in contact
    // with the other ones. Later we will decide which stop codon to drop.
    // For S and U (halves), only put ~n/2 points and start from 1-du/2 instead of 0+u0.
    // These edges have only one stop codon placed at the end.
    vector<dvec2> coords;
    coords.push_back(cp_left);
    for (int k = 0; k < n_; ++k) {
      coords.push_back(bezier(ej[0], ej[1], ej[2], ej[3], u0_ + k * du_));
    }
    coords.push_back(cp_right);

    edges.push_back(coords);
  }

  return edges;
}

bool is_undirected_mirror(const int ih_number, const int edge_idx) {
  /**
   * Manual logic that indicates if a tile's edge is an undirected
   * I edge.  In this case, we should treat it like a 
   * "straight U" edge.  There is nothing fundamental about this; 
   * this is just how Tactile is programmed.
   */
  if (ih_number == 17) {
    if (edge_idx == 0) {
      return true;
    }
  } else if (ih_number == 20) {
    if (edge_idx == 0) {
      return true;
    }
  } else if (ih_number == 26) {
    if (edge_idx == 2) {
      return true;
    }
  } else if (ih_number == 29) {
    if (edge_idx == 0) {
      return true;
    }
  } else if (ih_number == 40) {
    if (edge_idx == 1) {
      return true;
    }
  } else if (ih_number == 67) {
    if (edge_idx == 1 || edge_idx == 2) {
      return true;
    }
  } else if (ih_number == 72) {
    if (edge_idx == 0 || edge_idx == 1) {
      return true;
    }
  } else if (ih_number == 76) {
    if (edge_idx == 0) {
      return true;
    }
  } else if (ih_number == 82) {
    if (edge_idx == 1) {
      return true;
    }
  } else if (ih_number == 91) {
    if (edge_idx == 1) {
      return true;
    }
  } else if (ih_number == 93) {
    if (edge_idx == 0) {
      return true;
    }
  }

  return false;
}

void Colloid::load(const string filename) {
  /**
   * Load a colloid from a JSON file.
   *
   * @param filename Name of file to read from.
   *
   * @throws customException if anything goes wrong.
   */

  ifstream in(filename);
  json j;
  in >> j;

  if (j.contains("Motif")) {
    try {
      Motif m;
      m.setCoords(j["Motif"]["coords"].get<vector<vector<double>>>(),
                  j["Motif"]["parameters"][2].get<double>());
      m.setTypes(j["Motif"]["types"].get<vector<string>>());
      m.setParameters(
          j["Motif"]["parameters"]
              .get<vector<double>>());  // Unnecessary, but for good measure
      m.setSymmetry(j["Motif"]["symmetry"].get<string>());
      setMotif(m);
    } catch (const exception& e) {
      throw(customException("unable to load Motif"));
    }
  }

  if (j.contains("Tile")) {
    try {
      IsohedralTiling t(j["Tile"]["ih_type"].get<int>());
      vector<double> p = j["Tile"]["parameters"].get<vector<double>>();
      double params[t.numParameters()];
      for (int i = 0; i < t.numParameters(); ++i) {
        params[i] = p[i];
      }
      t.setParameters(params);
      setTile(t);
    } catch (const exception& e) {
      throw(customException("unable to load Tile"));
    }
  }

  if (j.contains("Properties")) {
    try {
      setDU(j["Properties"]["edge_du"].get<double>());
      setU0(j["Properties"]["edge_u0"].get<vector<double>>());
      setTileScale(j["Properties"]["tile_scale"].get<double>());
      setDform(j["Properties"]["edge_df"].get<vector<double>>());
      boundary_ids_ = j["Properties"]["boundary_ids"].get<vector<int>>();
      boundary_coords_ =
          j["Properties"]["boundary_coords"].get<vector<vector<double>>>();
      tile_control_points_ =
          j["Properties"]["tile_control_points"].get<vector<vector<double>>>();
      params_ = j["Properties"]["params"].get<vector<double>>();
      built_ = j["Properties"]["built"].get<bool>();
    } catch (const exception& e) {
      throw(customException("unable to load Properties"));
    }
  }
}

void Colloid::dump(const string filename) {
  /**
   * Dump the colloid to a JSON file.
   *
   * @param filename Name of file to write to. Will overwrite by default.
   *
   * @throws customException if anything goes wrong.
   */

  json j;

  try {
    if (motif_assigned_) {
      j["Motif"] = {{"coords", m_.getCoords()},
                    {"types", m_.getTypes()},
                    {"parameters", m_.getParameters()},
                    {"symmetry", m_.getSymmetry()}};
    }

    if (tile_assigned_) {
      double params[tile_.numParameters()];
      tile_.getParameters(params);
      vector<double> p;
      for (int i = 0; i < tile_.numParameters(); ++i) {
        p.push_back(params[i]);
      }
      j["Tile"] = {{"ih_type", static_cast<int>(tile_.getTilingType())},
                   {"parameters", p}};
    }

    if (built_) {
      j["Properties"] = {{"boundary_ids", boundary_ids_},
                         {"boundary_coords", boundary_coords_},
                         {"edge_df", getDform()},
                         {"edge_du", getDU()},
                         {"edge_u0", getU0()},
                         {"tile_scale", getTileScale()},
                         {"tile_control_points", tile_control_points_},
                         {"params", getParameters()},
                         {"built", built_}};
    }

    std::ofstream out(filename);
    out << std::setw(4) << j << std::endl;
  } catch (const exception& e) {
    throw(customException("unable to dump colloid"));
  }
}

void Colloid::dumpXYZ(const string filename, const bool full = false) {
  /**
   * Dump the colloid to an XYZ file. The `full` option can improve
   * visualization. This lifts the motif in the z-direction "above" the tile,
   * and also prints the control points (type "CP") to the file.
   *
   * @param filename Name of file to write to. Will overwrite by default.
   * @param full If true, print motif "above" the tile and show control points.
   *
   * @throws customException if anything goes wrong.
   */

  vector<vector<double>> motif_coords = m_.getCoords();
  vector<string> motif_types = m_.getTypes();

  if (motif_coords.empty()) {
    throw(customException("motif coordinates have not been assigned"));
  }
  if (motif_types.empty()) {
    throw(customException("motif identities have not been assigned"));
  }
  if (boundary_coords_.empty()) {
    throw(customException("tile boundary coordinates have not been assigned"));
  }
  if (boundary_ids_.empty()) {
    throw(customException("tile boundary identities have not been assigned"));
  }

  try {
    ofstream xyz(filename);

    int add = 0, shift = 0;
    if (full) {
      add = tile_control_points_.size();
      shift = 1;
    }

    xyz << boundary_coords_.size() + motif_coords.size() + add << endl;
    xyz << endl;

    for (size_t i = 0; i < boundary_coords_.size(); ++i) {
      xyz << boundary_ids_[i] << "\t" << boundary_coords_[i][0] << "\t"
          << boundary_coords_[i][1] << "\t" << 0 << endl;
    }

    for (size_t i = 0; i < motif_coords.size(); ++i) {
      xyz << motif_types[i] << "\t" << motif_coords[i][0] << "\t"
          << motif_coords[i][1] << "\t" << shift << endl;
    }
    if (full) {
      for (size_t i = 0; i < tile_control_points_.size(); ++i) {
        xyz << "CP"
            << "\t" << tile_control_points_[i][0] << "\t"
            << tile_control_points_[i][1] << "\t" << 0 << endl;
      }
    }
  } catch (const exception& e) {
    throw(customException("unable to write to xyz file"));
  }

  return;
}

void Colloid::unitCell(vector<vector<double>>* coords, vector<string>* types,
                       vector<vector<double>>* box, const int nx,
                       const int ny, const double tol) {
  /**
   * Make a unit cell (nx by ny) out of the motif. All particles are wrapped
   * into the box. This is useful when checking the final symmetry of the
   * system.
   *
   * @param[in] nx Number of copies to make in the "t1" direction.
   * @param[in] ny Number of copies to make in the "t2" direction.
   * @param[out] coords (x,y) coordinates of particles in unit cell.
   * @param[out] types Chemical type of each particle.
   * @param[out] box Box ( (v1.x, v1.y), (v2.x, v2.y) ) translation vectors.
   */

  if (!motif_assigned_ || !tile_assigned_) {
    throw customException("must assign assign tile and motif first");
  }
  types->clear();
  coords->clear();
  box->clear();

  // Get unit cell by combining the relevant "aspects" (Tactile terminology)
  const int n = tile_.numAspects();
  vector<vector<double>> mc = m_.getCoords();
  vector<string> mt = m_.getTypes();
  glm::dvec2 t1 = tile_.getT1() * getTileScale(), t2 = tile_.getT2() * getTileScale();
  vector<dvec2> glm_coords;
  for (int i = 0; i < n; ++i) {
    dmat3 T = tile_.getAspectTransform(i);
    vector<string> t = m_.getTypes();
    for (size_t j = 0; j < mc.size(); ++j) {
      dvec3 c = T * dvec3(mc[j][0], mc[j][1], getTileScale());
      glm_coords.push_back(c);
    }
  }

  box->push_back({t1.x, t1.y});
  box->push_back({t2.x, t2.y});

  // Wrap and make copies by translating
  dmat2 H = {t1, t2};
  dmat2 H_T = transpose(H);
  dmat2 H_T_inv = inverse(H_T);
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (size_t i = 0; i < glm_coords.size(); ++i) {
        dvec2 a = glm_coords[i] * H_T_inv;
        dvec2 x = dvec2(mod(a.x, 1.0), mod(a.y, 1.0)) * H_T;
        vector<double> y = {x.x + ix * t1.x + iy * t2.x,
                            x.y + ix * t1.y + iy * t2.y};

        // Do not print duplicates
        bool duplicate = false;
        for (unsigned int k = 0; k < coords->size(); ++k) {
          const double d2 = pow(coords->at(k)[0] - y[0], 2) + pow(coords->at(k)[1] - y[1], 2);
          if (d2 < tol*tol) {
            duplicate = true;
            break;
          }
        }
        if (!duplicate) {
          coords->push_back(y);
          types->push_back(mt[i % mc.size()]);
        }
      }
    }
  }
}
