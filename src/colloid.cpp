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

bool pip(const vector<vector<double>>& polygon, const vector<double>& point) {
  /**
   * Check if a point is inside a polygon.
   *
   * In particular, this code is based on the discussion found at
   * https://stackoverflow.com/questions/11716268/point-in-polygon-algorithm.
   *
   * @param polygon Matrix of (x,y) coordinates of polygon's vertices.
   * @param point Coordinate to check.
   *
   * @returns Boolean indicating if `point` is inside `polygon`.
   */

  bool c = false;
  for (size_t i = 0, j = polygon.size() - 1; i < polygon.size(); j = i++) {
    if (((polygon[i][1] >= point[1]) != (polygon[j][1] >= point[1])) &&
        (point[0] <= (polygon[j][0] - polygon[i][0]) *
                             (point[1] - polygon[i][1]) /
                             (polygon[j][1] - polygon[i][1]) +
                         polygon[i][0])) {
      c = !c;
    }
  }

  return c;
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

Colloid::Colloid(Motif m, IsohedralTiling t, double tile_u0, bool debug)
    : tile_(IsohedralTiling(1)) {
  /**
   * Instantiate a new colloid if all parameters are known (preferred method).
   *
   * The Tactile library requires an isohedral tile number be given at
   * instantiation so the colloid is initialized, by default, to have tile IH1.
   * This can be changed later.
   */

  defaults_();
  setU0(tile_u0);

  setMotif(m);
  setTile(t);
  init(debug);  // This will find a good value for tile_scale_
}

void Colloid::init(bool debug) {
  buildBoundary_();
  initMotif_(10.0, 0.1, 1000, 20, debug);
  built_ = true;
}  // Initialize the colloid.

void Colloid::defaults_() {
  /**
   * Assign default paramter values.
   */

  tile_assigned_ = false;
  motif_assigned_ = false;
  built_ = false;
  setTileScale(1.0);
  setU0(0.0);
  sphere_deform_ = 0.25;
  setDU(0.1);
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

void Colloid::setParameters(const vector<double>& params) {
  /**
   * Assign all parameters defining the colloid.
   *
   * The parameters are an unrolled (double) vector useful for optimization
   * schemes. They are: [{motif scaled_com_x, motif scaled_com_y}, motif_theta,
   * {v0, v1, etc. for tile}, edge_u0, tile_scale]. The scaled com's are values
   * [0, 1] which describe the location in terms of the bounding box around the
   * tile. Motif COM coordinates are only used for tiles which are fundamental
   * domains.
   *
   * The motif will track its coordinates in absolute units, but the colloid
   * uses reduced units since it knows about the tile and means that these
   * variables can be given reasonable bounds to an optimizer in advance.
   *
   * @param params Parameter vector described above.
   *
   * @throws customException if tile or motif has not been assigned yet.
   */

  if (!motif_assigned_ ||
      !tile_assigned_) {  // Makes sure the tile and colloid have been assigned
    throw customException(
        "must assign tile and motif before setting new parameters");
  }

  if (params.size() != (3 + tile_.numParameters() + 2)) {
    throw customException(
        "incorrect number of parameters provided");
  }

  // Tile
  double tile_params[tile_.numParameters()];
  for (int i = 3; i < tile_.numParameters() + 3; ++i) {
    tile_params[i - 3] = params[i];
  }
  tile_.setParameters(tile_params);

  setU0(params[3 + tile_.numParameters()]);
  setTileScale(params[3 + tile_.numParameters() + 1]);

  // Build boundary (need updated tile_control_points_) before computing scaled
  // coordinates
  buildBoundary_();

  // Motif - convert scaled to absolute coordinates
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

  if (ih_number == 64) {
    // Mirror line in the middle of tile. Tactile has it as a (fixed) vertical
    // line. Motif is assumed to have at least one mirror plane defined by
    // x-axis when at motif.theta_ = 0
    prefix = "d";
    induced = 1;  // S(P|M) = d1

    p0[0] = (tile_control_points_[1][0] - tile_control_points_[0][0]) / 2.0 +
            tile_control_points_[0][0];
    p0[1] = (tile_control_points_[1][1] + tile_control_points_[0][1]) / 2.0;

    p1[0] = (tile_control_points_[3][0] - tile_control_points_[2][0]) / 2.0 +
            tile_control_points_[2][0];
    p1[1] = (tile_control_points_[3][1] + tile_control_points_[2][1]) / 2.0;
  } else {
    throw(customException("unrecognized tile type"));
  }

  // Perform the revision
  results = revise_(p0, p1, orig_coords, motif_params->at(2), prefix, induced);

  // 3. Re-assign
  motif_params->at(0) = results[0];
  motif_params->at(1) = results[1];
  motif_params->at(2) = results[2];
}

const vector<double> Colloid::revise_(const vector<double>& p0,
                                      const vector<double>& p1,
                                      const vector<double>& orig_coords,
                                      const double current_theta,
                                      const string suffix, const int induced) {
  /**
   * Strategy:
   * Induced = c(n), place motif at rotation center, no forced rotation
   * necessary Induced = d(n>1), place motif at mirror intersection Induced =
   * d(n=1), place motif along the mirror line
   *
   * Proposal:
   * For d(n>1) could enforce that p0 should be provided such that this should
   * be the motif COM. Similarly, this convention could be helpful for c(n)
   * where p0 should provide the motif COM.
   */
  vector<double> unscaled_coords(2, 0), scaled_coords(2, 0);
  double absolute_theta = 0.0;

  if (suffix.compare("d") == 0) {
    // 1. Put motif COM on mirror lines
    if (induced == 1) {
      // If only 1 mirror line we have a DoF in terms of where on that line.
      unscaled_coords = project_to_line(p0, p1, orig_coords);
      scaled_coords = scale_coords_(unscaled_coords);
    } else {
      // With multiple mirrors, their intersection defines the COM of the motif
      throw customException("not implemented");
    }

    // 2. Round orientation to nearest allowable absolute theta value
    // Simply allow rotations that make mirror lines coincide (mod 2 leaves
    // motif unchanged)
    const int n = m_.symmetrySuffix(suffix);

    if ((n % induced != 0) ||
        (n < induced)) {  // n = 0 if no mirror so this catches that
      throw customException(
          "motif's reflection symmetry is incompatible with the tile");
    }
    const double curr_theta = thetaBounds(current_theta);
    double best_angle = 0.0, min_diff = pow(2.0 * M_PI, 2), diff = 0.0,
           angle = 0.0;
    for (int i = 0; i <= 2 * n;
         ++i) {  // Include 2*pi so rounding works correctly
      angle = i * M_PI / n;
      diff = pow(angle - curr_theta, 2);
      if (diff < min_diff) {
        min_diff = diff;
        if (i == (2 * n)) {
          best_angle = 0;  // Correct for 2*pi = 0
        } else {
          best_angle = angle;
        }
      }
    }
    absolute_theta = best_angle + tile_mirror_alignment(p0, p1);
  } else {
    throw customException("not implemented");

    // 1. Assign motif COM to rotation center

    // 2. No rotation is required
  }

  // 3. Update and return
  vector<double> results(3, 0.0);
  results[0] = scaled_coords[0];
  results[1] = scaled_coords[1];
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
   * {v0, v1, etc. for tile}, edge_u0, tile_scale]. See `setParameters()` for an
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
    params_.push_back(dummy[i]);  // Unscaled COM coordinates
  }
  rescaled = scale_coords_(params_);
  params_[0] = rescaled[0];  // scaled_com_x
  params_[1] = rescaled[1];  // scaled_com_y

  params_.push_back(dummy[2]);  // theta

  double tile_dummy[tile_.numParameters()];
  tile_.getParameters(tile_dummy);
  for (size_t i = 0; i < tile_.numParameters(); ++i) {
    params_.push_back(tile_dummy[i]);
  }
  params_.push_back(edge_u0_);
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

bool Colloid::isMotifInside(const int N = 20) {
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

  const double du = (1.0 - 0.0) / N;
  vector<int> boundary_ids;
  vector<vector<double>> polygon;
  vector<vector<double>> tile_control_points;

  perimeter_(0.0, du, N - 1, getTileScale(), &boundary_ids, &polygon,
             &tile_control_points);

  // Check each point in motif
  const vector<vector<double>> c = m_.getCoords();
  for (size_t i = 0; i < c.size(); ++i) {
    if (!pip(polygon, c[i])) {
      return false;
    }
  }

  return true;
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

  const double du = (1.0 - 0.0) / N;
  double outside = 0.0;

  vector<int> boundary_ids;
  vector<vector<double>> polygon;
  vector<vector<double>> tile_control_points;

  perimeter_(0.0, du, N - 1, getTileScale(), &boundary_ids, &polygon,
             &tile_control_points);

  // Check each point in motif
  const vector<vector<double>> c = m_.getCoords();
  for (size_t i = 0; i < c.size(); ++i) {
    if (!pip(polygon, c[i])) {
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

  perimeter_(edge_u0_, edge_du_, 1 + 4 + 1, getTileScale(), &boundary_ids_,
             &boundary_coords_, &tile_control_points_);
}

void Colloid::perimeter_(double u0, double du, int n, double scale,
                         vector<int>* boundary_ids,
                         vector<vector<double>>* boundary_coords,
                         vector<vector<double>>* tile_control_points) {
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
   * be in (0,1).
   * @param[in] du Gap along Bezier curve between boundary points. Should be
   * < 1.
   * @param[in] n Total number of points to place along each edge (incl. stop codons).
   * @param[in] scale The default Tactile tile is isotropically scaled by this
   * factor.
   * @param[out] boundary_ids Chemical identities (integers > 0) of boundary
   * points.  These are listed in ascending order from 1 up, without gap.
   * @param[out] boundary_coords Coordinates of points on tile's boundary.
   * @param[out] tile_control_points Control points from the ends of Bezier curves 
   * defining the edges. 
   */

  // Put points on perimeter (this returns stop codons AND control points)
  vector<vector<dvec2>> edges = perimeter_edges_(u0, du, n, scale);

  // Iterate over the edges of a single tile, asking the tiling to
  // tell you about the geometric information needed to transform
  // the edge shapes into position.  Note that this iteration is over
  // tiling "parts" which goes over entire edges for I and J, and each "half"
  // for S and U.
  vector<int> identity, cp_idx;
  vector<dvec2> shape;
  for (auto i : tile_.parts()) {
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
    } else {
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
  // order of the points.

  // Get control points at the "ends" of edges in counterclockwise order
  int iter = 0, edge_idx = 0;
  vector<vector<int>> cp_edges;
  for (auto i : tile_.parts()) {  // Iterator goes in counterclockwise edge order from Tactile
    vector<int> edges;
    if (i->getShape() == 1 || i->getShape() == 2) { // S or U edge
      if ((static_cast<int>(i->getId()) == edge_idx) && (iter > 0)) { // Should be the same as i->isSecondPart()?
        vector<int> idx = {cp_idx[iter-2], cp_idx[iter-1], cp_idx[iter], cp_idx[iter+1]};
        edges = unique_(shape, idx); // "Delete" the overlapping set of 2 to get edges
        assert(edges.size() == 2);
        cp_edges.push_back(edges); // Save the outer 2 points
      }
    } else { // J or I edge
      edges = {cp_idx[iter], cp_idx[iter+1]};
      cp_edges.push_back(edges); // Save the outer 2 points
    }
    edge_idx = static_cast<int>(i->getId());
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
        ++skip;
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
        if (d2 < eps) {
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

void Colloid::initMotif_(double max_scale_factor = 5.0,
                         double min_scale_factor = 0.2, int n_scale_incr = 100,
                         int N = 20, bool debug = false) {
  /**
   * Initialize the motif.
   *
   * This attempts to place the motif entirely inside the tile's boundary after
   * it has been initialized. The default tile settings from the Tactile library
   * are generally convenient and convex so it is recommended that these be kept
   * until this initialization is run. If a motif does not fit, then the tile is
   * expanded until it "just fits"; similarly, if it initially fits, the tile is
   * shrunk until it cannot be shrunk any more.
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

    buildBoundary_();  // Re-build based on final tile_scale_
  } else {
    // 1. Place motif COM and align
    // For certain IH tiles, the object must have certain conditions

    std::cout << "SKIPPING" << std::endl;
    // throw(customException("this tile type is not yet supported"));
  }

  return;
}

vector<vector<dvec2>> Colloid::perimeter_edges_(double u0, double du, int n,
                                                double scale) {
  /**
   * Create the Bezier curves that will serve as tile edges.
   *
   * These curves canonically go from (0,0) to (1,0) with some curvature between
   * them, and are uniformly scaled.  They will be rotated into place later on
   * by the Tactile library.
   *
   * @param u0 Starting point for boundary points along Bezier curve. Should be
   * in (0,1). S and U edges start from their center and essentially ignore this.
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
   * CP_right]. Points ascend from CP_left to CP_right.
   */

  // Create a vector to hold some edge shapes.  The tiling tells you
  // how many distinct edge shapes you need, but doesn't know anything
  // about how those shapes might be represented.  It simply assumes
  // that each one will be a curve from (0,0) to (1,0).  The tiling
  // provides tools to let you map those curves into position around
  // the outline of a tile.  All the curves below have exactly four
  // control points.
  vector<vector<dvec2>> edges;

  double u0_, du_;
  int n_;

  assert(n > 0);

  // Canonical Tactile coordinates
  dvec2 cp_left(0.0, 0.0), cp_right(1.0 * scale, 0.0), dummy;

  // Generate edge shapes.
  for (U8 idx = 0; idx < tile_.numEdgeShapes(); ++idx) {
    vector<dvec2> ej;

    // Define Bezier Curve that is sort of like a sphere impacting
    // This could definitely be changed, but will introduce more
    // free parameters.
    ej.push_back(dvec2(0.0, 0.0));
    ej.push_back(dvec2(1 / 3. * scale, sphere_deform_ * scale));
    ej.push_back(dvec2(2 / 3. * scale, sphere_deform_ * scale));
    ej.push_back(dvec2(1.0 * scale, 0.0));

    // Now, depending on the edge shape class, enforce symmetry
    // constraints on edges.
    switch (tile_.getEdgeShape(idx)) {
      case J:
        u0_ = u0;
        du_ = du;
        n_ = n;
        break;
      case U: // double dU since edge is half the size
        du_ = -2.0*du; // Walk "away" from midpoint
        u0_ = 1.0+du_/2.0; // Start U from x=1 = midpoint and walk away
        assert(n >= 2);
        n_ = (n-2)/2+1; // Only put 1 stop codon at the end

        // Reverse CP order
        dummy = cp_left;
        cp_left = cp_right;
        cp_right = dummy;
        break;
      case S: // double dU since edge is half the size
        du_ = -2.0*du; // Walk "away" from midpoint
        u0_ = 1.0+du_/2.0; // Start S from x=1 = midpoint and walk away
        assert(n >= 2);
        n_ = (n-2)/2+1; // Only put 1 stop codon at the end

        // Reverse CP order
        dummy = cp_left;
        cp_left = cp_right;
        cp_right = dummy;
        break;
      case I:
        ej[1].y = 0.0;
        ej[2].y = 0.0;
        u0_ = u0;
        du_ = du;
        n_ = n;
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
      setMotif(m);
    } catch (const exception& e) {
      throw(customException("unable to load Motif"));
    }
  }

  if (j.contains("Tile")) {
    try {
      IsohedralTiling t(j["Tile"]["ih_type"].get<int>());
      vector<double> p = j["Tile"]["parameters"].get<vector<double>>();
      double params[tile_.numParameters()];
      for (int i = 0; i < tile_.numParameters(); ++i) {
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
      edge_du_ = j["Properties"]["edge_du"].get<double>();
      edge_u0_ = j["Properties"]["edge_u0"].get<double>();
      setTileScale(j["Properties"]["tile_scale"].get<double>());
      sphere_deform_ = j["Properties"]["sphere_deform"].get<double>();
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
                    {"parameters", m_.getParameters()}};
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
                         {"sphere_deform", sphere_deform_},
                         {"edge_du", edge_du_},
                         {"edge_u0", edge_u0_},
                         {"tile_scale", getTileScale()},
                         {"tile_control_points", tile_control_points_},
                         {"params", params_},
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
                       vector<vector<double>>* box, const int nx = 1,
                       const int ny = 1) {
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

  // Get unit cell by combining the relevant "aspects"
  const int n = tile_.numAspects();
  vector<vector<double>> mc = m_.getCoords();
  vector<string> mt = m_.getTypes();
  glm::dvec2 t1 = tile_.getT1() * getTileScale(), t2 = tile_.getT2() * getTileScale();
  vector<dvec2> glm_coords;
  for (int i = 0; i < n; ++i) {
    dmat3 T = tile_.getAspectTransform(i);
    vector<string> t = m_.getTypes();
    for (size_t j = 0; j < mc.size(); ++j) {
      dvec3 c = T * dvec3(mc[j][0], mc[j][1], 0.0);
      glm_coords.push_back(c);
    }
  }

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
        coords->push_back(y);
        types->push_back(mt[i % mc.size()]);
      }
    }
  }
}
