/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 */

#include "colloid.hpp"

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

Colloid::Colloid(Motif m, IsohedralTiling t, double tile_u0)
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
  init();  // This will find a good value for tile_scale_
}

void Colloid::defaults_() {
  /**
   * Assign default paramter values.
   */

  tile_assigned_ = false;
  motif_assigned_ = false;
  built_ = false;
  tile_scale_ = 1.0;
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

  // Tile
  bool fundamental = isTileFundamental();
  int adjust = 0;
  if (!fundamental) {
    adjust = 1;
  }
  double tile_params[tile_.numParameters()];
  for (int i = 3 - 3 * adjust; i < tile_.numParameters() + 3 - 3 * adjust;
       ++i) {
    tile_params[i - 3 + 3 * adjust] = params[i + 1 * adjust];
  }
  tile_.setParameters(tile_params);

  setU0(params[3 - 2 * adjust + tile_.numParameters()]);
  tile_scale_ = params[3 - 2 * adjust + tile_.numParameters() + 1];

  // Build boundary (need tile_control_points_) before computing scaled
  // coordinates
  buildBoundary_();

  // Motif - convert scaled to absolute coordinates
  if (fundamental) {
    const vector<double> scaled_coords = {params[0], params[1]};
    const vector<double> us = unscale_coords_(scaled_coords);
    const vector<double> motif_params = {us[0], us[1], params[2]};
    m_.setParameters(motif_params);
  } else {
    // Will have to update motif's COM and angle according to specifics of each
    // tile.

    // Logic goes here ...

    throw(customException("non-fundamental tiles not yet supported"));
  }
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
  if (isTileFundamental()) {  // If tile is FD motif COM coordinates are valid
                              // DoFs
    for (size_t i = 0; i < 2; ++i) {
      params_.push_back(dummy[i]);  // Unscaled COM coordinates
    }
    rescaled = scale_coords_(params_);
    params_[0] = rescaled[0];  // scaled_com_x
    params_[1] = rescaled[1];  // scaled_com_y
  }
  params_.push_back(dummy[2]);  // theta

  double tile_dummy[tile_.numParameters()];
  tile_.getParameters(tile_dummy);
  for (size_t i = 0; i < tile_.numParameters(); ++i) {
    params_.push_back(tile_dummy[i]);
  }

  params_.push_back(edge_u0_);
  params_.push_back(tile_scale_);

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

const Motif Colloid::getMotif() {
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

  perimeter_(0.0, du, N - 1, tile_scale_, &boundary_ids, &polygon,
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

  perimeter_(0.0, du, N - 1, tile_scale_, &boundary_ids, &polygon,
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

  perimeter_(edge_u0_, edge_du_, 1 + 4 + 1, tile_scale_, &boundary_ids_,
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
   * library. The curve's control points are adjusted to be consistent with
   * symmetry. Each edge assumed to be "impacted" by a sphere to create a
   * curvature of a single sign along its edge. The identities of each
   * symmetrically unique point along the bounday are identified in ascending
   * order, without gap. 0 is used for stop codons, positive integers for
   * interacting points.
   *
   * @param[in] u0 Starting point for boundary points along Bezier curve. Should
   * be in (0,1).
   * @param[in] du Gap along Bezier curve between boundary points. Should be
   * < 1.
   * @param[in] n Number of points to place along each edge.
   * @param[in] scale The default Tactile tile is isotropically scaled by this
   * factor.
   * @param[out] boundary_ids Chemical identities (integers > 0) of boundary
   * points.
   * @param[out] boundary_coords Coordinates of points on tile's boundary.
   * @param[out] tile_control_points Control points on Bezier curves defining
   * edges.
   */

  vector<vector<dvec2>> edges = perimeter_edges_(u0, du, n, scale);

  // Use a vector to hold the control points of the final tile outline.
  vector<dvec2> shape;

  // Iterate over the edges of a single tile, asking the tiling to
  // tell you about the geometric information needed to transform
  // the edge shapes into position.  Note that this iteration is over
  // whole tiling edges.  It's also to iterator over partial edges
  // (i.e., halves of U and S edges) using t.parts() instead of t.shape().
  vector<int> identity;
  int total_edges = 0;
  for (auto i : tile_.shape()) {
    // 1 + 4 colors + 1 stop codon
    vector<int> pattern = {0,
                           i->getId() * 4 + 1,
                           i->getId() * 4 + 2,
                           i->getId() * 4 + 3,
                           i->getId() * 4 + 4,
                           0};  // 0 for stop codons

    // Interacting points on U and S edges are centro-symmetrically labelled
    if (i->getShape() == 1 || i->getShape() == 2) {
      pattern[3] = pattern[2];
      pattern[4] = pattern[1];
    }

    // Get the relevant edge shape created above using i->getId().
    const vector<dvec2> ed = edges[i->getId()];

    // Also get the transform that maps to the line joining consecutive
    // tiling vertices.
    const glm::dmat3& T = i->getTransform();

    // If i->isReversed() is true, we need to run the parameterization
    // of the path backwards.
    if (i->isReversed()) {
      for (size_t idx = 1; idx < ed.size(); ++idx) {
        shape.push_back(T * dvec3(ed[ed.size() - 1 - idx], 1.0 * scale));
      }
      for (size_t idx = 0; idx < pattern.size(); ++idx) {
        identity.push_back(pattern[pattern.size() - 1 - idx]);
      }
      identity.push_back(-1);  // -1 for control points
    } else {
      for (size_t idx = 1; idx < ed.size(); ++idx) {
        shape.push_back(T * dvec3(ed[idx], 1.0 * scale));
      }
      for (size_t idx = 0; idx < pattern.size(); ++idx) {
        identity.push_back(pattern[idx]);
      }
      identity.push_back(-1);  // -1 for control points
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

  boundary_ids->clear();
  boundary_coords->clear();
  tile_control_points->clear();

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
    } else {
      // Save control points to compute area, etc.
      tile_control_points->push_back(c);
    }
  }
}

void Colloid::initMotif_(double max_scale_factor = 5.0,
                         double min_scale_factor = 0.2, int n_scale_incr = 100,
                         int N = 20) {
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
   *
   * @throws customException if initialization fails for any reason.
   */

  if (isTileFundamental()) {
    // If no restrictions, place motif COM on tile COM and shrinkwrap
    // Tactile default parameters provide convenient, convex, starting points
    // at least one control point is at (0,0).

    // 1. Place motif so its COM is on tile COM
    // tile COM is estimated from the functional points on the perimeter
    vector<double> motif_com = m_.getCOM(), tile_com = boundaryCOM(), dx;
    for (size_t i = 0; i < motif_com.size(); ++i) {
      dx.push_back(tile_com[i] - motif_com[i]);
    }
    m_.translate(dx);

    // 2. Expand/contract the tile until motif "just" fits
    double min_scale = tile_scale_ * min_scale_factor;
    double max_scale = tile_scale_ * max_scale_factor;
    double orig_scale = tile_scale_;
    double last_scale = tile_scale_;
    if (isMotifInside(N)) {
      // Shrink the tile to fit
      bool found = false;
      for (int i = 0; i <= n_scale_incr; ++i) {
        tile_scale_ =
            (orig_scale - (orig_scale - min_scale) / n_scale_incr * i);
        vector<double> dd = {(tile_scale_ / last_scale - 1.0) * tile_com[0],
                             (tile_scale_ / last_scale - 1.0) * tile_com[1]};
        m_.translate(dd);

        if (!isMotifInside(N)) {
          // Move motif back to last position
          dd[0] = -dd[0];
          dd[1] = -dd[1];
          m_.translate(dd);

          // Save last good scale
          tile_scale_ = last_scale;
          found = true;
          break;
        }
        last_scale = tile_scale_;
      }
      if (!found) {
        throw(customException("unable to initialize motif inside tile"));
      }
    } else {
      // Expand the tile to fit
      bool found = false;
      for (int i = 0; i <= n_scale_incr; ++i) {
        tile_scale_ =
            (orig_scale + (max_scale - orig_scale) / n_scale_incr * i);
        vector<double> dd = {(tile_scale_ / last_scale - 1.0) * tile_com[0],
                             (tile_scale_ / last_scale - 1.0) * tile_com[1]};
        m_.translate(dd);
        if (isMotifInside(N)) {
          // Motif and tile now in acceptable positions
          found = true;
          break;
        }
        last_scale = tile_scale_;
      }
      if (!found) {
        throw(customException("unable to initialize motif inside tile"));
      }
    }

    buildBoundary_();  // Re-build based on final tile_scale_
  } else {
    // 1. Place motif COM and align
    // For certain IH tiles, the object must have certain conditions

    throw(customException("this tile type is not yet supported"));
  }

  return;
}

vector<vector<dvec2>> Colloid::perimeter_edges_(double u0, double du, int n,
                                                double scale) {
  /**
   * Create the Bezier curves that will serve as tile edges.
   *
   * These curves are generally from (0,0) to (1,0) with some curvature between
   * them, and are uniformly scaled.  They will be rotated into place later on
   * by the Tactile library.
   *
   * @param u0 Starting point for boundary points along Bezier curve. Should be
   * in (0,1).
   * @param du Gap along Bezier curve between boundary points. Should be < 1.
   * @param n Number of points to place along each edge.
   * @param scale The default Tactile tile is isotropically scaled by this
   * factor.
   *
   * @returns Discretized points along each edge, including Bezier control
   * points at the ends.
   */

  // Create a vector to hold some edge shapes.  The tiling tells you
  // how many distinct edge shapes you need, but doesn't know anything
  // about how those shapes might be represented.  It simply assumes
  // that each one will be a curve from (0,0) to (1,0).  The tiling
  // provides tools to let you map those curves into position around
  // the outline of a tile.  All the curves below have exactly four
  // control points.
  vector<vector<dvec2>> edges;

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
        break;
      case U:
        ej[2].x = 1.0 * scale - ej[1].x;
        ej[2].y = ej[1].y;
        break;
      case S:
        ej[2].x = 1.0 * scale - ej[1].x;
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
    coords.push_back(dvec2(0.0, 0.0));
    for (int k = 0; k < n; ++k) {
      coords.push_back(bezier(ej[0], ej[1], ej[2], ej[3], u0 + k * du));
    }
    coords.push_back(dvec2(1.0 * scale, 0));

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
      tile_scale_ = j["Properties"]["tile_scale"].get<double>();
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
                         {"tile_scale", tile_scale_},
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
