/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 *
 * This file contains functions to build and manipulate colloids. 
 *
 * A colloid is composed of an isohedral tile boundary and internal motif.
 * Tile information is inherited from the Tactile library, while a Motif
 * class is described in motif.hpp. This should be self-contained so
 * that, for example, external optimizers can manipulate the colloid's
 * parameters controlling its shape, etc. and compute properties from
 * that.
 */

#ifndef SRC_COLLOID_HPP_
#define SRC_COLLOID_HPP_

#include <bits/stdc++.h>

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "src/json.hpp"
#include "src/motif.hpp"
#include "tiling.hpp" // Cannot provide absolute reference due to complications during compile time with Tactile library files
#include "src/utils.hpp"

using glm::dmat2;
using glm::dmat3;
using glm::dvec2;
using glm::dvec3;
using glm::inverse;
using glm::transpose;

using std::ifstream;
using std::map;
using std::ofstream;
using std::stringstream;
using std::string;
using std::vector;
using std::to_string;

using csk::I;
using csk::IsohedralTiling;
using csk::J;
using csk::S;
using csk::U;
using csk::U8;

using json = nlohmann::json;

class Colloid {
 public:
  Colloid();
  Colloid(Motif m, IsohedralTiling t, vector<double> tile_u0, vector<double> edge_df, const double du = 0.1, bool debug=false, double min_angle = 0.0, double min_gap = 0.0);
  ~Colloid();

  vector<double> unscale_coords_(const vector<double>& scaled_coords);
  vector<double> scale_coords_(const vector<double>& unscaled_coords);
  const vector<double> revise_motif_params_(const vector<double>& p0,
                               const vector<double>& p1,
                               const vector<double>& orig_coords,
                               const double current_theta, const string suffix,
                               const int induced);

  void setParameters(const vector<double>& params, const double df_min = 0.1);
  const vector<double> getParameters();

  void setMinAngle(const double theta);
  double getMinAngle() const { return min_angle_; }

  void setMinGap(const double gap);
  double getMinGap() const { return min_gap_; }

  void setMotif(Motif m);
  Motif getMotif();

  void init(bool debug=false);

  void setTile(IsohedralTiling t);
  const IsohedralTiling getTile();
  void setTileScale(const double s) {
    tile_scale_ = s;
  }  // Assign the tile_scale_.
  double getTileScale() const { return tile_scale_; }
  double tileArea();

  bool isMotifInside(const int N);
  double fractionMotifInside(const int N = 20, const double skin = 0.0);

  bool isTileFundamental();

  void load(const string filename);
  void dump(const string filename);
  void dumpXYZ(const string filename, const bool full);

  vector<double> boundaryCOM();

  void setU0(const vector<double> u0) { edge_u0_ = u0; }
  const vector<double> getU0();

  void setDform(const vector<double> df) { edge_df_ = df; }
  const vector<double> getDform();

  void setDU(const double du) { edge_du_ = du; }
  const double getDU();

  void unitCell(vector<vector<double>>* coords, vector<string>* types,
                vector<vector<double>>* box, const int nx = 1, const int ny = 1, 
                const double tol = 1.0e-8, const bool unique = false,
                const bool boundary = false);

  vector<vector<double>> getBoundaryCoords() const { return boundary_coords_; }
  vector<int> getBoundaryIds() const { return boundary_ids_; }
  vector<vector<double>> getTileControlPoints() const { return tile_control_points_; }

  const int countIntersections(const int N = 10, int* bad_angles = NULL, double* min_gap = NULL);
  const vector<vector<dvec2>> buildTilePolygon(const int N);

 private:
  void defaults_();
  void buildBoundary_();
  void constrain_(vector<double>* motif_params);
  void initMotif_(double max_scale_factor, double min_scale_factor,
                  int n_scale_incr, int N, bool debug);
  void perimeter_(vector<double> u0, vector<double> df, double du, int n, double scale,
                  vector<int>* boundary_ids,
                  vector<vector<double>>* boundary_coords,
                  vector<vector<double>>* tile_control_points,
                  vector<vector<dvec2>>* polygon);
  vector<vector<dvec2>> perimeter_edges_(vector<double> u0, vector<double> df, double du, int n,
                                         double scale);

  bool tile_assigned_;   // Has the tile been assigned yet?
  bool motif_assigned_;  // Has the motif been assigned yet?
  bool built_;           // Has the colloid been constructed at least once?

  vector<double> edge_df_; // How much each edge is deflected (Bezier)
  vector<double>
      edge_u0_;  // Parameterized (Bezier) starting point for boundary points.

  double edge_du_;       // Parameterized (Bezier) gap between boundary points.
  double tile_scale_;  // The default Tactile tile is isotropically scaled by
                       // this factor.
  double min_angle_;   // Smallest angle the tile's boundary is allow to have.
  double min_gap_;     // Smallest distance allowable between points on non-adjacent edges

  vector<int> boundary_ids_;  // Chemical identities of boundary points.
  vector<double> params_;     // Unrolled parameter vector.

  vector<vector<double>>
      boundary_coords_;  // Coordinates of points on tile's boundary.
  vector<vector<double>>
      tile_control_points_;  // Control points on Bezier curves which are tile
                             // vertices.

  IsohedralTiling tile_;  // Isohedral tile from Tactile library.
  Motif m_;               // The colloid's motif.
};

vector<int> unique_(const vector<dvec2>& shape, const vector<int>& idx, const double eps = 1.0e-12);

void new_idx_(const vector<dvec2>& shape, vector<int>& last_edge, vector<int>& curr_edge, unsigned int& ov_, unsigned int& new_, const double eps = 1.0e-12);

bool is_undirected_mirror(const int ih_number, const int edge_idx) ;

#endif  // SRC_COLLOID_HPP_
