/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 *
 * A colloid is composed of an isohedral tile boundary and internal motif.
 *
 * This file contains functions to build and manipulate colloids. Tile
 * information is inherited from the Tactile library, while a Motif
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
#include "src/tiling.hpp"
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
using std::string;
using std::vector;

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
  Colloid(Motif m, IsohedralTiling t, double tile_u0);
  ~Colloid();

  vector<double> unscale_coords_(const vector<double>& scaled_coords);
  vector<double> scale_coords_(const vector<double>& unscaled_coords);
  const vector<double> revise_(const vector<double>& p0,
                               const vector<double>& p1,
                               const vector<double>& orig_coords,
                               const double current_theta, const string suffix,
                               const int induced);

  void setParameters(const vector<double>& params);
  const vector<double> getParameters();

  void setMotif(Motif m);
  Motif getMotif();

  void setTile(IsohedralTiling t);
  const IsohedralTiling getTile();
  void setTileScale(const double s) {
    tile_scale_ = s;
  }  // Assign the tile_scale_.
  double tileArea();

  bool isMotifInside(const int N);
  double fractionMotifInside(const int N);

  bool isTileFundamental();

  void init() {
    buildBoundary_();
    initMotif_(5.0, 0.2, 1000, 20);
    built_ = true;
  }  // Initialize the colloid.

  void load(const string filename);
  void dump(const string filename);
  void dumpXYZ(const string filename, const bool full);

  vector<double> boundaryCOM();

  void setU0(const double u0) { edge_u0_ = u0; }
  void setDU(const double du) { edge_du_ = du; }

  void unitCell(vector<vector<double>>* coords, vector<string>* types,
                vector<vector<double>>* box, const int nx, const int ny);

 private:
  void defaults_();
  void buildBoundary_();
  void constrain_(vector<double>* motif_params);
  void initMotif_(double max_scale_factor, double min_scale_factor,
                  int n_scale_incr, int N);
  void perimeter_(double u0, double du, int n, double scale,
                  vector<int>* boundary_ids,
                  vector<vector<double>>* boundary_coords,
                  vector<vector<double>>* tile_control_points);
  vector<vector<dvec2>> perimeter_edges_(double u0, double du, int n,
                                         double scale);

  bool tile_assigned_;   // Has the tile been assigned yet?
  bool motif_assigned_;  // Has the motif been assigned yet?
  bool built_;           // Has the colloid been constructed at least once?

  double sphere_deform_;  // Normalized amount a sphere "deforms" the edge.
  double edge_du_;        // Parameterized (Bezier) gap between boundary points.
  double
      edge_u0_;  // Parameterized (Bezier) starting point for boundary points.
  double tile_scale_;  // The default Tactile tile is isotropically scaled by
                       // this factor.

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

#endif  // SRC_COLLOID_HPP_
