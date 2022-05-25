/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 *
 * This file contains functions to build and manipulate motifs.
 *
 * A motif for the colloid is what is "inside" the tile.
 */

#ifndef SRC_MOTIF_HPP_
#define SRC_MOTIF_HPP_

#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "src/json.hpp"
#include "src/utils.hpp"

using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::string;
using std::vector;
using std::stringstream;

using json = nlohmann::json;

class Motif {
 public:
  Motif();
  ~Motif();

  Motif(Motif &other);
  void copy(Motif &other);

  void loadXYZ(const string filename);
  void dumpXYZ(const string filename);

  void setCoords(const vector<vector<double>> &coords, double theta);
  const vector<vector<double>> getCoords() const { return coords_; }

  void setTypes(const vector<string> &types);
  const vector<string> getTypes() const { return types_; }

  void setParameters(const vector<double> &params);
  const vector<double> getParameters();

  vector<double> getCOM() {
    computeCOM_();
    return com_;
  }
  void rotate(const double theta);
  void translate(const vector<double> &dx);

  void load(const string filename);

  void setSymmetry(const string s) { symmetry_ = s; }
  const string getSymmetry() const { return symmetry_; }
  const int symmetrySuffix(const string prefix);

  const double minDistance();

 private:
  void computeCOM_();

  double theta_;  // Absolute orientation (right handed, counterclockwise
                  // convention).
  string
      symmetry_;  // Point symmetry the motif has (must be determined by user).
  vector<string> types_;  // Chemical identity of each particle in the motif.
  vector<double> com_;    // (x,y) center of mass.

  vector<vector<double>>
      coords_;  // (x,y) coordinates of motif's constituent particles.
};


const int compatibility_check(Motif& m, const string prefix, const int induced);

#endif  // SRC_MOTIF_HPP_
