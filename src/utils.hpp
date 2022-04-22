/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 *
 * General utilities and variables with global scope.
 */

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_

#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "tiling.hpp" // Cannot provide absolute reference due to complications during compile time with Tactile library files

using std::endl;
using std::exception;
using std::ofstream;
using std::string;
using std::vector;
using std::min;
using std::max;

using glm::dvec2;

#define EXIT_FAILURE 1
#define PIP_INF 1000000
#define INF_SYMM 180 // 180 chosen so that when d_inf happens, we round allowed_angles to integer degrees

// Isohedral tile types that are fundamental domains of crystals.
const int FD_TYPES[46] = {1,  2,  3,  4,  5,  6,  7,  21, 22, 23, 24, 25,
                          27, 28, 30, 31, 33, 38, 39, 41, 42, 43, 44, 45,
                          46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 77,
                          78, 79, 80, 81, 83, 84, 85, 86, 87, 88};

class
    customException {  // Adapted from
                       // https://www.oreilly.com/library/view/c-cookbook/0596007612/ch09s02.html
 public:
  explicit customException(const string& msg) : msg_(msg) {}
  ~customException() {}

  string getMessage() const { return (msg_); }

 private:
  string msg_;
};

void dumpXYZ(const vector<vector<double>>& coords, const vector<string>& types,
             const string filename);

vector<double> project_to_line(const vector<double>& p0,
                               const vector<double>& p1,
                               const vector<double>& coords);

const double tile_mirror_alignment(const vector<double>& p0,
                                   const vector<double>& p1);

double thetaBounds(const double theta);

const int commonVertices(const dvec2 p1, const dvec2 q1, const dvec2 p2, const dvec2 q2, const double tol = 1.0e-12);
const bool onSegment(const dvec2 p, const dvec2 q, const dvec2 r);
const int orientation(const dvec2 p, const dvec2 q, const dvec2 r);
const bool doIntersect(const dvec2 p1, const dvec2 q1, const dvec2 p2, const dvec2 q2);
const bool isInside(const vector<vector<dvec2>>& polygon, const dvec2& point);

#endif  // SRC_UTILS_HPP_
