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

using std::endl;
using std::exception;
using std::ofstream;
using std::string;
using std::vector;

#define EXIT_FAILURE 1

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

#endif  // SRC_UTILS_HPP_
