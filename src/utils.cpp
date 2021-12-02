/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 */

#include "utils.hpp"

void dumpXYZ(const vector<vector<double>> &coords, const vector<string> &types, const string filename) {
  /**
   * Dump (x,y) coordinates to an XYZ file. 
   *
   * @param coords (x,y) coordinates of particles in unit cell.
   * @param types Chemical type of each particle.
   * @param filename Name of file to write to. Will overwrite by default.
   *
   * @throws customException if anything goes wrong.
   */

  try {
    ofstream xyz(filename);

    xyz << coords.size() << endl;
    xyz << endl;

    for (size_t i = 0; i < coords.size(); ++i) {
      xyz << types[i] << "\t" << coords[i][0] << "\t"
          << coords[i][1] << "\t" << 0 << endl;
    }
  } catch (const exception& e) {
    throw(customException("unable to write to xyz file"));
  }

  return;
}
