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

vector<double> project_to_line(const vector<double>& p0, const vector<double>& p1, const vector<double>& coords) {
  /**
   * Orthogonal projection of a point to a line.
   *
   * @param p0 One point on the line to project onto.
   * @param p1 Another point on the line to project onto.
   * @param coords The original point to project.
   *
   * @returns Vector of projected (x, y) coordinates.
   */

  // https://math.stackexchange.com/questions/62633/orthogonal-projection-of-a-point-onto-a-line
  // https://openstax.org/books/intermediate-algebra-2e/pages/4-6-solve-systems-of-equations-using-determinants

  const double eps = 1.0e-12;
  vector<double> res(2, 0.0);
  const double dx = (p1[0] - p0[0]), dy = (p1[1] - p0[1]);

  if ((std::abs(dx) < eps) && (std::abs(dy) < eps)) {
    throw(customException("points are identical"));
  } else if (std::abs(dx) < eps) {
    // Line is vertical
    res[0] = (p1[0] + p0[0])/2.0; // Should be essentially identical
    res[1] = coords[1]; // y-coordinate unaffected
  } else if (std::abs(dy) < eps) {
    // Line is horizontal 
    res[0] = coords[0]; // x-coordinate unaffected
    res[1] = (p1[1] + p0[1])/2.0; // Should be essentially identical
  } else {
    const double m1 = dy / dx, b1 = p0[1] - m1*p0[0];
    const double m2 = -1.0/m1, b2 = coords[1] - coords[0]*m2;
    const double D = (-m1 + m2), Dx = (b1 - b2), Dy = (-b2*m1 + m2*b1);
    res[0] = Dx/D;
    res[1] = Dy/D;
  }

  return res;
}

const double tile_mirror_alignment(const vector<double>& p0, const vector<double>& p1) {
  /**
   * Compute the counterclockwise angle a mirror forms with the x-axis.
   *
   * This mirror is defined by two points which belong to the line. Note that there
   * are 2 angles (this and another +/- pi radians) that would work as well because 
   * by definition the mirror line is a line.  Either angle should suffice.
   *
   * @param p0 One point on the mirror line.
   * @param p1 Another point on the mirror line.
   *
   * @returns theta Counterclockwise angle the mirror forms with the x-axis.
   */
  const double dx = (p1[0] - p0[0]), dy = (p1[1] - p0[1]);
  return thetaBounds(atan2(dy, dx));
}

double thetaBounds(const double theta) {
  /*
   * Bound theta inside [0, 2*pi).
   *
   * @param theta Angle.
   *
   * @return Equivalent angle
   */
  double eq_theta = theta;
  while(eq_theta < 0) { // If < 0, make it positive first
    eq_theta += 2*M_PI;
  }
  return fmod(eq_theta, 2.0*M_PI);
}
