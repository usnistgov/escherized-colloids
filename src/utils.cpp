/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 */

#include "src/utils.hpp"

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
  while (eq_theta < 0) { // If < 0, make it positive first
    eq_theta += 2*M_PI;
  }
  return fmod(eq_theta, 2.0*M_PI);
}

const int commonVertices(const dvec2 p1, const dvec2 q1, const dvec2 p2, const dvec2 q2, const double tol) {
  /**
   * Count the number of points (vertices) that are the same between different line segments.
   */
  int count = 0;
  vector<double> d2;

  d2.push_back(pow(p1.x-p2.x, 2) + pow(p1.y-p2.y, 2)); // p1 vs. p2
  d2.push_back(pow(p1.x-q2.x, 2) + pow(p1.y-q2.y, 2)); // p1 vs. q2

  d2.push_back(pow(q1.x-p2.x, 2) + pow(q1.y-p2.y, 2)); // q1 vs. p2
  d2.push_back(pow(q1.x-q2.x, 2) + pow(q1.y-q2.y, 2)); // q1 vs. q2

  for (unsigned int i=0; i < d2.size(); ++i) {
    if (d2[i] < tol*tol) {
      count++;
    }
  }

  return count;
}

const bool onSegment(const dvec2 p, const dvec2 q, const dvec2 r) {
  /**
   * Given three collinear points p, q, r, the function checks if
   * point q lies on line segment 'pr'. Thi is from:
   * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/.
   */
  if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y)) {
     return true;
  }

  return false;
}

const int orientation(const dvec2 p, const dvec2 q, const dvec2 r) {
  /**
   * To find orientation of ordered triplet (p, q, r).
   * The function returns following values:
   * 0 --> p, q and r are collinear
   * 1 --> Clockwise
   * 2 --> Counterclockwise
   * This is from: https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/.
   */

  // See https://www.geeksforgeeks.org/orientation-3-ordered-points/
  // for details of below formula.
  double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

  if (fabs(val) < 1.0e-12) {
    return 0;  // collinear
  }

  return (val > 0)? 1: 2; // clock- or counterclock- wise
}

const bool doIntersect(const dvec2 p1, const dvec2 q1, const dvec2 p2, const dvec2 q2) {
  /**
   * The main function that returns true if line segment 'p1q1'
   * and 'p2q2' intersect. This is from:
   * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/.
   */

  // Find the four orientations needed for general and
  // special cases
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);

  // General case
  if (o1 != o2 && o3 != o4)
      return true;

  // Special Cases
  // p1, q1 and p2 are collinear and p2 lies on segment p1q1
  if (o1 == 0 && onSegment(p1, p2, q1)) {
    return true;
  }

  // p1, q1 and q2 are collinear and q2 lies on segment p1q1
  if (o2 == 0 && onSegment(p1, q2, q1)) {
    return true;
  }

  // p2, q2 and p1 are collinear and p1 lies on segment p2q2
  if (o3 == 0 && onSegment(p2, p1, q2)) {
    return true;
  }

  // p2, q2 and q1 are collinear and q1 lies on segment p2q2
  if (o4 == 0 && onSegment(p2, q1, q2)) {
    return true;
  }

  return false; // Doesn't fall in any of the above cases
}

const bool isInside(const vector<vector<dvec2>>& polygon, const dvec2& point) {
  /**
   * Check if a point is inside a polygon.
   * Adapted from: https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
   *
   * @param polygon Vector of edges, where each edge has a vector of points.
   * Edges are not in a fixed orientational order, but points along an edge
   * are ordered.
   * @param point Point to check.
   */

  assert(polygon.size() >= 3); // There must be at least 3 edges in polygon
  const dvec2 extreme = {PIP_INF, point.y}; // Create a point for line segment from point to infinite

  int count = 0;
  for (unsigned int edge_idx=0; edge_idx < polygon.size(); ++edge_idx) {
    for (unsigned int i=0; i < polygon[edge_idx].size()-1; ++i) {
      // Check if the line segment from 'point' to 'extreme' intersects
      // with the line segment from 'polygon[edge][i]' to 'polygon[edge][i+1]'.
      // Note that polygon's edge points include "ends" (Bezier 0 -> 1).
      if (doIntersect(polygon[edge_idx][i], polygon[edge_idx][i+1], point, extreme)) {
        // If the point 'point' is collinear with line segment 'i->(i+1)',
        // then check if it lies on segment. If it lies, return true,
        // otherwise false.
        if (orientation(polygon[edge_idx][i], point, polygon[edge_idx][i+1]) == 0) {
          return onSegment(polygon[edge_idx][i], point, polygon[edge_idx][i+1]);
        }

        count++;
      }
    }
  }

  return (count%2 == 1); // Return true if count is odd, false otherwise
}
