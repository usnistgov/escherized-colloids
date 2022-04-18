#ifndef TEST_NONFDTILES_HPP_
#define TEST_NONFDTILES_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "tiling.hpp"
#include "../src/colloid.hpp"

using namespace std;
using namespace csk;

class NonFundamentalTest : public ::testing::Test {
 protected:
  void SetUp() override {
    prefix = "./fd_ref/";

  }

  vector<int> ih_type;
  string prefix;

  Motif m;
  IsohedralTiling *t;
  Colloid *c;

  void TearDown() override {
    delete t;
    delete c;
  }
};

TEST_F(NonFundamentalTest, IH64) {
  // IH 64 has S(P) = pm, S(P|M) = d1
  t = new IsohedralTiling(64);
  c = new Colloid();

  // Bad motif symmetry
  try {
    m.load("../motif_library/c1_random.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's reflection symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d1_vitruvian.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  
  // Double check control points are numbered as expected (counterclockwise)
  vector<vector<double>> cp = c->getTileControlPoints();
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]), p3(cp[3][0], cp[3][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p3), 2);
  ASSERT_EQ(orientation(p2, p3, p0), 2);
  ASSERT_EQ(orientation(p3, p0, p1), 2);

  ASSERT_FLOAT_EQ(p0.x, p3.x); // p0 is top right of rectangle
  ASSERT_FLOAT_EQ(p0.y, p1.y);
  ASSERT_FLOAT_EQ(p2.x, p1.x); // p2 is bottom left of rectangle
  ASSERT_FLOAT_EQ(p2.y, p3.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // centers to vertical line at center of tile
  ASSERT_FLOAT_EQ(new_params[1], 1); // doesn't affect y

  // Try to rotate
  params = c->getParameters();

  params[2] = 0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = 2*M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3.*M_PI/2.);

  params[2] = 3*M_PI/2.;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3.*M_PI/2.);

  params[2] = 2*M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3.*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH12) {
  // IH 12 has S(P) = cm, S(P|M) = d1
  t = new IsohedralTiling(12);
  c = new Colloid();

  // Bad motif symmetry
  try {
    m.load("../motif_library/c1_random.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's reflection symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d1_vitruvian.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  
  // Double check control points are numbered as expected (counterclockwise)
  vector<vector<double>> cp = c->getTileControlPoints();
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]), p3(cp[3][0], cp[3][1]), p4(cp[4][0], cp[4][1]), p5(cp[5][0], cp[5][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p3), 2);
  ASSERT_EQ(orientation(p2, p3, p4), 2);
  ASSERT_EQ(orientation(p3, p4, p5), 2);
  ASSERT_EQ(orientation(p4, p5, p0), 2);
  ASSERT_EQ(orientation(p5, p0, p1), 2);

  ASSERT_FLOAT_EQ(p1.x, p5.x); 
  ASSERT_FLOAT_EQ(p1.y, p2.y);
  ASSERT_FLOAT_EQ(p2.x, p4.x);
  ASSERT_FLOAT_EQ(p4.y, p5.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // centers to vertical line at center of tile
  ASSERT_FLOAT_EQ(new_params[1], 1); // doesn't affect y

  // Try to rotate
  params = c->getParameters();

  params[2] = 0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = 2*M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3.*M_PI/2.);

  params[2] = 3*M_PI/2.;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3.*M_PI/2.);

  params[2] = 2*M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3.*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH14) {
  // IH 14 has S(P) = cm, S(P|M) = d1
  t = new IsohedralTiling(14);
  c = new Colloid();

  // Bad motif symmetry
  try {
    m.load("../motif_library/c1_random.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's reflection symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d1_vitruvian.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  
  // Double check control points are numbered as expected (counterclockwise)
  vector<vector<double>> cp = c->getTileControlPoints();
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]), p3(cp[3][0], cp[3][1]), p4(cp[4][0], cp[4][1]), p5(cp[5][0], cp[5][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p3), 2);
  ASSERT_EQ(orientation(p2, p3, p4), 2);
  ASSERT_EQ(orientation(p3, p4, p5), 2);
  ASSERT_EQ(orientation(p4, p5, p0), 2);
  ASSERT_EQ(orientation(p5, p0, p1), 2);

  ASSERT_FLOAT_EQ(p5.x, p3.x); 
  ASSERT_FLOAT_EQ(p5.y, p0.y);
  ASSERT_FLOAT_EQ(p0.x, p2.x);
  ASSERT_FLOAT_EQ(p2.y, p3.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 1); // doesn't affect x
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // centers to horizontal line at center of tile

  // Try to rotate
  params = c->getParameters();

  params[2] = 0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 0);

  params[2] = M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 0);

  params[2] = M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI);

  params[2] = 3*M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI);

  params[2] = 3*M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 0);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 0);
}

#endif  // TEST_NONFDTILES_HPP_
