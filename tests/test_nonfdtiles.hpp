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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

TEST_F(NonFundamentalTest, IH68) {
  // IH 68 has S(P) = cm, S(P|M) = d1
  t = new IsohedralTiling(68);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p0.y, p2.y); 
  ASSERT_FLOAT_EQ(p1.x, p3.x);

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

TEST_F(NonFundamentalTest, IH13) {
  // IH 13 has S(P) = pmg, S(P|M) = d1
  t = new IsohedralTiling(13);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p5.x, p1.x); 
  ASSERT_FLOAT_EQ(p5.y, p4.y);
  ASSERT_FLOAT_EQ(p4.x, p2.x);
  ASSERT_FLOAT_EQ(p2.y, p1.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[1], 1); // doesn't affect y
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // centers to vertical line at center of tile

  // Try to rotate
  params = c->getParameters();

  params[2] = 0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH15) {
  // IH 15 has S(P) = cm, S(P|M) = d1
  t = new IsohedralTiling(15);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

TEST_F(NonFundamentalTest, IH66) {
  // IH 66 has S(P) = pmg, S(P|M) = d1
  t = new IsohedralTiling(66);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

TEST_F(NonFundamentalTest, IH69) {
  // IH 69 has S(P) = pmg, S(P|M) = d1
  t = new IsohedralTiling(69);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p0.x, p2.x); 
  ASSERT_FLOAT_EQ(p1.y, p3.y);

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

TEST_F(NonFundamentalTest, IH26) {
  // IH 26 has S(P) = cmm, S(P|M) = d1
  t = new IsohedralTiling(26);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]), p3(cp[3][0], cp[3][1]), p4(cp[4][0], cp[4][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p3), 2);
  ASSERT_EQ(orientation(p2, p3, p4), 2);
  ASSERT_EQ(orientation(p3, p4, p0), 2);
  ASSERT_EQ(orientation(p4, p0, p1), 2);

  ASSERT_FLOAT_EQ(p1.x, p2.x); 
  ASSERT_FLOAT_EQ(p0.x, p3.x);
  ASSERT_FLOAT_EQ((p0.y+p3.y)/2., p4.y);
  ASSERT_FLOAT_EQ((p1.y+p2.y)/2., p4.y);

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

TEST_F(NonFundamentalTest, IH67) {
  // IH 67 has S(P) = cmm, S(P|M) = d1
  t = new IsohedralTiling(67);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p0.y, p3.y);
  ASSERT_FLOAT_EQ(p1.y, p2.y);

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

TEST_F(NonFundamentalTest, IH91) {
  // IH 91 has S(P) = cmm, S(P|M) = d1
  t = new IsohedralTiling(91);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p0), 2);
  ASSERT_EQ(orientation(p2, p0, p1), 2);

  ASSERT_FLOAT_EQ(p0.y, p1.y);
  ASSERT_FLOAT_EQ((p1.x+p0.x)/2., p2.x);

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

TEST_F(NonFundamentalTest, IH16) {
  // IH 16 has S(P) = p31m, S(P|M) = d1
  t = new IsohedralTiling(16);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p4.x, p1.x); 
  ASSERT_FLOAT_EQ(p2.y, p0.y);
  ASSERT_FLOAT_EQ(p3.y, p5.y);

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

TEST_F(NonFundamentalTest, IH36) {
  // IH 36 has S(P) = p31m, S(P|M) = d1
  t = new IsohedralTiling(36);
  c = new Colloid();

  // Bad motif symmetry
  try {
    m.load("../motif_library/c1_random.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d1_vitruvian.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.y, p2.y); 
  ASSERT_FLOAT_EQ(p1.x, p3.x);

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

  params[2] = M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH29) {
  // IH 29 has S(P) = p4g, S(P|M) = d1
  t = new IsohedralTiling(29);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]), p3(cp[3][0], cp[3][1]), p4(cp[4][0], cp[4][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p3), 2);
  ASSERT_EQ(orientation(p2, p3, p4), 2);
  ASSERT_EQ(orientation(p3, p4, p0), 2);
  ASSERT_EQ(orientation(p4, p0, p1), 2);

  ASSERT_FLOAT_EQ(p3.y, p4.y); 
  ASSERT_FLOAT_EQ((p3.x+p4.x)/2., p1.x);
  ASSERT_FLOAT_EQ(p2.y, p0.y); 

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

  params[2] = M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH71) {
  // IH 71 has S(P) = p4g, S(P|M) = d1
  t = new IsohedralTiling(71);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p1.y, p2.y); 
  ASSERT_FLOAT_EQ(p0.y, p3.y);
  ASSERT_FLOAT_EQ(p1.x, p0.x);
  ASSERT_FLOAT_EQ(p3.x, p2.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // back to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // back to center of tile

  // Try to rotate
  params = c->getParameters();

  params[2] = 0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 7*M_PI/4.);

  params[2] = M_PI/4*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 7*M_PI/4.);

  params[2] = M_PI/4.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/4.);

  params[2] = 5*M_PI/4.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/4.);

  params[2] = 5*M_PI/4.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 7*M_PI/4.);

  params[2] = 7*M_PI/4.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 7*M_PI/4.);

  params[2] = 7*M_PI/4.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 7*M_PI/4.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 7*M_PI/4.);
}

TEST_F(NonFundamentalTest, IH82) {
  // IH 82 has S(P) = p4m, S(P|M) = d1
  t = new IsohedralTiling(82);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p0), 2);
  ASSERT_EQ(orientation(p2, p0, p1), 2);

  ASSERT_FLOAT_EQ(p2.y, p0.y); 
  ASSERT_FLOAT_EQ((p2.x+p0.x)/2., p1.x);

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

  params[2] = M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH32) {
  // IH 32 has S(P) = p6m, S(P|M) = d1
  t = new IsohedralTiling(32);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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

  ASSERT_FLOAT_EQ(p0.y, p2.y); 
  ASSERT_FLOAT_EQ(p1.x, p3.x);

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

TEST_F(NonFundamentalTest, IH40) {
  // IH 40 has S(P) = p6m, S(P|M) = d1
  t = new IsohedralTiling(40);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p0), 2);
  ASSERT_EQ(orientation(p2, p0, p1), 2);

  ASSERT_FLOAT_EQ(p2.y, p0.y); 
  ASSERT_FLOAT_EQ((p2.x+p0.x)/2., p1.x);

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

  params[2] = M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], M_PI/2.);

  params[2] = M_PI*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*0.99;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);

  params[2] = 3*M_PI/2.*1.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
 
  params[2] = -0.01;
  c->setParameters(params);
  ASSERT_FLOAT_EQ(c->getParameters()[2], 3*M_PI/2.);
}

TEST_F(NonFundamentalTest, IH72) {
  // IH 72 has S(P) = pmm, S(P|M) = d2
  t = new IsohedralTiling(72);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d2_dumbbell.json");
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

  ASSERT_FLOAT_EQ(p1.y, p0.y);
  ASSERT_FLOAT_EQ(p2.y, p3.y);
  ASSERT_FLOAT_EQ(p2.x, p1.x);
  ASSERT_FLOAT_EQ(p3.x, p0.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 2;
  const double delta = 0.01;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta));
    }
  }
}

TEST_F(NonFundamentalTest, IH17) {
  // IH 17 has S(P) = cmm, S(P|M) = d2
  t = new IsohedralTiling(17);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d2_dumbbell.json");
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

  ASSERT_FLOAT_EQ(p1.y, p2.y);
  ASSERT_FLOAT_EQ(p4.y, p5.y);
  ASSERT_FLOAT_EQ(p0.y, p3.y);
  ASSERT_FLOAT_EQ(p5.x, p1.x);
  ASSERT_FLOAT_EQ(p2.x, p4.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 2;
  const double delta = 0.01;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta));
    }
  }
}

TEST_F(NonFundamentalTest, IH74) {
  // IH 74 has S(P) = cmm, S(P|M) = d2
  t = new IsohedralTiling(74);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d2_dumbbell.json");
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

  ASSERT_FLOAT_EQ(p0.y, p2.y);
  ASSERT_FLOAT_EQ(p3.x, p1.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 2;
  const double delta = 0.01;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta));
    }
  }
}

TEST_F(NonFundamentalTest, IH73) {
  // IH 73 has S(P) = p4g, S(P|M) = d2
  t = new IsohedralTiling(73);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.x, p1.x);
  ASSERT_FLOAT_EQ(p0.y, p3.y);
  ASSERT_FLOAT_EQ(p1.y, p2.y);
  ASSERT_FLOAT_EQ(p2.x, p3.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 2;
  const double delta = 0.01;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta));
    }
  }
}

TEST_F(NonFundamentalTest, IH37) {
  // IH 37 has S(P) = p6m, S(P|M) = d2
  t = new IsohedralTiling(37);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.y, p2.y);
  ASSERT_FLOAT_EQ(p1.x, p3.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 2;
  const double delta = 0.01;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta));

      params[2] = i*dtheta + delta;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta));
    }
  }
}

TEST_F(NonFundamentalTest, IH18) {
  // IH 18 has S(P) = p31m, S(P|M) = d3
  t = new IsohedralTiling(18);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d3_triangle.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.y, p5.y);
  ASSERT_FLOAT_EQ(p1.y, p4.y);
  ASSERT_FLOAT_EQ(p2.y, p3.y);
  ASSERT_FLOAT_EQ(p0.x, p0.x);
  ASSERT_FLOAT_EQ(p3.x, p5.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 3;
  const double delta = 0.01, offset = M_PI/6.0;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta + offset));

      params[2] = i*dtheta + delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta + offset));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta + offset));

      params[2] = i*dtheta + delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta + offset));
    }
  }
}

TEST_F(NonFundamentalTest, IH93) {
  // IH 93 has S(P) = p6m, S(P|M) = d3
  t = new IsohedralTiling(93);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d3_triangle.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p0), 2);
  ASSERT_EQ(orientation(p2, p0, p1), 2);

  ASSERT_FLOAT_EQ(p0.x, p2.x);
  ASSERT_FLOAT_EQ(p1.y, (p0.y+p2.y)/2.);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 1/3.); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 3;
  const double delta = 0.01, offset = 0.0;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta + offset));

      params[2] = i*dtheta + delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta + offset));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta + offset));

      params[2] = i*dtheta + delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta + offset));
    }
  }
}

TEST_F(NonFundamentalTest, IH76) {
  // IH 76 has S(P) = p4m, S(P|M) = d4
  t = new IsohedralTiling(76);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d3_triangle.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/d4_square.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.x, p1.x);
  ASSERT_FLOAT_EQ(p2.x, p3.x);
  ASSERT_FLOAT_EQ(p0.y, p3.y);
  ASSERT_FLOAT_EQ(p1.y, p2.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate
  params = c->getParameters();

  const int n_mirrors = 4;
  const double delta = 0.01, offset = 0.0;
  double dtheta = (M_PI/n_mirrors/2.);
  for (unsigned int i = 0; i <= 4*n_mirrors; ++i) {
    if (i%2 == 0) { // +/- delta around allowed angle  
      params[2] = i*dtheta - delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta + offset));

      params[2] = i*dtheta + delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta + offset));
    } else { // +/- delta around halfway point
      params[2] = i*dtheta - delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i-1)*dtheta + offset));

      params[2] = i*dtheta + delta + offset;
      c->setParameters(params);
      ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds((i+1)*dtheta + offset));
    }
  }
}

TEST_F(NonFundamentalTest, IH8) {
  // IH 8 has S(P) = p2, S(P|M) = c2
  t = new IsohedralTiling(8);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p4.x, p2.x);
  ASSERT_FLOAT_EQ(p5.x, p1.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH57) {
  // IH 57 has S(P) = p2, S(P|M) = c2
  t = new IsohedralTiling(57);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.x, p3.x);
  ASSERT_FLOAT_EQ(p1.x, p2.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH09) {
  // IH 9 has S(P) = pgg, S(P|M) = c2
  t = new IsohedralTiling(9);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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
  ASSERT_EQ(orientation(p3, p5, p0), 2);
  ASSERT_EQ(orientation(p5, p0, p1), 2);

  ASSERT_FLOAT_EQ(p1.x, p5.x);
  ASSERT_FLOAT_EQ(p4.x, p2.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH59) {
  // IH 59 has S(P) = pgg, S(P|M) = c2
  t = new IsohedralTiling(59);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.y, p2.y);
  ASSERT_FLOAT_EQ(p1.x, p3.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH58) {
  // IH 58 has S(P) = pmg, S(P|M) = c2
  t = new IsohedralTiling(58);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p1.x, p2.x);
  ASSERT_FLOAT_EQ(p0.x, p3.x);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH61) {
  // IH 61 has S(P) = p4, S(P|M) = c2
  t = new IsohedralTiling(61);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p1.x, p0.x);
  ASSERT_FLOAT_EQ(p2.x, p3.x);
  ASSERT_FLOAT_EQ(p1.y, p2.y);
  ASSERT_FLOAT_EQ(p0.y, p3.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH34) {
  // IH 34 has S(P) = p6, S(P|M) = c2
  t = new IsohedralTiling(34);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p1.x, p3.x);
  ASSERT_FLOAT_EQ(p0.y, p2.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH10) {
  // IH 10 has S(P) = p3, S(P|M) = c3
  t = new IsohedralTiling(10);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c3_swirl.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.x, p2.x);
  ASSERT_FLOAT_EQ(p3.x, p5.x);
  ASSERT_FLOAT_EQ(p1.y, p4.y);
  ASSERT_FLOAT_EQ(p0.y, p5.y);
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
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH90) {
  // IH 90 has S(P) = p6, S(P|M) = c3
  t = new IsohedralTiling(90);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c3_swirl.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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
  dvec2 p0(cp[0][0], cp[0][1]), p1(cp[1][0], cp[1][1]), p2(cp[2][0], cp[2][1]);
  ASSERT_EQ(orientation(p0, p1, p2), 2);
  ASSERT_EQ(orientation(p1, p2, p0), 2);
  ASSERT_EQ(orientation(p2, p0, p1), 2);

  ASSERT_FLOAT_EQ(p0.x, p2.x);
  ASSERT_FLOAT_EQ(p1.y, (p0.y+p2.y)/2.);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 1/3.); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH62) {
  // IH 62 has S(P) = p4, S(P|M) = c4
  t = new IsohedralTiling(62);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c3_swirl.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c4_swirl.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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

  ASSERT_FLOAT_EQ(p0.x, p1.x);
  ASSERT_FLOAT_EQ(p3.x, p2.x);
  ASSERT_FLOAT_EQ(p0.y, p3.y);
  ASSERT_FLOAT_EQ(p1.y, p2.y);

  // Check COM is put back on mirror line
  vector<double> params = c->getParameters(), new_params;
  params[0] = 1; // Try to place motif in top right corner
  params[1] = 1;
  c->setParameters(params);
  new_params = c->getParameters();
  for (unsigned int i=2; i < params.size(); ++i) {
    ASSERT_FLOAT_EQ(params[i], new_params[i]);
  }
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

TEST_F(NonFundamentalTest, IH11) {
  // IH 11 has S(P) = p6, S(P|M) = c6
  t = new IsohedralTiling(11);
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
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/d2_dumbbell.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c2_taiji.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c3_swirl.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Bad motif symmetry
  try {
    m.load("../motif_library/c4_swirl.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.25);
    c->setU0(u0);
    c->setDform(df);
    c->setDU(0.1);
    c->init(true);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  // Good initialization
  try {
    m.load("../motif_library/c6_swirl_L.json");
    c->setMotif(m);
    c->setTile(*t);
    vector<double> u0(t->numEdgeShapes(), 0.1);
    vector<double> df(t->numEdgeShapes(), 0.1);
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
  ASSERT_EQ(orientation(p3, p5, p0), 2);
  ASSERT_EQ(orientation(p5, p0, p1), 2);

  ASSERT_FLOAT_EQ(p0.x, p2.x);
  ASSERT_FLOAT_EQ(p3.x, p5.x);
  ASSERT_FLOAT_EQ(p0.y, p5.y);
  ASSERT_FLOAT_EQ(p1.y, p4.y);
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
  ASSERT_FLOAT_EQ(new_params[0], 0.5); // to center of tile
  ASSERT_FLOAT_EQ(new_params[1], 0.5); // to center of tile

  // Try to rotate - should have no affect
  params = c->getParameters();

  const int n = 10;
  double dtheta = (2*M_PI/n);
  for (unsigned int i = 0; i <= n; ++i) {
    params[2] = i*dtheta;
    c->setParameters(params);
    ASSERT_FLOAT_EQ(c->getParameters()[2], thetaBounds(i*dtheta));
  }
}

#endif  // TEST_NONFDTILES_HPP_
