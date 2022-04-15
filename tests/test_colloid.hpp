#ifndef TEST_COLLOID_HPP_
#define TEST_COLLOID_HPP_

#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "../src/colloid.hpp"

class ColloidTest : public ::testing::Test {
 protected:
  void SetUp() override {
     m.load("../motif_library/d1_vitruvian.json");
     c = new Colloid();
     t = new IsohedralTiling(7);
     c->setMotif(m);
     c->setTile(*t);
     vector<double> u0(t->numEdgeShapes(), 0.1);
     vector<double> df(t->numEdgeShapes(), 0.25);
     c->setU0(u0);
     c->setDform(df);
     //c->setDU(0.1);
     //c->init();
  }
  Motif m;
  Colloid* c;
  IsohedralTiling* t;
};

TEST_F(ColloidTest, BuildTilePolygon) {
  vector<double> df(t->numEdgeShapes(), 0.0); // Make straight edges
  c->setDform(df);

  const int N = 10;
  vector<vector<dvec2>> polygon = c->buildTilePolygon(N);

  for (unsigned int i=0; i < polygon.size(); ++i) {
    ASSERT_EQ(polygon[i].size(), N); // Check each edge has N points
    double d2 = 0.0;
    for (unsigned int j=0; j < polygon[i].size()-1; ++j) { // Check points are evenly spaced when df=0
      double x2 = pow(polygon[i][j].x - polygon[i][j+1].x, 2) + pow(polygon[i][j].y - polygon[i][j+1].y, 2);
      if (j == 0) {
        d2 = x2;
      } else {
        EXPECT_FLOAT_EQ(x2, d2);
      }
    }
  }
}

TEST_F(ColloidTest, MotifInside) {
  c->setTileScale(1.0);

  const int N = 10;
  vector<vector<dvec2>> polygon = c->buildTilePolygon(N);

  ASSERT_EQ(isInside(polygon, {-0.12, 0.32}), false);
  ASSERT_EQ(isInside(polygon, {-0.09, 0.32}), true);
  ASSERT_EQ(isInside(polygon, {0.5, 0.32}), true);
  ASSERT_EQ(isInside(polygon, {0.88, 0.32}), true);
  ASSERT_EQ(isInside(polygon, {0.92, 0.32}), false);

  ASSERT_EQ(isInside(polygon, {0.5, -0.3}), false);
  ASSERT_EQ(isInside(polygon, {0.5, -0.28}), true);
  ASSERT_EQ(isInside(polygon, {0.5, 0.5}), true);
  ASSERT_EQ(isInside(polygon, {0.5, 0.84}), true);
  ASSERT_EQ(isInside(polygon, {0.5, 0.88}), false);
}

#endif  // TEST_COLLOID_HPP_

