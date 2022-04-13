#ifndef TEST_FDTILES_HPP_
#define TEST_FDTILES_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "tiling.hpp"
#include "../src/colloid.hpp"

using namespace std;
using namespace csk;

class FundamentalTileTest : public ::testing::Test {
 protected:
  void SetUp() override {
    prefix = "./fd_ref/";

    for (unsigned int j = 0; j < sizeof(tiling_types) / sizeof(tiling_types[0]); ++j) { // Supported types in Tactile
      int type = static_cast<int>(tiling_types[j]);
      for (unsigned int i = 0; i < sizeof(FD_TYPES) / sizeof(FD_TYPES[0]); ++i) {
        if (type == FD_TYPES[i]) {
          ih_type.push_back(type);
        }
      }
    }
  }

  vector<int> ih_type;
  string prefix;

  void check_default(const int idx) {
    // Load reference
    stringstream ss;
    ss << prefix << "colloid_" << ih_type[idx] << ".json";
    Colloid *reference = new Colloid();
    reference->load(ss.str());

    // Construct a new default tiling of the given type.
    IsohedralTiling t( ih_type[idx] );
    Motif m;
    stringstream tt;
    tt << prefix << "c1_random.json";
    m.load(tt.str());
    Colloid c(m, t, 0.25);

    ASSERT_EQ(reference->isTileFundamental(), true);
    ASSERT_EQ(c.isTileFundamental(), true);

    // Compare parameters
    vector<double> p1 = c.getParameters(), p2 = reference->getParameters();
    ASSERT_EQ(p1.size(), p2.size());
    for (unsigned int i=0; i < p1.size(); ++i) {
      EXPECT_FLOAT_EQ(p1[i], p2[i]);
    }

    // Check tile information
    vector<vector<double>> b1 = c.getBoundaryCoords(), b2 = reference->getBoundaryCoords();
    ASSERT_EQ(b1.size(), b2.size());
    for (unsigned int i=0; i < b1.size(); ++i) {
      ASSERT_EQ(b1[i].size(), b2[i].size());
      for (unsigned int j=0; j < b1[i].size(); ++j) {
        EXPECT_FLOAT_EQ(b1[i][j], b2[i][j]);
      }
    }

    b1 = c.getTileControlPoints();
    b2 = reference->getTileControlPoints();
    ASSERT_EQ(b1.size(), b2.size());
    for (unsigned int i=0; i < b1.size(); ++i) {
      ASSERT_EQ(b1[i].size(), b2[i].size());
      for (unsigned int j=0; j < b1[i].size(); ++j) {
        EXPECT_FLOAT_EQ(b1[i][j], b2[i][j]);
      }
    }

    vector<int> d1 = c.getBoundaryIds(), d2 = reference->getBoundaryIds();
    ASSERT_EQ(d1.size(), d2.size());
    for (unsigned int i=0; i < d1.size(); ++i) {
      EXPECT_EQ(d1[i], d2[i]);
    }

    // Check motif information
    b1 = c.getMotif().getCoords();
    b2 = reference->getMotif().getCoords();
    ASSERT_EQ(b1.size(), b2.size());
    for (unsigned int i=0; i < b1.size(); ++i) {
      ASSERT_EQ(b1[i].size(), b2[i].size());
      for (unsigned int j=0; j < b1[i].size(); ++j) {
        EXPECT_FLOAT_EQ(b1[i][j], b2[i][j]);
      }
    }

    vector<string> s1 = c.getMotif().getTypes(), s2 = reference->getMotif().getTypes();
    ASSERT_EQ(s1.size(), s2.size());
    for (unsigned int i=0; i < s1.size(); ++i) {
      EXPECT_EQ(s1[i].compare(s2[i]), 0);
    }

    EXPECT_EQ(c.getMotif().getSymmetry().compare(reference->getMotif().getSymmetry()), 0);
  }
};

TEST_F(FundamentalTileTest, History) {
  // Check default tiles are being generated the same as they were before.
  for (unsigned int i=0; i < ih_type.size(); ++i) {
    check_default(i);
  }
}

#endif  // TEST_FDTILES_HPP_
