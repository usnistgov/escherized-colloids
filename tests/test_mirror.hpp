#ifndef TEST_MIRROR_HPP_
#define TEST_MIRROR_HPP_

#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "../src/utils.hpp"

using std::vector;
using std::string;

class MirrorAlignTest : public ::testing::Test {
 protected:
  void SetUp() override {
     p0.push_back(0.0);
     p0.push_back(0.0);
     p1.push_back(0.0);
     p1.push_back(0.0);
  }
  std::vector<double> p0, p1;
};

TEST_F(MirrorAlignTest, Q1) {
  p1[0] = 0.5;
  p1[1] = 0.5;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), M_PI/4.0);
}

TEST_F(MirrorAlignTest, Q12) {
  p1[0] = 0.0;
  p1[1] = 0.5;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), M_PI/2.0);
}

TEST_F(MirrorAlignTest, Q2) {
  p1[0] = -0.5;
  p1[1] = 0.5;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), 3.0*M_PI/4.0);
}

TEST_F(MirrorAlignTest, Q23) {
  p1[0] = -1.0;
  p1[1] = 0.0;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), M_PI);
}

TEST_F(MirrorAlignTest, Q3) {
  p1[0] = -0.5;
  p1[1] = -0.5;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), 5.0*M_PI/4.0);
}

TEST_F(MirrorAlignTest, Q34) {
  p1[0] = 0.0;
  p1[1] = -0.5;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), 6.0*M_PI/4.0);
}

TEST_F(MirrorAlignTest, Q4) {
  p1[0] = 0.5;
  p1[1] = -0.5;
  EXPECT_FLOAT_EQ(tile_mirror_alignment(p0, p1), 7.0*M_PI/4.0);
}

#endif  // TEST_MIRROR_HPP_
