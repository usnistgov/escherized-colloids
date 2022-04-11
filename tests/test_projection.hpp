#ifndef TEST_PROJECTION_HPP_
#define TEST_PROJECTION_HPP_

#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "../src/utils.hpp"

using std::vector;
using std::string;

class ProjectToLineTest : public ::testing::Test {
 protected:
  void SetUp() override {
     p0.push_back(1.1);
     p0.push_back(1.1);
     p1.push_back(1.1);
     p1.push_back(1.1);
     coords.push_back(1.234);
     coords.push_back(1.234);
  }
  std::vector<double> p0, p1, coords, result;
};

TEST_F(ProjectToLineTest, Identical) {
  // Check that if identical points throw an exception
  try {
    result = project_to_line(p0, p1, coords);
  } catch (customException &e) {
    EXPECT_EQ(e.getMessage(), string("points are identical"));
    return;
  } catch (...) {
    FAIL() << "expected a customException";
  }
  FAIL() << "missed exception.";
}

TEST_F(ProjectToLineTest, Vertical) {
  // Check projection to vertical line
  p1[1] = 2.1;
  try {
    result = project_to_line(p0, p1, coords);
  } catch (...) {
    FAIL() << "unexpected exception";
  }
  EXPECT_EQ(result[0], p0[0]); // x same as line
  EXPECT_EQ(result[1], coords[1]); // y unaffected
}

TEST_F(ProjectToLineTest, Horizontal) {
  // Check projection to horizontal line
  p1[0] = 2.1;
  try {
    result = project_to_line(p0, p1, coords);
  } catch (...) {
    FAIL() << "unexpected exception";
  }
  EXPECT_EQ(result[1], p0[1]); // y same as line
  EXPECT_EQ(result[0], coords[0]); // x unaffected
}

TEST_F(ProjectToLineTest, PositiveSlope) {
  // Line has positive slope
  p1[0] = 2.1;
  p1[1] = 2.1;

  // Coord on the "right"
  coords[0] = 2.2;
  coords[1] = 0.0;
  try {
    result = project_to_line(p0, p1, coords);
  } catch (...) {
    FAIL() << "unexpected exception";
  }  
  EXPECT_EQ(result[0], 1.1);
  EXPECT_EQ(result[1], 1.1); 

  // Coord on the "left"
  coords[1] = 2.2;
  coords[0] = 0.0;
  try {
    result = project_to_line(p0, p1, coords);
  } catch (...) {
    FAIL() << "unexpected exception";
  }  
  EXPECT_EQ(result[0], 1.1);
  EXPECT_EQ(result[1], 1.1); 
}

TEST_F(ProjectToLineTest, NegativeSlope) {
  // Line has positive slope
  p1[0] = 2.2;
  p1[1] = 0.0;

  // Coord on the "right"
  coords[0] = 2.2;
  coords[1] = 2.2;
  try {
    result = project_to_line(p0, p1, coords);
  } catch (...) {
    FAIL() << "unexpected exception";
  }  
  EXPECT_EQ(result[0], 1.1);
  EXPECT_EQ(result[1], 1.1); 

  // Coord on the "left"
  coords[1] = 0.0;
  coords[0] = 0.0;
  try {
    result = project_to_line(p0, p1, coords);
  } catch (...) {
    FAIL() << "unexpected exception";
  }  
  EXPECT_EQ(result[0], 1.1);
  EXPECT_EQ(result[1], 1.1); 
}

#endif  // TEST_PROJECTION_HPP_
