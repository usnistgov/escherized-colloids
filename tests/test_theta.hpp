#ifndef TEST_THETA_HPP_
#define TEST_THETA_HPP_

#include <vector>
#include <string>
#include <gtest/gtest.h>
#include "../src/utils.hpp"

class ThetaBoundsTest : public ::testing::Test {
 protected:
  void SetUp() override {
     theta = 0.0;
  }
  double theta;
};

TEST_F(ThetaBoundsTest, UnderOne) {
  theta = -0.1;
  EXPECT_FLOAT_EQ(thetaBounds(theta), 2.0*M_PI-0.1);
}

TEST_F(ThetaBoundsTest, UnderMany) {
  theta = -16.0*M_PI - 0.1;
  EXPECT_FLOAT_EQ(thetaBounds(theta), 2.0*M_PI-0.1);
}

TEST_F(ThetaBoundsTest, OverOne) {
  theta = 2.0*M_PI + 0.1;
  EXPECT_FLOAT_EQ(thetaBounds(theta), 0.1);
}

TEST_F(ThetaBoundsTest, OverMany) {
  theta = 16.0*M_PI + 0.1;
  EXPECT_FLOAT_EQ(thetaBounds(theta), 0.1);
}

TEST_F(ThetaBoundsTest, Inside) {
  theta = 0.1234;
  EXPECT_FLOAT_EQ(thetaBounds(theta), 0.1234);
}

#endif  // TEST_THETA_HPP_
