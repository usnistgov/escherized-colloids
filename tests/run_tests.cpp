#include <gtest/gtest.h>

#include "test_colloid.hpp"
#include "test_fdtiles.hpp"
#include "test_mirror.hpp"
#include "test_nonfdtiles.hpp"
#include "test_projection.hpp"
#include "test_symmcompat.hpp"
#include "test_theta.hpp"

int main(int argc, char **argv) {
    ::testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS();
}
