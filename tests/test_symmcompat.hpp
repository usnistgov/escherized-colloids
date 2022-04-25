#ifndef TEST_SYMMCOMPAT_HPP_
#define TEST_SYMMCOMPAT_HPP_

#include <iostream>
#include <string>
#include "../src/motif.hpp"

using namespace std;

class SymmetryCompatibility : public ::testing::Test {
 protected:
  void SetUp() override {
    n = 0;
    c1.load("../motif_library/c1_random.json");
    c2.load("../motif_library/c2_taiji.json");
    c3.load("../motif_library/c3_swirl.json");
    c4.load("../motif_library/c4_swirl.json");
    c6L.load("../motif_library/c6_swirl_L.json");
    c6D.load("../motif_library/c6_swirl_D.json");

    d1.load("../motif_library/d1_vitruvian.json");
    d2.load("../motif_library/d2_dumbbell.json");
    d3.load("../motif_library/d3_triangle.json");
    d4.load("../motif_library/d4_square.json");
    d6.load("../motif_library/d6_hexagon.json");
    dinf.load("../motif_library/dinf_circle.json");
  }

  Motif c1, c2, c3, c4, c6L, c6D, d1, d2, d3, d4, d6, dinf;
  int n;
};

TEST_F(SymmetryCompatibility, ReflectionSame) {
  // Check d(n)
  try {
    n = compatibility_check(d1, "d", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 1);

  // Check d(n)
  try {
    n = compatibility_check(d2, "d", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 2);

  // Check d(n)
  try {
    n = compatibility_check(d3, "d", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 3);

  // Check d(n)
  try {
    n = compatibility_check(d4, "d", 4);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 4);

  // Check d(n)
  try {
    n = compatibility_check(d6, "d", 6);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 6);

  // Check d(n)
  try {
    n = compatibility_check(dinf, "d", 6);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, INF_SYMM);
}

TEST_F(SymmetryCompatibility, InvalidReflectionSubgroup) {
  // Motif doesn't have enough mirrors
  try {
    n = compatibility_check(d1, "d", 2);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(d4, "d", 6);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(d3, "d", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }
}

TEST_F(SymmetryCompatibility, InvalidReflectionSupergroup) {
  // Number of mirrors is incompatible
  try {
    n = compatibility_check(d4, "d", 3);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(d6, "d", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(d6, "d", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }
}

TEST_F(SymmetryCompatibility, ValidReflectionSupergroup) {
  // Number of mirrors IS compatible
  try {
    n = compatibility_check(d6, "d", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d6, "d", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d4, "d", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d3, "d", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d4, "d", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  for (int i = 1; i <= 6; ++i) {
    try {
      n = compatibility_check(dinf, "d", i);
    } catch (const customException &e) {
      std::cerr << e.getMessage() << std::endl;
      ASSERT_TRUE(false);
    }
  }
}

TEST_F(SymmetryCompatibility, RotationSame) {
  // Check c(n)
  try {
    n = compatibility_check(c1, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 1);

  // Check c(n)
  try {
    n = compatibility_check(c2, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 2);

  // Check c(n)
  try {
    n = compatibility_check(c3, "c", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 3);

  // Check c(n)
  try {
    n = compatibility_check(c4, "c", 4);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 4);

  // Check c(n)
  try {
    n = compatibility_check(c6L, "c", 6);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 6);

  // Check c(n)
  try {
    n = compatibility_check(c6D, "c", 6);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  ASSERT_EQ(n, 6);
}

TEST_F(SymmetryCompatibility, InvalidRotationSubgroup) {
  // Motif doesn't have enough rotation
  try {
    n = compatibility_check(c1, "c", 2);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(c4, "c", 6);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(c3, "c", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }
}

TEST_F(SymmetryCompatibility, InvalidRotationSupergroup) {
  // Number of rotations is incompatible
  try {
    n = compatibility_check(c4, "c", 3);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(c6L, "c", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(c6D, "c", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }
}

TEST_F(SymmetryCompatibility, ValidRotationSupergroup) {
  // Number of rotations IS compatible
  try {
    n = compatibility_check(c6L, "c", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(c6D, "c", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(c6L, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(c6D, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(c4, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(c3, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(c4, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
}

TEST_F(SymmetryCompatibility, InvalidCross) {
  // c(n) can never successfully induce any d(n)
  for (int i = 1; i < 6; ++i) {
    try {
      n = compatibility_check(c1, "d", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 1; i < 6; ++i) {
    try {
      n = compatibility_check(c2, "d", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 1; i < 6; ++i) {
    try {
      n = compatibility_check(c3, "d", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 1; i < 6; ++i) {
    try {
      n = compatibility_check(c4, "d", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 1; i < 6; ++i) {
    try {
      n = compatibility_check(c6L, "d", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 1; i < 6; ++i) {
    try {
      n = compatibility_check(c6D, "d", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  // d(n) can induce c(n) but not more than n or with incompatible lower n
  for (int i = 2; i < 6; ++i) {
    try {
      n = compatibility_check(d1, "c", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 3; i < 6; ++i) {
    try {
      n = compatibility_check(d2, "c", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }

  for (int i = 4; i < 6; ++i) {
    try {
      n = compatibility_check(d3, "c", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }
  try {
    n = compatibility_check(d3, "c", 2);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  for (int i = 5; i < 6; ++i) {
    try {
      n = compatibility_check(d4, "c", i);
    } catch (const customException &e) {
      ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
    }
  }
  try {
    n = compatibility_check(d4, "c", 3);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }
  try {
    n = compatibility_check(d4, "c", 6);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }

  try {
    n = compatibility_check(d6, "c", 4);
  } catch (const customException &e) {
    ASSERT_EQ(e.getMessage().compare("motif's symmetry is incompatible with the tile"), 0);
  }
}

TEST_F(SymmetryCompatibility, ValidCross) {
  // d(n) can induce c(n) with same or smaller n as long as compatible
  for (int i = 1; i <= 6; ++i) {
    try {
      n = compatibility_check(dinf, "d", i);
    } catch (const customException &e) {
      std::cerr << e.getMessage() << std::endl;
      ASSERT_TRUE(false);
    }
  }

  try {
    n = compatibility_check(d6, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d6, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d6, "c", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d6, "c", 6);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d4, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d4, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d4, "c", 4);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d3, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d3, "c", 3);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d2, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
  try {
    n = compatibility_check(d2, "c", 2);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }

  try {
    n = compatibility_check(d1, "c", 1);
  } catch (const customException &e) {
    std::cerr << e.getMessage() << std::endl;
    ASSERT_TRUE(false);
  }
}

#endif  // TEST_SYMMCOMPAT_HPP_
