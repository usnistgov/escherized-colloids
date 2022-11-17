#!/bin/bash
# Compile the example.
make clean; make;

# Run the code.
# S(P|M) = d1 Examples
./pso -n 64 -m ../../motif_library/d1_vitruvian.json -l tile_64_vm;
#./pso -n 12 -m ../../motif_library/d1_vitruvian.json -l tile_12_vm;
#./pso -n 14 -m ../../motif_library/d1_vitruvian.json -l tile_14_vm;
#./pso -n 68 -m ../../motif_library/d1_vitruvian.json -l tile_68_vm;
#./pso -n 13 -m ../../motif_library/d1_vitruvian.json -l tile_13_vm;
#./pso -n 15 -m ../../motif_library/d1_vitruvian.json -l tile_15_vm;
#./pso -n 66 -m ../../motif_library/d1_vitruvian.json -l tile_66_vm;
#./pso -n 69 -m ../../motif_library/d1_vitruvian.json -l tile_69_vm;
#./pso -n 26 -m ../../motif_library/d1_vitruvian.json -l tile_26_vm;
#./pso -n 67 -m ../../motif_library/d1_vitruvian.json -l tile_67_vm;
#./pso -n 91 -m ../../motif_library/d1_vitruvian.json -l tile_91_vm;
#./pso -n 16 -m ../../motif_library/d1_vitruvian.json -l tile_16_vm;
#./pso -n 36 -m ../../motif_library/d1_vitruvian.json -l tile_36_vm;
#./pso -n 29 -m ../../motif_library/d1_vitruvian.json -l tile_29_vm;
#./pso -n 71 -m ../../motif_library/d1_vitruvian.json -l tile_71_vm;
#./pso -n 82 -m ../../motif_library/d1_vitruvian.json -l tile_82_vm;
#./pso -n 32 -m ../../motif_library/d1_vitruvian.json -l tile_32_vm;
#./pso -n 40 -m ../../motif_library/d1_vitruvian.json -l tile_40_vm;

# S(P|M) = d2 Examples
#./pso -n 72 -m ../../motif_library/d2_dumbbell.json -l tile_72_db;
#./pso -n 17 -m ../../motif_library/d2_dumbbell.json -l tile_17_db;
#./pso -n 74 -m ../../motif_library/d2_dumbbell.json -l tile_74_db;
#./pso -n 73 -m ../../motif_library/d2_dumbbell.json -l tile_73_db;
#./pso -n 37 -m ../../motif_library/d2_dumbbell.json -l tile_37_db;

# S(P|M) = d3 Examples
#./pso -n 18 -m ../../motif_library/d3_triangle.json -l tile_18_tr;
#./pso -n 93 -m ../../motif_library/d3_triangle.json -l tile_93_tr;

# S(P|M) = d4 Examples
#./pso -n 76 -m ../../motif_library/d4_square.json -l tile_76_sq;

# S(P|M) = d6 Examples
#./pso -n 20 -m ../../motif_library/d6_hexagon.json -l tile_20_hx;

# S(P|M) = c2 Examples
#./pso -n 8 -m ../../motif_library/c2_taiji.json -l tile_08_tj;
#./pso -n 57 -m ../../motif_library/c2_taiji.json -l tile_57_tj;
#./pso -n 9 -m ../../motif_library/c2_taiji.json -l tile_09_tj;
#./pso -n 59 -m ../../motif_library/c2_taiji.json -l tile_59_tj;
#./pso -n 58 -m ../../motif_library/c2_taiji.json -l tile_58_tj;
#./pso -n 61 -m ../../motif_library/c2_taiji.json -l tile_61_tj;
#./pso -n 34 -m ../../motif_library/c2_taiji.json -l tile_34_tj;

# S(P|M) = c3 Examples
#./pso -n 10 -m ../../motif_library/c3_swirl.json -l tile_10_sw3;
#./pso -n 90 -m ../../motif_library/c3_swirl.json -l tile_90_sw3;

# S(P|M) = c4 Examples
#./pso -n 62 -m ../../motif_library/c4_swirl.json -l tile_62_sw4;

# S(P|M) = c6 Examples
#./pso -n 11 -m ../../motif_library/c6_swirl_L.json -l tile_11_sw6;

# Visualize the output (success.xyz) with any standard molecular visualization program.

