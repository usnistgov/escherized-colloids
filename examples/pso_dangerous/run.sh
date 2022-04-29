#!/bin/bash
# Compile the example.
#make clean; make;

# 1. Run the code; 2. Examine unit cell to check the answer.

# S(P|M) = d1 using S(M) = d2 Examples
# ------------------------------------
# Dangerous (d2 not explicitly forbidden)
#./pso -n 13 -m ../../motif_library/d2_dumbbell.json -l d1_d2_IH13_p2mg_dangerous; # works - warned by watch_out_for()
#./pso -n 16 -m ../../motif_library/d2_dumbbell.json -l d1_d2_IH16_p31m_dangerous; # subjective - warned by watch_out_for()

# Forbidden (d2 is forbidden supergroup)
#./pso -n 64 -m ../../motif_library/d2_dumbbell.json -l d1_d2_IH64_p1m1_forbidden;
#./pso -n 12 -m ../../motif_library/d2_dumbbell.json -l d1_d2_IH12_c1m1_forbidden;

# S(P|M) = d3/c3 using S(M) = d6 Examples
# ---------------------------------------
# Dangerous (d6 not explicitly forbidden)
#./pso -n 93 -m ../../motif_library/d6_hexagon.json -l d3_d6_IH93_p6mm_dangerous; # works - no warnings
#./pso -n 90 -m ../../motif_library/d6_hexagon.json -l c3_d6_IH90_p6_dangerous; # works - warned by watch_out_for()

# Forbidden (d6 is a forbidden supergroup)
#./pso -n 18 -m ../../motif_library/d6_hexagon.json -l d3_d6_IH18_p31m_forbidden;
#./pso -n 12 -m ../../motif_library/d6_hexagon.json -l d1_d6_IH12_c1m1_forbidden;

# S(P|M) = d6/d2 using S(M) = dinf Examples
# --------------------------------------
# Dangerous (d6 not explicitly forbidden)
#./pso -n 20 -m ../../motif_library/dinf_circle.json -l d6_dinf_IH20_p6mm_dangerous; # works - no warnings
#./pso -n 72 -m ../../motif_library/dinf_circle.json -l d2_dinf_IH72_p2mm_dangerous; # doesn't - correctly warned about p4mm in watch_out_for()

# Forbidden (dinf is a forbidden supergroup)
#./pso -n 1 -m ../../motif_library/dinf_circle.json -l c1_dinf_IH01_p1_forbidden;
#./pso -n 10 -m ../../motif_library/dinf_circle.json -l c3_dinf_IH10_p3_forbidden;

# S(P|M) = c1/2/3 using S(M) = c6 Examples
# ------------------------------------
# Dangerous (c6 not explicitly forbidden)
#./pso -n 61 -m ../../motif_library/c6_swirl_D.json -l c2_c6_IH61_p4_dangerous; # works - no warnings
#./pso -n 58 -m ../../motif_library/c6_swirl_D.json -l c2_c6_IH58_p2mg_dangerous; # works - no warnings

# Forbidden (c6 is a forbidden supergroup)
#./pso -n 41 -m ../../motif_library/c6_swirl_D.json -l c1_c6_IH41_p1_forbidden;
#./pso -n 10 -m ../../motif_library/c6_swirl_D.json -l c3_c6_IH10_p3_forbidden;

