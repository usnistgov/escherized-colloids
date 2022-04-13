#!/bin/bash
# Compile the example.
make clean; make;

# Run
# If you do not specify any pX values, defaults will be used; just specify tile number, n.
./init -n 7 -p0 0.6 -p1 0.19 -m ../../motif_library/d1_vitruvian.json;

# You can visualize the resulting colloid.xyz using any standard molecular simulation visualization tool.
