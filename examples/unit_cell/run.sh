#!/bin/bash
# Compile the example
make;

# Run
./unit_cell -n 7 -p0 0.6 -p1 0.19 -m ../../motif_library/d_1_vitruvian.json

# You can visualize the resulting colloid.xyz and unit_cell.xyz using any standard molecular simulation visualization tool
