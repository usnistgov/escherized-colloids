#!/bin/bash
# Compile the example.
make clean; make;

# Run the code.
#./pso -n 64 -m ../../motif_library/d1_vitruvian.json;
#./pso -n 12 -m ../../motif_library/d1_vitruvian.json;
#./pso -n 14 -m ../../motif_library/d1_vitruvian.json;
./pso -n 68 -m ../../motif_library/d1_vitruvian.json;

# Visualize the output (success.xyz) with any standard molecular visualization program.
