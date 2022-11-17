#!/bin/bash
# Compile the example.
make clean; make;

# Vitruvian Man Examples
./unit_cell -n 7 -p0 0.6 -p1 0.19 -m ../../motif_library/d1_vitruvian.json;
mv colloid.xyz colloid_7_vm.xyz; 
mv unit_cell.xyz unit_cell_7_vm.xyz;

./unit_cell -n 28 -p0 0.4528 -p1 0.5 -m ../../motif_library/d1_vitruvian.json;
mv colloid.xyz colloid_28_vm.xyz; 
mv unit_cell.xyz unit_cell_28_vm.xyz;

./unit_cell -n 21 -p0 0.1045 -p1 0.65 -m ../../motif_library/d1_vitruvian.json;
mv colloid.xyz colloid_21_vm.xyz; 
mv unit_cell.xyz unit_cell_21_vm.xyz;

./unit_cell -n 4 -m ../../motif_library/d1_vitruvian.json;
mv colloid.xyz colloid_4_vm.xyz; 
mv unit_cell.xyz unit_cell_4_vm.xyz;

# Random Patch Examples
./unit_cell -n 7 -p0 0.6 -p1 0.19 -m ../../motif_library/c1_random.json;
mv colloid.xyz colloid_7_rp.xyz; 
mv unit_cell.xyz unit_cell_7_rp.xyz;

./unit_cell -n 28 -p0 0.4528 -p1 0.5 -m ../../motif_library/c1_random.json;
mv colloid.xyz colloid_28_rp.xyz; 
mv unit_cell.xyz unit_cell_28_rp.xyz;

./unit_cell -n 21 -p0 0.1045 -p1 0.65 -m ../../motif_library/c1_random.json;
mv colloid.xyz colloid_21_rp.xyz; 
mv unit_cell.xyz unit_cell_21_rp.xyz;

./unit_cell -n 4 -m ../../motif_library/c1_random.json;
mv colloid.xyz colloid_4_rp.xyz; 
mv unit_cell.xyz unit_cell_4_rp.xyz;

# You can visualize the resulting colloid.xyz and unit_cell.xyz using any standard molecular simulation visualization tool.
