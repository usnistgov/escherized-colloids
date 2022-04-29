#!/bin/bash
# Compile the example.
make clean; make;

du=0.1 # Default spacing

# "Safe" tiling options
for tile in 64 12 14 68 13 15 66 69 26 67 91 16 36 29 71 82 32 40; do
	if [ $tile -eq 40 ]; then
		du=0.05;
	elif [ $tile -eq 69]; then
		du=0.05;
	elif [ $tile -eq 13]; then
		du=0.15
	fi	

	# Using interaction between tile and motif
	./pso -n $tile -m ../../motif_library/d1_vitruvian.json -l energized_$tile -u $du -e 1.0;
	cat energized_"$tile"_*.xyz >> energized_"$tile".xyz; # Summarize the results

	# No interaction between tile and motif
	./pso -n $tile -m ../../motif_library/d1_vitruvian.json -l deenergized_$tile -u $du -e 0.0;
        cat deenergized_"$tile"_*.xyz >> deenergized_"$tile".xyz; # Summarize the results

	# To easily compare the two versions of the colloid
	cat *energized_"$tile"_success.xyz >> compare_"$tile".xyz;
done;
