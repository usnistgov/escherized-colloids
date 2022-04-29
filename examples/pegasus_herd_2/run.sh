#!/bin/bash
# Compile the example.
make clean; make;

du=0.1 # Default spacing

# "Safe" tiling options
for tile in 1 41 2 3 43 44 42 22 45 83 4 23 46 47 84 5 6 25 27 51 52 53 86 24 49 50 85 54 78 7 33 30 38 28 55 79 56 81 21 31 39 88 77; do
	#if [ $tile -eq 40 ]; then
	#	du=0.05;
	#elif [ $tile -eq 69]; then
	#	du=0.05;
	#elif [ $tile -eq 13]; then
	#	du=0.15
	#fi	

	# Using interaction between tile and motif
	./pso -n $tile -m ../../motif_library/c1_pegasus_2.json -l energized_$tile -u $du -e 0.1;
	cat energized_"$tile"_*.xyz >> energized_"$tile".xyz; # Summarize the results

	# No interaction between tile and motif
	./pso -n $tile -m ../../motif_library/c1_pegasus_2.json -l deenergized_$tile -u $du -e 0.0;
        cat deenergized_"$tile"_*.xyz >> deenergized_"$tile".xyz; # Summarize the results

	# To easily compare the two versions of the colloid
	cat *energized_"$tile"_success.xyz >> compare_"$tile".xyz;
done;
