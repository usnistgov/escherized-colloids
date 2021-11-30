#++ -Wall -std=c++11 -O3 -march=native -ffp-contract=fast -I/usr/include/ -I/home/nam/NIST/workspace/magic_colloids/escherized_colloids/optim/header_only_version pso.cpp -o pso.out

INCLUDES = -I../../tactile/deps/glm \
	-I../../tactile/deps/gl3w/build/include \
	-I../../tactile/deps/glfw-3.2.1/include \
	-I../../src

CC = gcc 
CPP = g++
CFLAGS = $(INCLUDES) -I../../tactile/src -g -Wall 
CPPFLAGS = $(INCLUDES) -I../../tactile/src -I../../src/ -g -Wall -std=c++11 -O3
TACTILE_CPPFLAGS = -I../../tactile/deps/glm -I../../tactile/src -g -Wall -std=c++11 -O3

LIBS = -L../../tactile/deps/glfw-3.2.1/build/src \
	-framework OpenGL -framework IOKit -framework CoreVideo

pso: pso.o colloid.o tiling.o motif.o
	$(CPP) -o $@ $^ 

gl3w.o: ../../tactile/deps/gl3w/build/src/gl3w.c
	$(CC) $(CFLAGS) $^ -c -o $@

tiling.o: ../../tactile/src/tiling.cpp ../../tactile/src/tiling_arraydecl.inc
	$(CPP) $(TACTILE_CPPFLAGS) ../../tactile/src/tiling.cpp -c -o $@

motif.o: ../../src/motif.cpp
	$(CPP) $(CPPFLAGS) ../../src/motif.cpp -c -o $@

colloid.o: ../../src/colloid.cpp ../../tactile/src/tiling.cpp
	$(CPP) $(CPPFLAGS) ../../src/colloid.cpp -c -o $@

pso.o: pso.cpp 
	$(CPP) $(CPPFLAGS) -march=native -ffp-contract=fast -I/usr/include/ -I../../optim/header_only_version $^ -c -o $@

clean:
	rm -f *.o pso
