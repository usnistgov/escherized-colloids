## This isn't a newfangled, cross-platform, all-purpose Makefile.  You'll 
## probably have to customize it to get it to build on other computers,
## especially non-Mac computers.

SIMP_INC = -I../tactile/deps/glm 

INCLUDES = $(SIMP_INC) \
	-I../tactile/deps/gl3w/build/include \
	-I../tactile/deps/glfw-3.2.1/include \
	-I../src

CC = gcc 
CPP = g++
CFLAGS = $(INCLUDES) -I../tactile/src -g -Wall 
CPPFLAGS = $(INCLUDES) -I../tactile/src -I../src/ -g -Wall -std=c++11 
SIMP_CPPFLAGS = $(SIMP_INC) -I../tactile/src -g -Wall -std=c++11

## This is probably the only platform dependency.  Edit as necessary.
LIBS = -L../tactile/deps/glfw-3.2.1/build/src \
	-framework OpenGL -framework IOKit -framework CoreVideo

main: main.o colloid.o tiling.o motif.o
	$(CPP) -o $@ $^ 

## gl3w is even simpler
gl3w.o: ../tactile/deps/gl3w/build/src/gl3w.c
	$(CC) $(CFLAGS) $^ -c -o $@

## And my library is inspired by the simple ones above.  Throw it in!
## This compilation unit doesn't need all the other header files, just
## GLM.
tiling.o: ../tactile/src/tiling.cpp ../tactile/src/tiling_arraydecl.inc
	$(CPP) $(SIMP_CPPFLAGS) ../tactile/src/tiling.cpp -c -o $@

motif.o: ../src/motif.cpp
	$(CPP) $(CPPFLAGS) ../src/motif.cpp -c -o $@

colloid.o: ../src/colloid.cpp ../tactile/src/tiling.cpp
	$(CPP) $(CPPFLAGS) ../src/colloid.cpp -c -o $@

main.o: main.cpp 
	$(CPP) $(CPPFLAGS) $^ -c -o $@

clean:
	rm -f *.o main
