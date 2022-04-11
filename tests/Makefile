INCLUDES = -I../tactile/deps/glm \
	-I../tactile/deps/gl3w/build/include \
	-I../tactile/deps/glfw-3.2.1/include \
	-I../src

CC = gcc 
CPP = g++
CFLAGS = $(INCLUDES) -I../tactile/src -g -Wall 
CPPFLAGS = $(INCLUDES) -I../tactile/src -I../src/ -g -Wall -std=c++11 
TACTILE_CPPFLAGS = -I../tactile/deps/glm -I../tactile/src -g -Wall -std=c++11

LIBS = -L../tactile/deps/glfw-3.2.1/build/src \
	-framework OpenGL -framework IOKit -framework CoreVideo

run_tests: run_tests.o colloid.o tiling.o motif.o utils.o
	$(CPP) -o $@ $^ -lgtest -lpthread 

gl3w.o: ../tactile/deps/gl3w/build/src/gl3w.c
	$(CC) $(CFLAGS) $^ -c -o $@

tiling.o: ../tactile/src/tiling.cpp ../tactile/src/tiling_arraydecl.inc
	$(CPP) $(TACTILE_CPPFLAGS) ../tactile/src/tiling.cpp -c -o $@

motif.o: ../src/motif.cpp
	$(CPP) $(CPPFLAGS) ../src/motif.cpp -c -o $@

colloid.o: ../src/colloid.cpp ../tactile/src/tiling.cpp
	$(CPP) $(CPPFLAGS) ../src/colloid.cpp -c -o $@

utils.o: ../src/utils.cpp
	$(CPP) $(CPPFLAGS) ../src/utils.cpp -c -o $@

run_tests.o: run_tests.cpp 
	$(CPP) $(CPPFLAGS) $^ -c -o $@

clean:
	rm -f *.o run_tests
