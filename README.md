# Escherized Colloids

Use isohedral tiles to decorate colloidal components.

Installation
============

Getting Started
---------------

~~~code
$ git clone mahynski/escherized_colloids
$ cd escherized_colloids
<!--$ git clone https://github.com/isohedral/tactile-->
$ # Perform required installations, see below
$ # Perform optional installations, see below
$ cd tactile; git checkout 33d5b03 # This is the one used during development
$ cd ../;
~~~

Dependencies
------------

### Tactile (required, included)
This code relies on Prof. Craig Kaplan's [tactile](https://github.com/isohedral/tactile) library, which is cloned and shipped with this code.
You can view an interactive version of this library [here](https://isohedral.ca/software/tactile/).

### JSON (required, included)
JSON support for c++ is provided by [here](https://github.com/nlohmann/json). The lone necessary json.hpp file is shipped with this code in
src/, but can be updated from this repo. For example

~~~code
$ wget -O src/json.hpp https://github.com/nlohmann/json/releases/download/v3.10.4/json.hpp
~~~

### OptimLib (required)
Optimizations are performed using the [OptimLib](https://optimlib.readthedocs.io/en/latest/) library.

~~~code
$ sudo apt install libarmadillo-dev
$ git clone https://github.com/kthohr/optim ./optim
$ cd optim
$ export ARMA_INCLUDE_PATH=/usr/include/; 
$ ./configure -l arma --header-only-version
~~~

See [here](https://optimlib.readthedocs.io/en/latest/installation.html) for more details on installing it as a header only library.

When compiling your code, add `#define OPTIM_ENABLE_ARMA_WRAPPERS` and `#include "optim.hpp` to your cpp file, and set the include path to the head_only_version directory (e.g.,-I/path/to/optimlib/header_only_version) in your Makefile. See the `examples/` directory for example implementations.

You should also perform tests as described [here](https://optimlib.readthedocs.io/en/latest/examples_and_tests.html).

### GoogleTest (optional)
If you want to run tests to check the code (recommended) you need to install [GoogleTest](https://github.com/google/googletest).  This installation procedure also requires [CMake](http://www.cmake.org/) and superuser privileges.

~~~code
$ git clone https://github.com/google/googletest.git
$ cd googletest
$ mkdir build
$ cd build
$ cmake ..
$ make
$ sudo make install
~~~

Once installed, you can automatically execute all tests as follows.

~~~code
$ cd tests
$ bash auto.sh
~~~

### pymatgen (optional)
Symmetry checks and manipulations are performed using [pymatgen](https://pymatgen.org/) in the `design` directory.  This package must be installed if you intend to use the tools therein.

### pre-commit (optional)
This repo uses pre-commit to manage style, as described [here](https://github.com/bmorcos/pre-commit-hooks-cpp).  You will need to install [cpplint](https://pypi.org/project/cpplint/) and [clang-format](https://clang.llvm.org/docs/ClangFormat.html) if you want to contribute. If you do not already have [pre-commit](https://pre-commit.com/) installed, you will need to do that as well.

~~~code
$ pip install cpplint
$ sudo apt install clang-format
~~~

Examples
========
Some basic examples of using the code to perform some basic manipulations have been provided in the `examples/` directory. For example, initializing a colloid by placing a motif inside an isohedral tile, optimizing the tile to fit the motif (``shrinkwrapping''), and creating a unit cell out of motifs.

~~~code
$ cd examples/initialize_colloid/
$ make
$ bash run.sh # You can change parameters in this file as needed.
~~~

How To
======

In the `design/` directory is a jupyter notebook (`Design.ipynb`) which walks the user through the design process.  Briefly, it consists of:

* Step 1a: Build new motif or decide on one from the `motif_libary/`.  The Jupyter notebook in the `motif_library/` directory illustrates how these are built and saved as JSON files.
* Step 1b: It can be a good idea to check the point symmetry you have assigned to your motif using pymatgen.  This optional step is illustrated in the Design.ipynb notebook.
* Step 2: Given the motif's symmetry, determine which groups (tiles) are safe, dangerous, or forbidden.
* Step 3: Select a tile, create a colloid (e.g., see `examples/initialize_colloid/`) by combining it with the chosen motif.
* Step 4: Create a unit cell (2x2 or greater is usually best) to inspect the design for its symmetry.  If you made a **safe** choice, this should be what you wanted. If you
made a **dangerous** choice, you should double check.

To Do
=====
* Logic and coding for 47 non-FD tiles
* DoF for 2 Bezier CP on edges
