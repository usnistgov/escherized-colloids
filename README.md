# Escherized Colloids

Use isohedral tiles to produce patchy colloids which self-assemble into two-dimensional films with a desired symmetry group and void fraction.

Installation
============

Getting Started
---------------

~~~code
$ git clone https://github.com/mahynski/escherized-colloids.git
$ cd escherized_colloids
$ # 1. Perform required installations, see below
$ # 2. Perform optional installations, see below
~~~

Dependencies
------------

### Tactile (required, included)
This code relies on Prof. Craig Kaplan's [tactile](https://github.com/isohedral/tactile) library.
You can view an interactive version of this library [here](https://isohedral.ca/software/tactile/).
~~~code
$ git clone https://github.com/isohedral/tactile ./tactile
$ cd tactile; git checkout 33d5b03 # This is the commit used during development
$ cd ../;
~~~

### JSON (required, included)
JSON support for c++ is provided by [here](https://github.com/nlohmann/json). The lone necessary json.hpp file is shipped with this code in
src/, but can be updated from this repo. For example:

~~~code
$ wget -O src/json.hpp https://github.com/nlohmann/json/releases/download/v3.10.4/json.hpp
~~~

### OptimLib (required)
Optimizations are performed using the [OptimLib](https://optimlib.readthedocs.io/en/latest/) library.  On Ubuntu you can install the requirements as follows:

~~~code
$ sudo apt install libopenblas-dev
$ sudo apt install libomp-dev
$ sudo apt install libarmadillo-dev
$ git clone https://github.com/kthohr/optim ./optim
$ cd optim
$ git submodule update --init
$ export ARMA_INCLUDE_PATH=/usr/include/;
$ ./configure ---header-only-version
~~~

See [here](https://optimlib.readthedocs.io/en/latest/installation.html) for more details on installing it as a header only library.

When compiling your code, add `#define OPTIM_ENABLE_ARMA_WRAPPERS` and `#include "optim.hpp` to your cpp file, and set the include path to the head_only_version directory (e.g.,-I/path/to/optimlib/header_only_version) in your Makefile. See the `examples/` directory for example implementations.

You should also perform tests as described [here](https://optimlib.readthedocs.io/en/latest/examples_and_tests.html).

### GoogleTest (optional)
If you want to run tests to check the code (recommended) you need to install [GoogleTest](https://github.com/google/googletest).  This installation procedure also requires [CMake](http://www.cmake.org/) and superuser privileges. Note: this project currently is set up for c++11, however, googletest >= v1.14 has moved to c++14 by default, so use the last c++11 release (v1.12.1) for backward compatibility.

~~~code
$ git clone https://github.com/google/googletest.git --branch release-1.12.1 
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
This repo uses pre-commit to manage style, as described [here](https://github.com/bmorcos/pre-commit-hooks-cpp).  You will need to install [cpplint](https://pypi.org/project/cpplint/) if you want to contribute. If you do not already have [pre-commit](https://pre-commit.com/) installed, you will need to do that as well.

~~~code
$ pip install cpplint
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
A [conda](https://anaconda.org/) environment yaml file is included to reproduce the development environment.  You can install this as follows:

```code
$ conda env create -f conda-env.yml
$ conda activate escherc
$ python -m ipykernel install --user --name=escherc
```

The design procedure is as follows:

* Step 1: Build a new motif or decide on one from the `motif_libary/`.  The Jupyter notebook in the `motif_library/` directory (Create_Motifs.ipynb) illustrates how these are built and saved as JSON files.
* Step 2: Use the Design.ipynb notebook in `design/` to determine which groups (tiles) are safe, dangerous, or forbidden and decide what you want to create.  It is a good idea to check the point symmetry you have assigned to your motif using pymatgen.  This optional step is illustrated in the Design.ipynb notebook.
* Step 3: Select a tile, then create a colloid by placing the motif inside of it (e.g., see `examples/initialize_colloid/`); you can do this in an "optimal" way by optimizing the fit (e.g., see `examples/pso_auto/`).
* Final check: Create a unit cell (2x2 or greater is usually best) to inspect the design for its symmetry.  If you made a **safe** choice, this should be what you wanted. If you made a **dangerous** choice, you should double check (e.g., `examples/unit_cell/`).
* Step 4: See `simulations/` directory for tools that use [LAMMPS](https://www.lammps.org/) to simulate the self-assembly of these colloidal particles.  In particular, the Preparing.ipynb notebook illustrates how to build the simulation scripts (automatically) and also analyze the results.

Citation
========
Refer to the CITATION.cff file for information regarding the citation of this repository.  

Contact information
===================
* Nathan A. Mahynski, Material Measurement Laboratory, Chemical Sciences Division, Chemical Informatics Group. Contact: nathan.mahynski@nist.gov

License Information
===================
* See LICENSE.md for more information.
* Any mention of commercial products is for information only; it does not imply recommendation or endorsement by [NIST](https://www.nist.gov/).
