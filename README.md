# escherized_colloids

Use isohedral tiles to decorate colloidal components.

Installation
============

Getting Started
---------------

~~~code
$ git clone mahynski/escherized_colloids
$ cd escherized_colloids
$ git clone https://github.com/isohedral/tactile
$ cd tactile; git checkout 33d5b03 # This is the one used during development
$ cd ../;
~~~

For optimizations, we must install pagmo - see instructions below for details.

Dependencies
------------

Tactile
=======
This code relies on Prof. Craig Kaplan's [tactile](https://github.com/isohedral/tactile) library, which is cloned and shipped with this code.
You can view an interactive version of this library [here](https://isohedral.ca/software/tactile/).

JSON
====
JSON support for c++ is provided by [here](https://github.com/nlohmann/json). The lone necessary json.hpp file is shipped with this code in
src/, but can be updated from this repo. For example

~~~code
$ wget -O src/json.hpp https://github.com/nlohmann/json/releases/download/v3.10.4/json.hpp
~~~

OptimLib
========
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

Pre-commit
==========
This repo uses pre-commit to manage style, as described [here](https://github.com/bmorcos/pre-commit-hooks-cpp).  You will need to install [cpplint](https://pypi.org/project/cpplint/) and [clang-format](https://clang.llvm.org/docs/ClangFormat.html) if you want to contribute.
If you do not already have [pre-commit](https://pre-commit.com/) installed, you will need to do that as well.

~~~code
$ pip install cpplint
$ sudo apt install clang-format
~~~

Example
=======
~~~code
$ cd examples/initialize_colloid/
$ make
$ bash run.sh # You can change parameters in this file as needed.
~~~

To Do
=====
1. Add code standardization
2. Add unittests
3. Logic and coding for 47 non-FD tiles
