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
$ cd ../src/
$ make
~~~

Dependencies
------------
This code relies on Prof. Craig Kaplan's [tactile](https://github.com/isohedral/tactile) library, which is cloned and shipped with this code.
You can view an interactive version of this library [here](https://isohedral.ca/software/tactile/).

JSON support for c++ is provided by [here](https://github.com/nlohmann/json). The lone necessary json.hpp file is shipped with this code in
src/, but can be updated from this repo. For example

~~~code
$ wget -O src/json.hpp https://github.com/nlohmann/json/releases/download/v3.10.4/json.hpp
~~~

Example
=======
