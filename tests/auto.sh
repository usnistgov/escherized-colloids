#!/bin/bash

# 1. Compile tests.
make clean; make;

# 2. Run tests.
./run_tests --gtest_shuffle --gtest_random_seed=42;
