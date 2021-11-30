/**
 * Copyright 2022 Nathan A. Mahynski
 * @author Nathan A. Mahynski
 * 
 * General utilities and variables with global scope.
 */

#ifndef SRC_UTILS_HPP_
#define SRC_UTILS_HPP_

#include <iostream>
#include <exception>
#include <vector>
#include <string>

using std::string;
using std::exception;
using std::endl;

#define EXIT_FAILURE 1

// Isohedral tile types that are fundamental domains of crystals.
const int FD_TYPES[46] = {1, 2, 3, 4, 5, 6, 7, 21, 22, 23, 24, 25, 27, 28, 30, 31, 33, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 77, 78, 79, 80, 81, 83, 84, 85, 86, 87, 88};

class customException { // Adapted from https://www.oreilly.com/library/view/c-cookbook/0596007612/ch09s02.html
 public:
	explicit customException(const string& msg) : msg_(msg) {}
	~customException() {}

	string getMessage() const { return(msg_); }
 private:
	string msg_;
};

#endif // SRC_UTILS_HPP_
