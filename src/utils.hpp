#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <exception>

using namespace std;

#define EXIT_FAILURE 1

class customException {
 // Adapted from https://www.oreilly.com/library/view/c-cookbook/0596007612/ch09s02.html
public:
	customException(const string& msg) : msg_(msg) {}
	~customException() {}

	string getMessage() const {return(msg_);}
private:
	string msg_;
};

#endif
