#ifndef LENFIXEDCHARARRAY_H
#define LENFIXEDCHARARRAY_H

#include "includeLibs.hpp"

Int lenfixedchararray(char **arr) {
	Int maxlen = 100; // stop runaways
	Int len = 0;
	while (arr[len] != NULL && len < maxlen) ++len;
	return len;
}

#endif