#ifndef MYTHROW_H
#define MYTHROW_H

#include "includeLibs.hpp"

void mythrow(const char *message) {
	printf("%s\n", message);
	exit(1);
}

#endif
