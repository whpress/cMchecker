#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif

#include "../include/includeLibs.hpp"

void mythrow(const char *message) {
	printf("%s\n", message);
	exit(1);
}

#include "../include/cMchecker_body.hpp"
