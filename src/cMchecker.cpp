#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif

// #include "../include/includeLibs.hpp"
#include "../include/Parse.hpp"

int main(int argc, char **argv) { // production of standalone executable
	printf("cMchecker: Multivariate Gaussian model for genealogies, version %s\n", version);
	printf("Copyright (C) 2019 DNA Doe Project, Inc.\n\n");
	if (argc < 2) {
		Parse::print_help();
		exit(1);
	}
	char filin[1024], cwd[2048];
	Int flag;
	strcpy(filin, argv[1]);
	Parse parse;
	parse.parse_file(filin);
	return 0;
}
