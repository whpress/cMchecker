#ifndef CMCHECKER_H
#define CMCHECKER_H

#include "parse.hpp"
#include "DataSet.hpp"
#include "Genealogy.hpp"

//#include "non-production-mains.h"  // all excursions and off-line uses

int main(int argc, char **argv) { // production of standalone executable
	printf("cMchecker: Multivariate Gaussian model for genealogies, version %s\n", version);
	printf("Copyright (C) 2019 DNA Doe Project, Inc.\n\n");
	if (argc < 2) {
		printf("%s", helpstring);
		exit(1);
	}
	char filin[1024], cwd[2048];
	Int flag;
	Genealogy gg;
	DataSet dd(gg);
	strcpy(filin, argv[1]);
	FILE *INP = fopen(filin, "rb");
	if (!INP) {
		_getcwd(cwd, 1023);
		printf("Can\'t find input file %s in working directory %s\n", filin, cwd);
		mythrow("exiting.");
	}
	for (;;) {
		flag = parse(gg, dd, INP);
		//printf("flag=%d\n", flag);
		if (flag != 1) break;
	}
	return 0;
}

#endif
