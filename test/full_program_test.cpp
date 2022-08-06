#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif
#include "../nr3a.h"
#include "../ran.h"
#include "../ludcmp.h"
#include "../gamma.h"
#include "../incgammabeta.h"
#include "../hash.h"
#include "../sort.h"

void mythrow(const char *message) {
	printf("%s\n", message);
	exit(1);
}

#include "../cMchecker_body.h"
#include <gtest/gtest.h>

void simulate_body_main(int argc, char **argv, Genealogy &gg, DataSet &dd) {
  char filin[1024], cwd[2048];
	Int flag;
  // Genealogy gg;
	// DataSet dd(gg);
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
}

TEST(FullProgramTest, TestUncleRelationship) {
  
  int input_argc = 2;
  char *input_argv[] = {
    (char*)"",
    (char*)"test/input_files/uncle_newphew_relationship.txt",
    NULL
  };

  Genealogy gg;
  DataSet dd(gg);

  simulate_body_main(input_argc, input_argv, gg, dd);

  MultivariateModel mm(gg, trials_g, dd, measerr_g);
  std::cout << "RESULTING chisqprob: " << mm.chisqprob << endl;

  EXPECT_TRUE((mm.chisqprob > 0.80) && (mm.chisqprob < 0.90));
}