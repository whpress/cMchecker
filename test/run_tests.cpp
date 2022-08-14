#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif

#include "../include/includeLibs.hpp"
#include "../include/Parse.hpp"
#include "../include/DataSet.hpp"
#include "../include/Genealogy.hpp"
#include "../include/myThrow.hpp"

#include <gtest/gtest.h>

int main (int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

Parse simulate_body_main(int argc, char **argv) {
  char filin[1024], cwd[2048];
	Int flag;
  // Genealogy gg;
	// DataSet dd(gg);
  if (argc < 2) {
    Parse::print_help();
    exit(1);
  }
	strcpy(filin, argv[1]);
  Parse parse;
  parse.set_allow_output(false);
  parse.parse_file(filin);
  return parse;
	// FILE *INP = fopen(filin, "rb");
	// if (!INP) {
	// 	_getcwd(cwd, 1023);
	// 	printf("Can\'t find input file %s in working directory %s\n", filin, cwd);
	// 	mythrow("exiting.");
	// }
	// for (;;) {
	// 	flag = parse(gg, dd, INP);
	// 	//printf("flag=%d\n", flag);
	// 	if (flag != 1) break;
	// }
}

double obtain_chisquared(std::string input_file) {
  int input_argc = 2;
  char *input_argv[] = {
    (char*)"",
    (char*)input_file.c_str(),
    NULL
  };

  Parse parse = simulate_body_main(input_argc, input_argv);

  MultivariateModel mm = parse.get_multivariate_model();
  return mm.chisqprob;
}

#include "test_ideal_values.cpp"
#include "test_avg_values.cpp"
#include "test_bad_values.cpp"