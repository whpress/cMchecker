#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif

#include "../include/includeLibs.hpp"
#include "../include/parse.hpp"
#include "../include/DataSet.hpp"
#include "../include/Genealogy.hpp"
#include "../include/myThrow.hpp"

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

double obtain_chisquared(std::string input_file) {
  int input_argc = 2;
  char *input_argv[] = {
    (char*)"",
    (char*)input_file.c_str(),
    NULL
  };

  Genealogy gg;
  DataSet dd(gg);

  simulate_body_main(input_argc, input_argv, gg, dd);

  MultivariateModel mm(gg, trials_g, dd, measerr_g);
  return mm.chisqprob;
}

TEST(FullProgramTest, TestFatherSonRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/father_son_relationship.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared == 1);
}

TEST(FullProgramTest, TestUncleNephewRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/uncle_nephew_relationship.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared > 0.95);
}

TEST(FullProgramTest, TestSiblingRelationship) {
  double returned_chisquared = obtain_chisquared("test/input_files/sibling_relationship.txt");
  std::cout << "Resulting chi^2 probability: " << returned_chisquared << endl;
  EXPECT_TRUE(returned_chisquared > 0.95);
}
