#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif
#include "../libs/nr3a.h"
#include "../libs/ran.h"
#include "../libs/ludcmp.h"
#include "../libs/gamma.h"
#include "../libs/incgammabeta.h"
#include "../libs/hash.h"
#include "../libs/sort.h"

void mythrow(const char *message) {
	printf("%s\n", message);
	exit(1);
}

#include "../include/globals.hpp"
#include "../include/lenFixedCharArray.hpp"
#include "../include/Haploid.hpp"
#include "../include/Genome.hpp"
#include "../include/Person.hpp"
#include "../include/TrialsOutput.hpp"
#include "../include/setResolution.hpp"
#include "../include/Genealogy.hpp"
#include "../include/DataSet.hpp"
#include "../include/MultivariateModel.hpp"
#include "../include/parse.hpp"

#include "../include/cMchecker_body.hpp"
