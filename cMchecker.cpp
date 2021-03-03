#ifdef _WIN32
#define _CRT_SECURE_NO_WARNINGS
#include "direct.h"
#else
#include <unistd.h>
#define _getcwd getcwd
#endif
#include "nr3a.h"
#include "ran.h"
#include "ludcmp.h"
#include "gamma.h"
#include "incgammabeta.h"
#include "hash.h"
#include "sort.h"

void mythrow(const char *message) {
	printf("%s\n", message);
	exit(1);
}

#include "cMchecker_body.h"
