#ifndef GLOBALS_H
#define GLOBALS_H

#define version "0.018 (beta)"

#include "includeLibs.hpp"

// used in Genome Haploid and setResolution
const Int NCHR = 22; // autosomal only
#define WHICHMETHOD 1 // method 1 or 2, see https://isogg.org/wiki/Autosomal_DNA_statistics
// used only in setResolution
const Doub twinsvalue = 3400.; // renormalizes sum of cmperchr to this
// used in parse and full_program_test
const Doub measerr_g = 0.03; // default "measurement error"

#endif