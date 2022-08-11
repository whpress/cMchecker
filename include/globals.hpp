#ifndef GLOBALS_H
#define GLOBALS_H

#define version "0.018 (beta)"

#include "includeLibs.hpp"

// used in Genome Haploid and setResolution
const Int NCHR = 22; // autosomal only
#define WHICHMETHOD 1 // method 1 or 2, see https://isogg.org/wiki/Autosomal_DNA_statistics
// used only in setResolution
Doub cmperchr[22] = { 284.,269.,223.,214.,204.,192.,187.,168.,166.,
	181.,158.,175.,126.,119.,141.,134.,128.,117.,108.,108.,62.7,72.7};
	// using 23andMe values from https://isogg.org/wiki/CentiMorgan
// used only in setResolution
const Doub twinsvalue = 3400.; // renormalizes sum of cmperchr to this
// used in DataSet Genealogy MultivariateModel and setResolution
Doub totcm;
// used in Genealogy Genome Haploid and setResolution
Int totbins = 0;
// used in Haploid and setResolution
VecInt binsperchr;
// used in Genome Haploid and setResolution
VecDoub cmperbin;
// used only in setResolution
Doub meancmperbin;
// used in MultivariateModel and parse
Doub logprob_g = 0.;
// used in MultivariateModel and parse
Doub chsqprob_g = 1.;
// used in parse and full_program_test
Doub measerr_g = 0.03; // default "measurement error"
// used in Genealogy and parse
char *mpersons_g[1000]; // global array of mpersons names
// used in Genealogy and parse
Int npersons_g = 0; // global number of mpersons

#endif