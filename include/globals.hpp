#ifndef GLOBALS_H
#define GLOBALS_H

#define version "0.018 (beta)"

Int NCHR = 22; // autosomal only
#define WHICHMETHOD 1 // method 1 or 2, see https://isogg.org/wiki/Autosomal_DNA_statistics
Doub cmperchr[22] = { 284.,269.,223.,214.,204.,192.,187.,168.,166.,
	181.,158.,175.,126.,119.,141.,134.,128.,117.,108.,108.,62.7,72.7};
	// using 23andMe values from https://isogg.org/wiki/CentiMorgan
Doub twinsvalue = 3400.; // renormalizes sum of cmperchr to this
Doub totcm;
Int totbins = 0;
VecInt binsperchr;
VecDoub cmperbin;
Doub meancmperbin;
Doub logprob_g = 0.;
Doub chsqprob_g = 1.;
Doub measerr_g = 0.03; // default "measurement error"
char idtocode[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ?"; // ? is [62]
Ran ran;
char *mpersons_g[1000]; // global array of mpersons names
Int npersons_g = 0; // global number of mpersons

#endif