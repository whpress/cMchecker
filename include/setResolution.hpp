#ifndef SETRESOLUTION_H
#define SETRESOLUTION_H

#include "includeLibs.hpp"
#include "globals.hpp"

void setresolution(Int totb) {
	// given total number of bins in genome, get integer bins for each chr
	binsperchr.resize(NCHR);
	cmperbin.resize(NCHR);
	totbins = totb;
	totcm = 0.;
	for (int i = 0; i < NCHR; i++) totcm += cmperchr[i];
	for (int i = 0; i < NCHR; i++) cmperchr[i] *= (twinsvalue / totcm);
	totcm = twinsvalue;
	meancmperbin = totcm / Doub(totbins);
	Doub cmleft = totcm;
	Int binsleft = totbins;
	for (int i = NCHR - 1; i > 0; i--) { // yes, > is correct
		binsperchr[i] = int((cmperchr[i] * binsleft / cmleft) + 0.5); // round to nearest int
		binsleft -= binsperchr[i];
		cmleft -= cmperchr[i];
		cmperbin[i] = cmperchr[i] / binsperchr[i];
	}
	binsperchr[0] = binsleft;
	cmperbin[0] = cmleft / binsleft;
	if (0) for (int i = 0; i < NCHR; i++) printf("%2d: %.3f %4d %.3f %.3f\n",
		i + 1, cmperchr[i], binsperchr[i],cmperbin[i],meancmperbin);
}

#endif