#ifndef MULTIVARIATEMODEL_H
#define MULTIVARIATEMODEL_H

#include "includeLibs.hpp"
#include "DataSet.hpp"
#include "Genealogy.hpp"
#include "TrialsOutput.hpp"
#include "globals.hpp"
#include "myThrow.hpp"

struct MultivariateModel {
	// model brings together a genealogy, a set of trials (theoretical), and a set of data (measured)
	Genealogy &gg;
	TrialsOutput &trials;
	DataSet &dat;
	VecDoub mu, pctile;
	MatDoub sigma, sigmainv;
	Int mpr, ncp, nmp, ntrials, dofreduce; // mpr is number of pairs
	VecInt tri, trj; // persons index into trials
	Doub chisq, chisqprob, tsqsum, logprob;
	// used in MultivariateModel and parse refactored from globals
	Doub logprob_g = 0.;
	// used in MultivariateModel and parse refactored from globals
	Doub chsqprob_g = 1.;

	MultivariateModel(Genealogy &ggg, TrialsOutput &ttrials, DataSet &ddat, Doub measerr);
	void domodel(Doub measerr); // SQR(measerr) is added to covar diagonals
	void showmodel();
};

MultivariateModel::MultivariateModel(Genealogy &ggg, TrialsOutput &ttrials, DataSet &ddat, Doub measerr = 0.03):
	gg(ggg),
	dat(ddat),
	trials(ttrials),
	mpr(ddat.ndat),
	tri(1000),
	trj(1000)
{ // max 1000 pairs
	ncp = trials.block.size();
	nmp = trials.block[0].nrows();
	ntrials = trials.block[0].ncols();
	domodel(measerr);
	printf("New model with %d measurements (%d uninformative) and ntrials=%d data...\n",
		dat.ndat, dofreduce, ntrials);
}

void MultivariateModel::domodel(Doub measerr) { // SQR(measerr) is added to covar diagonals
	Int i, j, k, ii, jj, itry, ti, tj, ui, uj;
	Doub sum, sd, tval, val;
	// convert trials to sqrts (because more Normally distributed)
	TrialsOutput sqtrials(ncp, nmp, ntrials);
	for (i = 0; i < ncp; i++) for (j = 0; j < nmp; j++) for (k = 0; k < ntrials; k++) {
		sqtrials[i][j][k] = sqrt(MAX(0., trials[i][j][k]/totcm));
	}
	mu.resize(mpr);
	pctile.resize(mpr);
	sigma.resize(mpr, mpr);
	sigmainv.resize(mpr, mpr);
	VecDoub valminusmu(mpr);
	for (ii = 0; ii < mpr; ii++) { // where are the mpersons in trials block?
		i = trials.cperid.index(dat.pri[ii]); 
		j = trials.mperid.index(dat.prj[ii]);
		if (i < 0 || j < 0) {
			printf("attempting a model with someone who was not a person of interest when trials run\n");
			mythrow("exiting...");
		}
		tri[ii] = i;
		trj[ii] = j;
	}
	// calculate percentile of each measured value in its trials
	for (ii = 0; ii < mpr; ii++) {
		val = dat.valij[ii];
		ti = tri[ii];
		tj = trj[ii];
		for (i = 0, itry = 0; itry < ntrials; itry++) {
			if (sqtrials[ti][tj][itry] > val) ++i;
		}
		pctile[ii] = Doub(i) / Doub(ntrials);
	}
	// calculate sample means of the trials
	for (ii = 0; ii < mpr; ii++) {
		sum = 0.;
		ti = tri[ii];
		tj = trj[ii];
		for (itry = 0; itry < ntrials; itry++) {
			sum += sqtrials[ti][tj][itry];
		}
		mu[ii] = sum / ntrials;
	}
	// calculate sample covariances of the trials
	for (ii = 0; ii < mpr; ii++) {
		ti = tri[ii];
		tj = trj[ii];
		for (jj = ii; jj < mpr; jj++) { // upper triangle
			sum = 0.;
			ui = tri[jj];
			uj = trj[jj];
			for (itry = 0; itry < ntrials; itry++) {
				sum += (sqtrials[ti][tj][itry]-mu[ii])*(sqtrials[ui][uj][itry]-mu[jj]);
			}
			sigma[ii][jj] = sigma[jj][ii] = sum / ntrials;
		}
		sigma[ii][ii] += SQR(measerr);
	}
	// matrix invert covariances to get the model
	LUdcmp LUsigma(sigma);
	LUsigma.inverse(sigmainv);
	Doub detsigma = LUsigma.det();
	dofreduce = 0;
	tsqsum = 0.;
	for (ii = 0; ii < mpr; ii++) {
		valminusmu[ii] = dat.valij[ii] - mu[ii];
		if (dat.valij[ii] == 0 && mu[ii] == 0) ++dofreduce;
		// sum of squares of tvals is model chi-sq *without* accounting for dependencies
		sd = sqrt(MAX(0., sigma[ii][ii]));
		tval = (dat.valij[ii] - mu[ii]) / sd;
		tsqsum += SQR(tval);
	}
	if (mpr - dofreduce < 1) mythrow("Error: Forgot to enter data?  Or all data of unrelated persons?");
// calculate chi-sq *with* accounting for dependencies
	chisq = dotproduct(valminusmu,matmul(sigmainv, valminusmu));
	chisqprob = MAX(1.e-30,1. - Chisqdist(mpr-dofreduce).cdf(chisq));
	// calculate model logprob
	//logprob = (-0.5*chisq - 0.5*(mpr - dofreduce)*1.837 - 0.5*log(detsigma))/2.3026; // ln(2*pi), ln(10.)
	logprob = (-0.5*chisq + 0.5*(mpr - dofreduce))/2.3026; // see paper note
}
void MultivariateModel::showmodel() {
	Doub sd, tval;
	for (int ii = 0; ii < mpr; ii++) {
		char *istr = gg.persons[dat.pri[ii]]->idstr;
		char *jstr = gg.persons[dat.prj[ii]]->idstr;
		sd = sqrt(MAX(0., sigma[ii][ii]));
		tval = (dat.valij[ii] - mu[ii]) / sd;
		if (dat.valij[ii] != 0. || mu[ii] != 0.) {
			Doub xmeas=SQR(dat.valij[ii])*3400., xlo=SQR(MAX(0.,mu[ii]-sd))*3400.,
				xhi=SQR(mu[ii]+sd)*3400;
			printf("%3d: (%-12s %-12s) model=%.3f+-%.3f (%4.0f-%4.0f cM) data=%.3f (%4.0f cM) pct=%5.3f",
				ii, istr, jstr, mu[ii], sd, xlo, xhi, dat.valij[ii], xmeas, pctile[ii]);
			printf(" tval=%5.1f ", tval);
			Int nstars = MIN(5,MAX(0, int(abs(tval)) - 1));
			for (int i = 0; i < nstars; i++) printf("*");
			printf("\n");
		}
	}
	printf("model ChiSquare = %.2f (%.2f) with DoF=%d ChiSquareProbability = %.4g ",
		chisq, tsqsum, mpr-dofreduce, MAX(1.e-20,chisqprob));
	if (chisqprob > 0.001) printf("\n");
	else if (chisqprob > 1.e-6) printf("**\n");
	else printf("*****\n");
	if (logprob - logprob_g != 0.) {
		printf("Comparison:  Chi-Square tail probability ratio of this model to baseline is %.4g\n\
Comparison: MV-Gaussian Bayes likelihood ratio of this model to baseline is %.4g\n\n",
		chisqprob / (chsqprob_g + 1.e-20), pow(10., logprob - logprob_g));
	}
	//printmat(sigma, "%8.1e ");
}

#endif