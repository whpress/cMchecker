#ifndef HAPLOID_H 
#define HAPLOID_H

#include "includeLibs.hpp"
#include "globals.hpp"

struct Haploid {
	// random number generator
	Ran ran;
	const char idtocode[64] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ?";
	// each chromosome is broken into bins
	NRvector<VecInt> chromos; // matrix of [chromosome][bin]
	// used in Genealogy Genome Haploid
	Int &totbins;
	// used in Haploid
	VecInt &binsperchr;
	// used in Genome Haploid
	VecDoub &cmperbin;

	Haploid(Int &totbins, VecInt &binsperchr, VecDoub &cmperbin);
	Haploid(const Haploid& old);
	Haploid &operator=(const Haploid& other);
	void printme();
	void init(Int uniq); // uniq is normally either 2*id (fa) or 2*id+1 (mo)
	Haploid recombine(Haploid &other);
	Int stats3(Int uniq);  // number of runs to uniq
	VecInt stats2(Int uniq);  // run statistics for ancestor uniq (like below, but run to next start)
	VecInt stats(Int uniq);  // run statistics for ancestor uniq
	VecInt stats(Haploid &other);  // run statistics between two haploids
};

// const char Haploid::idtocode[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ?"; // ? is [62]

Haploid::Haploid(Int &totbins, VecInt &binsperchr, VecDoub &cmperbin):
	chromos(NCHR),
	totbins(totbins),
	binsperchr(binsperchr),
	cmperbin(cmperbin)
{
	for (int i = 0; i < NCHR; i++) chromos[i].resize(binsperchr[i]);
}

Haploid::Haploid(const Haploid& old):
	totbins(old.totbins),
	binsperchr(old.binsperchr),
	cmperbin(old.cmperbin)
{ // copy constructor
	for (int i = 0; i < NCHR; i++) {
		chromos[i].resize(binsperchr[i]);
		for (int j = 0; j < chromos[i].size(); j++) chromos[i][j] = old.chromos[i][j];
	}
}

Haploid & Haploid::operator=(const Haploid& other) {
	chromos = other.chromos;
	totbins = other.totbins;
	binsperchr = other.binsperchr;
	cmperbin = other.cmperbin;
	return *this;
}

void Haploid::printme() {
	if (totbins <= 500) {
		for (int i = 0; i < NCHR; i++) {
			if (i>0) printf("|");
			for (int j = 0; j < chromos[i].size(); j++) {
				char ch = idtocode[MIN(62, chromos[i][j]/2)]; // div by 2 to get id
				printf("%c", ch);
			}
		}
		printf("\n");
	}
	else printf("too large to print\n");
}

void Haploid::init(Int uniq) { // uniq is normally either 2*id (fa) or 2*id+1 (mo)
	for (int i = 0; i < NCHR; i++)
		for (int j = 0; j < chromos[i].size(); j++) chromos[i][j] = uniq;
}

Haploid Haploid::recombine(Haploid &other) {
	// does a recombination of two Haploids by the cM probabilities
	Haploid ans = Haploid(totbins, binsperchr, cmperbin);
	Doub prob;
	VecInt *aa, *bb;
	for (int i = 0; i < NCHR; i++) {
		prob = 1. - exp(-cmperbin[i] / 100.); // centiMorgans to prob
		aa = &(chromos[i]);
		bb = &(other.chromos[i]);
		if (ran.doub() < 0.5) SWAP(aa, bb);
		for (int j = 0; j < chromos[i].size(); j++) {
			ans.chromos[i][j] = aa->operator[](j);
			if (ran.doub() < prob) SWAP(aa, bb);
		}
	}
	return ans;
}

Int Haploid::stats3(Int uniq) { // number of runs to uniq
	Int count, sz;
	count = 0;
	for (int i = 0; i < NCHR; i++) {
		sz = chromos[i].size();
		for (int j = 1; j < sz; j++) {
			if (chromos[i][j] == uniq && chromos[i][j - 1] != uniq) { // start of a run
				++count;
			}
		}
	}
	return count;
}
VecInt Haploid::stats2(Int uniq) { // run statistics for ancestor uniq (like below, but run to next start)
	VecInt ans(totbins, 0);
	bool inrun;
	Int count, sz;
	for (int i = 0; i < NCHR; i++) {
		sz = chromos[i].size();
		inrun = false; // no run from start of chromo
		count = 0;
		for (int j = 1; j < sz; j++) {
			if (chromos[i][j] == uniq && chromos[i][j - 1] != uniq) { // start of a run
				if (inrun) {
					++ans[count];
				}
				inrun = true;
				count = 0;
			}
			else ++count;
		}
	}
	return ans;
}
VecInt Haploid::stats(Int uniq) { // run statistics for ancestor uniq
	VecInt ans(totbins, 0);
	bool inrun;
	Int count, sz;
	for (int i = 0; i < NCHR; i++) {
		sz = chromos[i].size();
		inrun = false;
		for (int j = 0; j < sz; j++) {
			if (chromos[i][j] != uniq || j+1 == sz) {
				if (inrun) { // run ended, record it
					++ans[count];
					inrun = false;
				}
			}
			else { // chromos[i][j] == uniq
				if (inrun) ++count;
				else {
					count = 1;
					inrun = true;
				}
			}
		}
	}
	return ans;
}
VecInt Haploid::stats(Haploid &other) { // run statistics between two haploids
	VecInt ans(totbins, 0);
	bool inrun;
	Int count, sz;
	for (int i = 0; i < NCHR; i++) {
		sz = chromos[i].size();
		inrun = false;
		for (int j = 0; j < sz; j++) {
			if (chromos[i][j] != other.chromos[i][j] || j+1 == sz) {
				if (inrun) { // run ended, record it
					++ans[count];
					inrun = false;
				}
			}
			else { // chromos[i][j] == uniq
				if (inrun) ++count;
				else {
					count = 1;
					inrun = true;
				}
			}
		}
	}
	return ans;
}

#endif