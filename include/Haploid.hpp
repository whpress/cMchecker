#ifndef HAPLOID_H 
#define HAPLOID_H

#include "globals.hpp"

struct Haploid {
	NRvector<VecInt> chromos;
	Haploid() : chromos(NCHR) {
		for (int i = 0; i < NCHR; i++) chromos[i].resize(binsperchr[i]);
	}
	Haploid(const Haploid& old) { // copy constructor
		for (int i = 0; i < NCHR; i++) {
			chromos[i].resize(binsperchr[i]);
			for (int j = 0; j < chromos[i].size(); j++) chromos[i][j] = old.chromos[i][j];
		}
	}
	void printme() {
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
	void init(Int uniq) { // uniq is normally either 2*id (fa) or 2*id+1 (mo)
		for (int i = 0; i < NCHR; i++)
			for (int j = 0; j < chromos[i].size(); j++) chromos[i][j] = uniq;
	}
	Haploid recombine(Haploid &other) {
		// does a recombination of two Haploids by the cM probabilities
		Haploid ans;
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
	Int stats3(Int uniq) { // number of runs to uniq
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
	VecInt stats2(Int uniq) { // run statistics for ancestor uniq (like below, but run to next start)
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
	VecInt stats(Int uniq) { // run statistics for ancestor uniq
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
	VecInt stats(Haploid &other) { // run statistics between two haploids
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
};

#endif