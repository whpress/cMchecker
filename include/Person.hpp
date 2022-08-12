#ifndef PERSON_H
#define PERSON_H

#include "includeLibs.hpp"
#include "Genome.hpp"

struct Person {
	// a person in (or mated to) the genealogy tree
	int idnum, idfa, idmo;
	char idstr[64];
	bool isTree, isCperson, isMperson; // C = candidate node, M = GEDmatched node
	// // used in Genealogy Genome Haploid
	// Int &totbins;
	// // used in Haploid
	// VecInt &binsperchr;
	// // used in Genome Haploid
	// VecDoub &cmperbin;
	Genome genome;

	Person(Int &totbins, VecInt &binsperchr, VecDoub &cmperbin);
};

Person::Person(Int &totbins, VecInt &binsperchr, VecDoub &cmperbin):
	idnum(-1),
	idfa(-1),
	idmo(-1),
	isTree(false),
	isCperson(false),
	isMperson(false),
	// totbins(totbins),
	// binsperchr(binsperchr),
	// cmperbin(cmperbin),
	genome(Genome(totbins, binsperchr, cmperbin))
	{}

#endif