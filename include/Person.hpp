#ifndef PERSON_H
#define PERSON_H

#include "includeLibs.hpp"
#include "Genome.hpp"

struct Person {
	// a person in (or mated to) the genealogy tree
	int idnum, idfa, idmo;
	char idstr[64];
	bool isTree, isCperson, isMperson; // C = candidate node, M = GEDmatched node
	Genome genome;
	Person();
};

Person::Person():
	idnum(-1),
	idfa(-1),
	idmo(-1),
	isTree(false),
	isCperson(false),
	isMperson(false)
	{}

#endif