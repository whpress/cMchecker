struct Person {
	// a person in (or mated to) the genealogy tree
	int idnum, idfa, idmo;
	char idstr[64];
	bool isTree, isCperson, isMperson; // C = candidate node, M = GEDmatched node
	Genome genome;
	Person() : idnum(-1), idfa(-1), idmo(-1), isTree(false), // nchildren(0),
		isCperson(false), isMperson(false) {}
};