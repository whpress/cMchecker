struct TrialsOutput {
	NRvector<MatDoub> block;
	VecInt cperid, mperid;
	bool iamsymmetric;
	TrialsOutput(Int ncpeople, Int nmpeople, Int ntrials) : block(ncpeople), iamsymmetric(false) {
		for (int i = 0; i < ncpeople; i++) block[i].resize(nmpeople, ntrials);
	}
	MatDoub operator [](int i) const {return block[i];}
    MatDoub & operator [](int i) {return block[i];}
	Int size() { return block.size(); }
};
TrialsOutput trials_g(0,0,0); // a global TrialsOutput