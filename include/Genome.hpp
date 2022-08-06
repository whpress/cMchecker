struct Genome {
	// hides all features of the genome; every person has one
	char ancestry[1024]; // ancestry as a human-readable string
	Haploid fa, mo;

	Doub match(Genome &other) {
	// returns centiMorgan value of match between this and other
		Doub sum = 0.;
		bool famatch, momatch;
		Int thisfa, thismo, otherfa, othermo;
		for (int i = 0; i < NCHR; i++) {
			for (int j = 0; j < fa.chromos[i].size(); j++) {
				thisfa = fa.chromos[i][j];
				otherfa = other.fa.chromos[i][j];
				thismo = mo.chromos[i][j];
				othermo = other.mo.chromos[i][j];
				famatch = (thisfa == otherfa || thisfa == othermo);
				momatch = (thismo == othermo || thismo == otherfa);
				if (WHICHMETHOD == 1) {
					if (famatch) sum += cmperbin[i];
					if (momatch) sum += cmperbin[i];
				}
				else if (WHICHMETHOD == 2) {
					if (famatch || momatch) sum += cmperbin[i];
				}
				else mythrow("this cannot happen in Genome");
			}
		}
		return sum;
	}
	void initgenome(Int id) {
	// initializes genome as from person id
		sprintf(ancestry, "%d", id);
		fa.init(2 * id);
		mo.init(2 * id + 1);
	}
	Genome mate(Genome &other) { // NB: always call as father.mate(mother)
	// returns a genome resulting from mating this and other
		Genome ans;
		sprintf(ans.ancestry, "(%s+%s)", ancestry, other.ancestry);
		ans.fa = fa.recombine(mo);
		ans.mo = other.fa.recombine(other.mo);
		return ans;
	}
	void showgenome() {
	// displays a genome summary (somehow)
		printf("%s\n", ancestry);
		fa.printme();
		mo.printme();
	}
	VecInt stats(Haploid &other) { // run statistics between a genome and a haploid (e.g. for uncle)
		VecInt ans(totbins, 0);
		bool inrun;
		Int count, sz;
		for (int i = 0; i < NCHR; i++) {
			sz = fa.chromos[i].size();
			inrun = false;
			for (int j = 0; j < sz; j++) {
				// count match to either haploid as continuing a run
				bool matches = (fa.chromos[i][j] == other.chromos[i][j]
					|| mo.chromos[i][j] == other.chromos[i][j]);
				if (! matches || j+1 == sz) {
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