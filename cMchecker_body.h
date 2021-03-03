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

Int lenfixedchararray(char **arr) {
	Int maxlen = 100; // stop runaways
	Int len = 0;
	while (arr[len] != NULL && len < maxlen) ++len;
	return len;
}

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

struct Person {
	// a person in (or mated to) the genealogy tree
	int idnum, idfa, idmo;
	char idstr[64];
	bool isTree, isCperson, isMperson; // C = candidate node, M = GEDmatched node
	Genome genome;
	Person() : idnum(-1), idfa(-1), idmo(-1), isTree(false), // nchildren(0),
		isCperson(false), isMperson(false) {}
};

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

struct Genealogy {
	// just a list of all persons defined, for creating, cleanup, etc.
	const static int maxpersons = 1000;
	Int npersons;
	NRvector<Person*> persons;
	Genealogy() : npersons(0), persons(maxpersons) {
		if (totbins == 0) setresolution(3000);
	}
	~Genealogy() {
		for (int i = 0; i < npersons; i++) delete persons[i];
	}
	Person* newPerson(const char *id) {
		// creates and stores a new person with new genome
		if (idfromstr(id) >= 0 && strcmp(id, "anonXS") != 0) { // anonXS is (by obscurity) a reserved name
			printf("attempting to define %s, but that person already defined\n",id);
			mythrow("exiting...");
		}
		Person* ans = new Person;
		ans->idnum = npersons;
		strcpy(ans->idstr, id);
		ans->genome.initgenome(npersons);
		persons[npersons++] = ans;
		return ans;
	}
	void newChild(const char *id, const char *idfa, const char *idmo) {
		// creates and stores new child of fa and mo
		Person* ans = newPerson(id); // , true);
		Int ifa = idfromstr(idfa), imo = idfromstr(idmo);
		if (ifa < 0 || imo < 0) {
			printf("father or mother (%s, %s) undefined for newChild %s\n",idfa,idmo,id);
			mythrow("exiting...");
		}
		ans->idfa = ifa;
		ans->idmo = imo;
		return;
	}
	void newChild(const char *id, const char *idparent, Int ngen=1) {
		// creates and stores new ngenth generation child of parent and a new anons
		Int igen, imo, ifa = idfromstr(idparent);
		if (ifa < 0) {
			printf("parent %s undefined in anon newChild %s\n",idparent,id);
			mythrow("exiting...");
		}
		for (igen = 0; igen < ngen-1; igen++) {
			Person* par = newPerson("anonXS"); // intermediate mother
			imo = par->idnum;
			Person* chld = newPerson("anonXS"); // intermediate child
			chld->idfa = ifa;
			chld->idmo = imo;
			ifa = chld->idnum; // child becomes the next father
		}
		Person* ans = newPerson(id);
		Person* par = newPerson("anonXS"); // final mother
		imo = par->idnum;
		if (ifa < 0 || imo < 0) mythrow("in newChild with ngen, this shouldn\'t happen");
		ans->idfa = ifa;
		ans->idmo = imo;
		return;
	}
	Int idfromstr(const char *str) { // return idnum from idstr
		for (int i = 0; i < npersons; i++) {
			if (strcmp(str, persons[i]->idstr) == 0) return persons[i]->idnum;
		}
		return -1;
	}
	void makeGED(char *filename) {
		Int nis = 0, nfam = 0, idmo, idfa, id;
		VecInt famc(1000), fams1(1000), fams2(1000), famid(1000), famdone(1000, 0);
		Person *per;
		FILE *OUTP = fopen(filename, "wb");
		if (OUTP == NULL) mythrow("error: can\'t open file in makeGED");
		fprintf(OUTP, "0 HEAD\r\n1 SOUR CMC\r\n2 NAME cMchecker\r\n2 VERS 0.0\r\n");
		fprintf(OUTP, "1 SUBM @SUBM@\r\n1 GEDC\r\n2 VERS 5.5.1\r\n");
		fprintf(OUTP, "2 FORM LINEAGE-LINKED\r\n0 @SUBM@ SUBM\r\n1 NAME cMchecker\r\n");
		// this pass populates the families tables
		for (int i = 0; i < npersons; i++) {
			per = persons[i];
			idmo = per->idmo;
			idfa = per->idfa;
			//printf("i=%d id=%d mo=%d fa=%d\r\n",i,per->idnum,idmo,idfa);
			if (idmo < 0 || idfa < 0) continue;
			famid[nfam] = 1000 * (MIN(idmo, idfa)) + MAX(idmo, idfa) + 1;
			famc[nfam] = per->idnum;
			fams1[nfam] = MIN(idmo, idfa);
			fams2[nfam] = MAX(idmo, idfa);
			//printf("%6d: famid=%6d   ch=%6d mo=%6d fa=%6d\r\n", nfam, famid[nfam], famc[nfam], fams1[nfam], fams2[nfam]);
			++nfam;
		}
		// this pass outputs individuals
		for (int i = 0; i < npersons; i++) {
			per = persons[i];
			id = per->idnum;
			fprintf(OUTP, "0 @I%d@ INDI\r\n", id + 1);
			fprintf(OUTP, "1 NAME %s\r\n", per->idstr);
			{ // scope
				Hash<Int, Int, Hashfn2> checkoff(1000, 1000);
				for (int j = 0; j < nfam; j++) {
					if (id == famc[j]) {
						if (checkoff.count(famid[j])) continue; // already done
						fprintf(OUTP, "1 FAMC @F%d@\r\n", famid[j]);
					}
				}
			} // end scope
			{ // scope
				Hash<Int, Int, Hashfn2> checkoff(1000, 1000);
				for (int j = 0; j < nfam; j++) {
					if (id == fams1[j] || id == fams2[j]) {
						if (checkoff.count(famid[j])) continue; // already done
						fprintf(OUTP, "1 FAMS @F%d@\r\n", famid[j]);
						checkoff.set(famid[j], 1);
					}
				}
			} // endscope
		}
		Hash<Int, Int, Hashfn2> checkoff(1000, 1000);
		bool more = true;
		while (more) {
			more = false;
			for (int j = 0; j < nfam; j++) {
				if (checkoff.count(famid[j])) continue;
				//printf("doing %d: famid=%d\n", j, famid[j]);
				// not previously processed :
				more = true; 
				checkoff.set(famid[j], 1);
				fprintf(OUTP, "0 @F%d@ FAM\r\n1 HUSB @I%d@\r\n1 WIFE @I%d@\r\n1 MARR\r\n",
					famid[j], fams1[j] + 1, fams2[j] + 1);
				fprintf(OUTP, "1 CHIL @I%d@\r\n", famc[j] + 1);
				for (int k = j + 1; k < nfam; k++) { // find other children same family
					//printf("trying %d: famid=%d\n",k,famid[k]);
					if (famid[k] == famid[j]) {
						//printf("   and doing it.\n");
						fprintf(OUTP, "1 CHIL @I%d@\r\n", famc[k] + 1);
					}
				}
			}
		}
		fprintf(OUTP, "0 @U1@ SUBM\r\n1 NAME Submitter\r\n0 TRLR\r\n");
		fclose(OUTP);
	}
	void instantiate() { // computes an instance of all the genomes
		Person *per, *permo, *perfa;
		for (int i = 0; i < npersons; i++) {
			per = persons[i];
			if (per->idfa < 0 || per->idmo < 0) continue; // only children get new genomes
			permo = persons[per->idmo];
			perfa = persons[per->idfa];
			per->genome = (perfa->genome).mate(permo->genome); // order is fa.mate(mo)
		}
	}
	void showallgenomes() {
		for (int i = 0; i < npersons; i++) {
			printf("%3d:%s:",i,persons[i]->idstr);
			(persons[i]->genome).showgenome();
		}
	}
	TrialsOutput runtrials(Int ntry, VecInt &cperid, VecInt &mperid) {
		// primary interface
		Int itry, icp, imp, ncp = cperid.size(), nmp = mperid.size();
		TrialsOutput ans(ncp, nmp, ntry);
		ans.cperid = cperid;
		ans.mperid = mperid;
		for (itry = 0; itry < ntry; itry++) {
			instantiate();
			for (icp = 0; icp < ncp; icp++) for (imp = 0; imp < nmp; imp++) {
				ans[icp][imp][itry] = (persons[cperid[icp]]->genome).match(persons[mperid[imp]]->genome);
			}
		}
		return ans;
	}
	TrialsOutput runtrials(Int ntry, char **cpersons, char **mpersons) {
		// interface with names
		Int i, ncp = lenfixedchararray(cpersons), nmp = lenfixedchararray(mpersons);
		VecInt cperid(ncp), mperid(nmp);
		for (i = 0; i < ncp; i++) {
			cperid[i] = idfromstr(cpersons[i]);
			if (cperid[i] < 0) {
				printf("bad Cperson %s in runtrials\n", cpersons[i]);
				mythrow("exiting");
			}
		}
		for (i = 0; i < nmp; i++) {
			mperid[i] = idfromstr(mpersons[i]);
			if (mperid[i] < 0) {
				printf("bad Mperson %s in runtrials\n", mpersons[i]);
				mythrow("exiting...");
			}
		}
		return runtrials(ntry, cperid, mperid);
	}
	TrialsOutput runtrials(Int ntry, char **mpersons) {
		// interface with identical lists (not needed?)
		TrialsOutput ans = runtrials(ntry, mpersons, mpersons);
		ans.iamsymmetric = true; // I don't think this is used
		return ans;
	}
	TrialsOutput runtrials(Int ntry) {
		//run trials using global mpersons_g
		mpersons_g[npersons_g] = NULL;
		return runtrials(ntry, mpersons_g);
	}

	void showtrials(TrialsOutput &trials, char **cpersons,
		char **mpersons, MatDoub *measurements = NULL) {
		Doub measerr = 40.; // out of totcm
		Int icp, imp, itry, ncp = trials.size(), nmp = trials[0].nrows(), ntry = trials[0].ncols();
		Doub sum, sumsq, xum, xumsq, val, tval, errtot, vval; 
		char filout[] = "D:\\Dropbox\\MyDocuments\\Projects_Current\\DNAdoe\\temp.txt";
		for (icp = 0; icp < ncp; icp++) {
			printf("Comparing %16s to:\n", cpersons[icp]);
			for (imp = 0; imp < nmp; imp++) { 
				printf("     (%3d) %16s: ", icp*nmp+imp, mpersons[imp]);
				sum = sumsq = 0.;
				xum = xumsq = 0.;
				for (itry = 0; itry < ntry; itry++) {
					val = trials[icp][imp][itry];
					sum += val;
					sumsq += SQR(val);
					xum += sqrt(val / totcm);
					xumsq += (val / totcm);
				}
				sum /= ntry;
				xum /= ntry;
				sumsq = sqrt(MAX(0.,(sumsq / ntry) - SQR(sum)));
				xumsq = sqrt(MAX(0.,(xumsq / ntry) - SQR(xum)));
				printf("%6.0f +- %4.0f    %6.3f +- %6.3f",sum,sumsq,xum,xumsq);
				if (measurements) {
					val = (*measurements)[icp][imp];
					printf("  v=%8.1f ", val);
					if (val > 0. && xumsq > 0.) {
						vval = sqrt(val / totcm);
						// add errors in quadrature: bio + meas
						errtot = sqrt(SQR(xumsq)+SQR(0.5*(measerr/totcm)/xum));
						tval = (vval - xum) / errtot;
						printf("t=%5.1f ", tval);
						Int nstars = MAX(0, int(abs(tval)) - 1);
						for (int i = 0; i < nstars; i++) printf("*");
						//fprintf(OUTP, "%.1f %.1f\n", sum, val);
					} else printf("t= infinity *****");
				}
				printf("\n");
			}
		}
	}
	void dumptrials(FILE *OUTP, TrialsOutput &trials) {
		Int icp, imp, itry, ncp = trials.size(), nmp = trials[0].nrows(), ntry = trials[0].ncols();
		//printf("DEBUG ncp=%d nmp=%d\n", ncp, nmp);
		for (itry = 0; itry < ntry; itry++) {
			for (icp = 0; icp < ncp; icp++) {
				for (imp = 0; imp < nmp; imp++) {
						fprintf(OUTP,"%.3f ",trials[icp][imp][itry]);
				}
			}
			fprintf(OUTP, "\n");
		}
	}

};

struct DataSet {
	// allows for entry and manipulation of autosomal shared data
	Int ndat;
	VecInt pri, prj;
	VecDoub valij;
	Genealogy &gg;
	DataSet(Genealogy &ggg) : gg(ggg), ndat(0), pri(1000), prj(1000), valij(1000) {}
	void addpair(const char* id1, const char* id2, Doub val) {
		Int ti = gg.idfromstr(id1), tj = gg.idfromstr(id2);
		if (ti < 0 || tj < 0) {
			printf("attempting to add data for %s and %s, but at least one is unknown\n", id1, id2);
			mythrow("exiting...");
		}
		// see if already in data, in which case only change val
		for (int idat = 0; idat < ndat; idat++) {
			if ((pri[idat] == ti && prj[idat] == tj) || (pri[idat] == tj && prj[idat] == ti)) {
				valij[idat] = sqrt(val/totcm);
				return;
			}
		}
		// make new entry
		pri[ndat] = ti;
		prj[ndat] = tj;
		valij[ndat++] = sqrt(val/totcm);
	}
	void moveto(const char *id1, const char *id2) {
		Int src = gg.idfromstr(id1), dest = gg.idfromstr(id2);
		if (src < 0 || dest < 0) {
			printf("attempting move, but %s or %s is an unknown person\n",id1,id2);
			mythrow("exiting...");
		}
		for (int ii = 0; ii < ndat; ii++) if (pri[ii] == dest || prj[ii] == dest) {
			printf("attempting to move data from %s to %s, but %s already has their own data\n",
				id1,id2,id2);
			mythrow("exiting...");
		}
		for (int ii = 0; ii < ndat; ii++) {
			if (pri[ii] == src) pri[ii] = dest;
			if (prj[ii] == src) prj[ii] = dest;
		}
	}
	void omit(const char *id) {
		// omit person id from dataset
		Int nn=0, pr = gg.idfromstr(id);
		if (pr < 0) {
			printf("attempting to omit an unknown person %s from the dataset\n",id);
			mythrow("exiting...");
		}
		// copy everyone except omitted person
		for (int ii = 0; ii < ndat; ii++) if (pri[ii] != pr && prj[ii] != pr) {
			pri[nn] = pri[ii];
			prj[nn] = prj[ii];
			valij[nn] = valij[ii];
			++nn;
		}
		ndat = nn;
	}
};

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

	MultivariateModel(Genealogy &ggg, TrialsOutput &ttrials, DataSet &ddat, Doub measerr = 0.03)
		: gg(ggg), dat(ddat), trials(ttrials), mpr(ddat.ndat), tri(1000), trj(1000) { // max 1000 pairs
		ncp = trials.block.size();
		nmp = trials.block[0].nrows();
		ntrials = trials.block[0].ncols();
		domodel(measerr);
		printf("New model with %d measurements (%d uninformative) and ntrials=%d data...\n",
			dat.ndat, dofreduce, ntrials);
	}

	void domodel(Doub measerr) { // SQR(measerr) is added to covar diagonals
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
	void showmodel() {
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
};

Int parse(Genealogy &gg, DataSet &dd, FILE* &INP) {
	Int nread, nconv;
	char lin[4096], line[4096], token[1024], token2[1024], token3[1024];
	char *tok, *p;
	if (INP == NULL) mythrow("call to parse with bad or unknown input file");
	nread = fscanf(INP, "%[^\r\n]%*[\r\n]", lin);
	strcpy(line, lin);
	if (nread == 1 && strlen(lin) > 0) {
		//printf("%s\n", lin);
		tok = strtok(lin, " ,\t");
		if (tok[0] == '#') {
			printf("%s\n", line);
			return 1;
		}
		strcpy(token, tok);
		for (p=token; *p; ++p) *p = tolower(*p); // idiom for lowercase
		//printf("command is %s\n", token);
		if (strcmp(token, "newperson") == 0 || strcmp(token, "person") == 0) {
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing name in newPerson command");
			strcpy(token, tok);
			gg.newPerson(token);
			printf("tree: %s is a person\n", token);
		}
		else if (strcmp(token, "newchild") == 0 || strcmp(token, "child") == 0) {
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing child name in newChild command");
			strcpy(token, tok);
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing first parent name in newChild command");
			strcpy(token2, tok);
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) {
				gg.newChild(token, token2);
				printf("tree: %s is a child of parent %s and anonymous\n", token, token2);
				return 1; // done, no third token
			}
			strcpy(token3, tok);
			if (isdigit(token3[0])) {
				Int ngen;
				sscanf(token3, "%d", &ngen);
				if (ngen < 2 || ngen > 20) mythrow("bad number of generations in newChild command");
				gg.newChild(token, token2, ngen);
				printf("tree: %s is a %d generation descendent of %s and anonymous\n", token, ngen, token2);
			}
			else {
				gg.newChild(token, token2, token3);
				printf("tree: %s is a child of parents %s and %s\n", token, token2, token3);
			}
		}
		else if (strcmp(token, "personofinterest") == 0 || strcmp(token, "personsofinterest") == 0) {
			while ((tok = strtok(NULL, " ,\t")) != NULL && strlen(tok) > 0) {
				mpersons_g[npersons_g] = new char[strlen(tok) + 1]; // will ignore memory leak
				strcpy(mpersons_g[npersons_g], tok);
				printf("pre-trials: %s is a person of interest\n", mpersons_g[npersons_g]);
				++npersons_g;
			}
		}
		else if (strcmp(token, "allzeros") == 0) {
			for (int i = 0; i < npersons_g; i++) {
				for (int j = i+1; j < npersons_g; j++) {
					dd.addpair(mpersons_g[i], mpersons_g[j], 0.);
				}
			}
		}
		else if (strcmp(token, "makeged") == 0) {
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing filename in makeGED command");
			printf("making GED file %s\n", tok);
			gg.makeGED(tok);
		}
		else if (strcmp(token, "runtrials") == 0) {
			Int ntrials;
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing number of trials in runTrials command");
			sscanf(tok, "%d", &ntrials);
			if (ntrials < 1 || ntrials > 100000) mythrow("number of trials too big or too small");
			printf("starting trials...");
			trials_g = gg.runtrials(ntrials);
			printf("done.\n\n");
		}
		else if (strcmp(token, "newdata") == 0 || strcmp(token, "data") == 0) {
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing first person in newData command");
			strcpy(token, tok);
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing second person in newData command");
			strcpy(token2, tok);
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing value in newData command");
			Doub val;
			nconv = sscanf(tok, "%lf", &val);
			if (nconv != 1 || val < 0. || val > 100000.) mythrow("bad value in newData command");
			dd.addpair(token, token2, val);
			printf("data: %6.1f cM is common between %s and %s\n", val, token, token2);
		}
		else if (strcmp(token, "seterror") == 0) {
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing value in setError comman");
			strcpy(token, tok);
			sscanf(token, "%lf", &measerr_g);
			printf("\nchanging error parameter to %.3f (usual value is 0.03)\n", measerr_g);
		}
		else if (strcmp(token, "done") == 0 || strcmp(token, "exit") == 0) {
			return -2;
		}
		else if (strcmp(token, "show") == 0 || strcmp(token, "showmodel") == 0) {
			printf("\nBaseline model:\n");
			MultivariateModel mm(gg, trials_g, dd, measerr_g);
			logprob_g = mm.logprob;
			chsqprob_g = mm.chisqprob;
			mm.showmodel();
		}
		else if (strcmp(token, "showmodelomitting") == 0) {
			Int nomit = 0;
			DataSet d1(dd); // new dataset is copy of dd
			printf("\n");
			while ((tok = strtok(NULL, " ,\t")) != NULL) {
				d1.omit(tok);
				++nomit;
				printf("Omitting person %s:\n", tok);
			}
			if (nomit == 0) mythrow("missing person in showModelOmitting command");
			MultivariateModel m1(gg, trials_g, d1, measerr_g);
			m1.showmodel();
		}
		else if (strcmp(token, "permanentlyremove") == 0) {
			Int nomit = 0;
			printf("\n");
			while ((tok = strtok(NULL, " ,\t")) != NULL) {
				dd.omit(tok);
				++nomit;
				printf("Removing person %s.\n", tok);
			}
			if (nomit == 0) mythrow("missing person in permanentlyRemove command");
		}
		else if (strcmp(token, "showmodelmoving") == 0) {
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing first person in showModelMoving command");
			strcpy(token, tok);
			tok = strtok(NULL, " ,\t");
			if (tok == NULL) mythrow("missing second person in showModelMoving command");
			strcpy(token2, tok);
			DataSet d2(dd); // new dataset is copy of dd
			d2.moveto(token, token2);
			printf("\nMoving person %s to %s:\n", token, token2);
			MultivariateModel m2(gg, trials_g, d2, measerr_g);
			m2.showmodel();
		}
		else {
			printf("with input line \"%s\"\n", line);
			mythrow("unrecognized command in input");
		}
		return 1; // is good
	}
	return -1; // is not good (e.g., end of file)
}

char helpstring[] = "\
Usage: cMchecker inputfilename\n\n\
Each line of the input file is a keyword and arguments (separated by space or comma).\n\
Keywords (usually in this order):\n\
newPerson personname\n\
newChild childname parentname1 [parentname2] <-- missing 2nd parent means a unique anon.\n\
newChild [great...]grandchildname ancestorname ngen <-- ngen = 2 for grandchild, etc.\n\
makeGED filename <-- outputs a GEDcom file of the tree as entered\n\
personOfInterest name1 [, name2 ...] <-- will have data or be an hypothesis\n\
runTrials ntrials <-- ntrials = 1000 is typical\n\
allZeros <-- initializes \"newData name1 name2 0.0\" for every pair of persons-of-interest\n\
newData name1 name2 value-in-cM <-- sets measured similarity in cM (Method II re 3400 cM for sibs)\n\
showModel <-- makes baseline model, tests data entered so far\n\
showModelOmitting name <-- model excursion omitting a single person (see if they are a problem)\n\
showModelMoving name1 name2 <-- model excusion moving name1's data to name2\n\
permanentlyRemove name <-- removes name from tree\n\
setError value <-- default is 0.03, but can increase to make model more tolerant (lower Chi-Square) \n\
";

//#include "non-production-mains.h"  // all excursions and off-line uses

int main(int argc, char **argv) { // production of standalone executable
	printf("cMchecker: Multivariate Gaussian model for genealogies, version %s\n", version);
	printf("Copyright (C) 2019 DNA Doe Project, Inc.\n\n");
	if (argc < 2) {
		printf("%s", helpstring);
		exit(1);
	}
	char filin[1024], cwd[2048];
	Int flag;
	Genealogy gg;
	DataSet dd(gg);
	strcpy(filin, argv[1]);
	FILE *INP = fopen(filin, "rb");
	if (!INP) {
		_getcwd(cwd, 1023);
		printf("Can\'t find input file %s in working directory %s\n", filin, cwd);
		mythrow("exiting.");
	}
	for (;;) {
		flag = parse(gg, dd, INP);
		//printf("flag=%d\n", flag);
		if (flag != 1) break;
	}
	return 0;
}


