#ifndef GENEALOGY_H 
#define GENEALOGY_H

#include "includeLibs.hpp"
#include "lenFixedCharArray.hpp"
#include "Person.hpp"
#include "TrialsOutput.hpp"
#include "globals.hpp"
#include "myThrow.hpp"

struct Genealogy {
	// just a list of all persons defined, for creating, cleanup, etc.
	const static int maxpersons = 1000;
	Int npersons;
	NRvector<Person*> persons;
	// used in Genealogy and parse refactored from global
	char *mpersons_g[1000]; // global array of mpersons names
	// used in Genealogy and parse refactored from global
	Int npersons_g = 0; // global number of mpersons
	// used in DataSet Genealogy MultivariateModel and setResolution refactored from globals
	Doub totcm;
	// used only in setResolution pointless???
	Doub cmperchr[22] = { 284.,269.,223.,214.,204.,192.,187.,168.,166.,
		181.,158.,175.,126.,119.,141.,134.,128.,117.,108.,108.,62.7,72.7};
		// using 23andMe values from https://isogg.org/wiki/CentiMorgan
	// used in Genealogy Genome Haploid and setResolution
	Int totbins = 0;
	// used in Haploid and setResolution
	VecInt binsperchr;
	// used in Genome Haploid and setResolution
	VecDoub cmperbin;
	// used only in setResolution
	Doub meancmperbin;

	Genealogy();
	~Genealogy();
	Person* newPerson(const char *id);
	void newChild(const char *id, const char *idfa, const char *idmo);
	void newChild(const char *id, const char *idparent, Int ngen=1);
	Int idfromstr(const char *str); // return idnum from idstr
	void makeGED(string filename);
	void instantiate(); // computes an instance of all the genomes
	void showallgenomes();
	TrialsOutput runtrials(Int ntry, VecInt &cperid, VecInt &mperid);
	TrialsOutput runtrials(Int ntry, char **cpersons, char **mpersons);
	TrialsOutput runtrials(Int ntry, char **mpersons);
	TrialsOutput runtrials(Int ntry);
	void showtrials(TrialsOutput &trials, char **cpersons, char **mpersons, MatDoub *measurements = NULL);
	void dumptrials(FILE *OUTP, TrialsOutput &trials);
	void setresolution(Int totb);
};

Genealogy::Genealogy():
	npersons(0),
	persons(maxpersons) 
{
	if (totbins == 0) setresolution(3000);
}
Genealogy::~Genealogy() {
	for (int i = 0; i < npersons; i++) delete persons[i];
}
Person* Genealogy::newPerson(const char *id) {
	// creates and stores a new person with new genome
	if (idfromstr(id) >= 0 && strcmp(id, "anonXS") != 0) { // anonXS is (by obscurity) a reserved name
		printf("attempting to define %s, but that person already defined\n",id);
		mythrow("exiting...");
	}
	Person* ans = new Person(totbins, binsperchr, cmperbin);
	ans->idnum = npersons;
	strcpy(ans->idstr, id);
	ans->genome.initgenome(npersons);
	persons[npersons++] = ans;
	return ans;
}
void Genealogy::newChild(const char *id, const char *idfa, const char *idmo) {
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
void Genealogy::newChild(const char *id, const char *idparent, Int ngen) {
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
Int Genealogy::idfromstr(const char *str) { // return idnum from idstr
	for (int i = 0; i < npersons; i++) {
		if (strcmp(str, persons[i]->idstr) == 0) return persons[i]->idnum;
	}
	return -1;
}
void Genealogy::makeGED(string filename) {
	Int nis = 0, nfam = 0, idmo, idfa, id;
	VecInt famc(1000), fams1(1000), fams2(1000), famid(1000), famdone(1000, 0);
	Person *per;
	FILE *OUTP = fopen(filename.c_str(), "wb");
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
void Genealogy::instantiate() { // computes an instance of all the genomes
	Person *per, *permo, *perfa;
	for (int i = 0; i < npersons; i++) {
		per = persons[i];
		if (per->idfa < 0 || per->idmo < 0) continue; // only children get new genomes
		permo = persons[per->idmo];
		perfa = persons[per->idfa];
		per->genome = (perfa->genome).mate(permo->genome); // order is fa.mate(mo)
	}
}
void Genealogy::showallgenomes() {
	for (int i = 0; i < npersons; i++) {
		printf("%3d:%s:",i,persons[i]->idstr);
		(persons[i]->genome).showgenome();
	}
}
TrialsOutput Genealogy::runtrials(Int ntry, VecInt &cperid, VecInt &mperid) {
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
TrialsOutput Genealogy::runtrials(Int ntry, char **cpersons, char **mpersons) {
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
TrialsOutput Genealogy::runtrials(Int ntry, char **mpersons) {
	// interface with identical lists (not needed?)
	TrialsOutput ans = runtrials(ntry, mpersons, mpersons);
	ans.iamsymmetric = true; // I don't think this is used
	return ans;
}
TrialsOutput Genealogy::runtrials(Int ntry) {
	//run trials using global mpersons_g
	mpersons_g[npersons_g] = NULL;
	return runtrials(ntry, mpersons_g);
}

void Genealogy::showtrials(TrialsOutput &trials, char **cpersons, char **mpersons, MatDoub *measurements) {
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
void Genealogy::dumptrials(FILE *OUTP, TrialsOutput &trials) {
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

void Genealogy::setresolution(Int totb) {
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