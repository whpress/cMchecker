#ifndef DATASET_H
#define DATASET_H

#include "includeLibs.hpp"
#include "Genealogy.hpp"
#include "globals.hpp"
#include "myThrow.hpp"

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

#endif