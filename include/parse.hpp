#ifndef PARSE_H
#define PARSE_H

#include "includeLibs.hpp"
#include "DataSet.hpp"
#include "Genealogy.hpp"
#include "MultivariateModel.hpp"
#include "globals.hpp"
#include "myThrow.hpp"

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

#endif