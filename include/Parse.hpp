#ifndef PARSE_H
#define PARSE_H

#include "TrialsOutput.hpp"
#include "DataSet.hpp"
#include "Genealogy.hpp"
#include "MultivariateModel.hpp"
#include "globals.hpp"
#include "myThrow.hpp"

#include <cctype>
#include <cstddef>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <unistd.h>
#define _getcwd getcwd

enum COMMANDS {
	newperson,
	newchild,
	makeged,
	personofinterest,
	runtrials,
	allzeros,
	newdata,
	showmodel,
	showmodelomitting,
	showmodelmoving,
	permanentlyremove,
	seterror,
	done,
	comment,
};

struct Parse {
	private:
		Genealogy gg;
		DataSet dd;
		TrialsOutput trials;
		double measerr = measerr_g;
		map<string, COMMANDS> map_commands;
		bool allow_output = true;

		// handler for 'newperson' command
		int new_person(vector<string> tokens) {
			gg.newPerson(tokens[1].c_str());	
			return 0;
		}

		// handler for 'newchild' command
		int new_child(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing child name in newChild command");
			if (tokens.size() < 3) mythrow("missing first parent name in newChild command");
			if (tokens.size() < 4) {
				gg.newChild(tokens[1].c_str(), tokens[2].c_str());
				if (allow_output) {
					printf("tree: %s is a child of parent %s and anonymous\n", tokens[1].c_str(), tokens[2].c_str());
				}
				return 0; // done, no third token
			}
			if (isdigit(tokens[3][0])) {
				int ngen;
				sscanf(tokens[3].c_str(), "%d", &ngen);
				if (ngen < 2 || ngen > 20) mythrow("bad number of generations in newChild command");
				gg.newChild(tokens[1].c_str(), tokens[2].c_str(), ngen);
				if (allow_output) {
					printf("tree: %s is a %d generation descendent of %s and anonymous\n", tokens[1].c_str(), ngen, tokens[2].c_str());
				}
			}
			else {
				gg.newChild(tokens[1].c_str(), tokens[2].c_str(), tokens[3].c_str());
				if (allow_output) {
					printf("tree: %s is a child of parents %s and %s\n", tokens[1].c_str(), tokens[2].c_str(), tokens[3].c_str());
				}
			}
			return 0;
		}

		// handler for 'personsofinterest'
		// personOfInterest name1 [, name2 ...] <-- will have data or be an hypothesis
		// mark each name as a person of interest (mpersons_g)
		int person_of_interest(vector<string> tokens) {
			for (int i=1; i < tokens.size(); i++) {
				cout << tokens[i] << endl;
				gg.mpersons_g[gg.npersons_g] = new char[tokens[i].length() + 1]; // will ignore memory leak
				strcpy(gg.mpersons_g[gg.npersons_g], tokens[i].c_str());
				if (allow_output) {
					printf("pre-trials: %s is a person of interest\n", gg.mpersons_g[gg.npersons_g]);
				}
				++gg.npersons_g;
			}
			return 0;
		}

		// handler for 'allzeros'
		int all_zeros(vector<string> tokens) {
			for (int i = 0; i < gg.npersons_g; i++) {
				for (int j = i+1; j < gg.npersons_g; j++) {
					dd.addpair(gg.mpersons_g[i], gg.mpersons_g[j], 0.);
				}
			}
			return 0;
		}

		// export tree to GED file
		int make_ged(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing filename in makeGED command");
			if (allow_output) {
				printf("making GED file %s\n", tokens[1].c_str());
			}
			gg.makeGED(tokens[1]);
			return 0;
		}

		// run trials handler
		int run_trials(vector<string> tokens) {
			int ntrials;
			if (tokens.size() < 2) mythrow("missing number of trials in runTrials command");
			sscanf(tokens[1].c_str(), "%d", &ntrials);
			if (ntrials < 1 || ntrials > 100000) mythrow("number of trials too big or too small");
			if (allow_output) {
				printf("starting trials...");
			}
			trials = gg.runtrials(ntrials);
			if (allow_output) {
				printf("done.\n\n");
			}
			return 0;
		}
		
		int new_data(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing first person in newData command");
			if (tokens.size() < 3) mythrow("missing second person in newData command");
			if (tokens.size() < 4) mythrow("missing value in newData command");
			double val;
			int nconv = sscanf(tokens[3].c_str(), "%lf", &val);
			if (nconv != 1 || val < 0. || val > 100000.) mythrow("bad value in newData command");
			dd.addpair(tokens[1].c_str(), tokens[2].c_str(), val);
			if (allow_output) {
				printf("data: %6.1f cM is common between %s and %s\n", val, tokens[1].c_str(), tokens[2].c_str());
			}
			return 0;
		}

		int set_error(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing value in setError comman");
			sscanf(tokens[1].c_str(), "%lf", &measerr);
			if (allow_output) {
				printf("\nchanging error parameter to %.3f (usual value is 0.03)\n", measerr);
			}
			return 0;
		}

		int show_model(vector<string> tokens) {
			if (allow_output) {
				printf("\nBaseline model:\n");
			}
			MultivariateModel mm(gg, trials, dd, measerr);
			mm.logprob_g = mm.logprob;
			mm.chsqprob_g = mm.chisqprob;
			mm.showmodel();
			return 0;
		}

		int show_model_omitting(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing person in showModelOmitting command");
			DataSet d1(dd); // new dataset is copy of dd
			if (allow_output) {
				printf("\n");
			}
			for (int i=1; i < tokens.size(); i++) {
				d1.omit(tokens[i].c_str());
				if (allow_output) {
					printf("Omitting person %s:\n", tokens[i].c_str());
				}
			}
			MultivariateModel m1(gg, trials, d1, measerr);
			m1.showmodel();
			return 0;
		}

		int permanently_remove(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing person in permanentlyRemove command");
			if (allow_output) {
				printf("\n");
			}
			for (int i=1; i < tokens.size(); i++) {
				dd.omit(tokens[i].c_str());
				if (allow_output) {
					printf("Removing person %s.\n", tokens[i].c_str());
				}
			}
			return 0;
		}
		
		int show_model_moving(vector<string> tokens) {
			if (tokens.size() < 2) mythrow("missing first person in showModelMoving command");
			if (tokens.size() < 3) mythrow("missing second person in showModelMoving command");
			DataSet d2(dd); // new dataset is copy of dd
			d2.moveto(tokens[1].c_str(), tokens[2].c_str());
			if (allow_output) {
				printf("\nMoving person %s to %s:\n", tokens[1].c_str(), tokens[2].c_str());
			}
			MultivariateModel m2(gg, trials, d2, measerr);
			m2.showmodel();
			return 0;
		}

		// helper to cast string to all lower
		string to_lower(string val) {
			string res = val;
			for (int i=0; i < res.size(); i++) {
				res[i] = tolower(res[i]);
			}
			return res;
		}
	
		// helper to split line from input file
		vector<string> split_line(string line, string delim = " ,\t") {
			vector<string> res;
			std::size_t start = line.find_first_not_of(delim);
			std::size_t end = line.find_first_of(delim, start);
			// while a valid substring start point
			while (start != string::npos) {
				// valid substring found
				if (end > start || (start != string::npos && end == string::npos)) {
					string ss = line.substr(start, end - start);
					res.push_back(ss);
					// if at end of string exit
					if (end == string::npos) {
						return res;
					}
				}
				// advance start and end
				start = line.find_first_not_of(delim, end);
				end = line.find_first_of(delim, start);
			}
			return res;
		}

	public:
		Parse(bool allow_output = true):
			gg(Genealogy()),
			dd(gg),
			trials(0, 0, 0),
			allow_output(allow_output)
		{
			// initialize map_commands 
			map_commands["newperson"] = newperson;
			map_commands["newchild"] = newchild;
			map_commands["makeged"] = makeged;
			map_commands["personofinterest"] = personofinterest;
			map_commands["runtrials"] = runtrials;
			map_commands["allzeros"] = allzeros;
			map_commands["newdata"] = newdata;
			map_commands["showmodel"] = showmodel;
			map_commands["showmodelomitting"] = showmodelomitting;
			map_commands["showmodelmoving"] = showmodelmoving;
			map_commands["permanentlyremove"] = permanentlyremove;
			map_commands["seterror"] = seterror;
			map_commands["done"] = done;
			map_commands["#"] = comment;
		};

		void set_allow_output(bool value) {
			allow_output = value;
		}

		MultivariateModel get_multivariate_model() {
			return MultivariateModel(gg, trials, dd, measerr);
		}

		int parse_line(string line) {
			vector<string> tokens = split_line(line); // split line
			if (tokens.size() == 0) return -1; // make sure line is not empty
			// switch on first token
			switch (map_commands[to_lower(tokens[0])]) {
				case newperson:
					return new_person(tokens);
				case newchild:
					return new_child(tokens);
				case makeged:
					return make_ged(tokens);
				case personofinterest:
					return person_of_interest(tokens);
				case runtrials:
					return run_trials(tokens);
				case allzeros:
					return all_zeros(tokens);
				case newdata:
					return new_data(tokens);
				case showmodel:
					return show_model(tokens);
				case showmodelomitting:
					return show_model_omitting(tokens);
				case showmodelmoving:
					return show_model_moving(tokens);
				case permanentlyremove:
					return permanently_remove(tokens);
				case seterror:
					return set_error(tokens);
				case done:
					return -2;
				case comment:
					return 0;
				default:
					if (allow_output) {
						printf("with input line \"%s\"\n", line.c_str());
						printf("unrecognized command in input\n");
					}
					return -1;
			}
		}

		int parse_file(string filename) {
			ifstream INP(filename); // open file
			if (!INP.is_open() && allow_output) {
				char cwd[1024];
				_getcwd(cwd, 1023);
				printf("Can\'t find input file %s in working directory %s\n", filename.c_str(), cwd);
				print_help();
			}
			string line; // variable to store line from file as a string
			// iterate over file lines
			while (getline(INP, line)) {
				parse_line(line);
			}
			INP.close(); // close file
			return 0; // all good
		}

		static void print_help() {
			printf("\
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
");
		}
};

#endif