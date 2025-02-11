/*
 *  This file is part of coatli.
 *
 *  coatli is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  coatli is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with coatli.  If not, see <http://www.gnu.org/licenses/>
 *
 *  Created on: 10.10.2013
 *      Author: Alexander Klassmann
 *
 *  Processes output of ms; (re-)prints (nested/disjoint) mutations
 *  relative to focal mutations of a specified size;
 */

#include <iostream>
#include <cstdio>
#include <iomanip>
#include <vector>
#include <cstdlib>
#include <limits>
#include "ntdef.h"

using namespace std;

void scan();

short int relation = ANY;
int s_min = 0;
int s_max = std::numeric_limits<int>::max();
int k = 0;
bool includeFocal = false;

void printUsage() {
	cout << "Usage: mscond k [-nested][-disjoint][-i]" << endl;
	cout << "Filters output of ms conditional on mutations of size k" << endl;
	cout << "Options: -disjoint prints only disjoint segregating sites" << endl;
	cout << "         -nested prints only nested segregating sites" << endl;
	cout << "         -i include focal mutation itself" << endl;
	cout << "         -smin <int> include only samples with at least smin segregating sites" << endl;
	cout << "         -smax <int> include only samples with at most smax segregating sites" << endl;
}

int main(int argc, char *argv[]) {

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "-disjoint") {
			relation = DISJOINT;
		} else if (string(argv[i]) == "-nested") {
			relation = NESTED;
		} else if (string(argv[i]) == "-i") {
			includeFocal = true;
		} else if (string(argv[i]) == "-smin" && i + 1 < argc) {
			i++;
			s_min = atoi(argv[i]);
		} else if (string(argv[i]) == "-smax" && i + 1 < argc) {
			i++;
			s_max = atoi(argv[i]);
		} else {
			k = atoi(argv[i]);
		}
	}

	if (k < 1) {
		printUsage();
		return (1);
	}
	std::ios_base::sync_with_stdio(true);
	cout << fixed << setprecision(4);

	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header) {
	cout << header << endl;
}

void processTree(const char* newick) {
	cout << newick << endl;
}

/*
 * adds mutations of size k to vector
 * returns true if at least one is found
 */
bool extractFocalMutations(int n, int s, bool** pol, vector<int>& focalMutations) {
	int mutationSize;
	focalMutations.clear();
	for (int j = 0; j < s; j++) {
		mutationSize = 0;
		for (int i = 0; i < n; i++) {
			if (pol[i][j]) {
				mutationSize++;
			}
		}

		if (mutationSize == k) {
			focalMutations.push_back(j);
		}
	}
	return focalMutations.size() > 0;
}

void printSubsample(const int n, vector<int>& remainingSites, const double* pos, bool** pol) {
	cout << "segsites: " << remainingSites.size() << endl;

	if (remainingSites.size() > 0) {
		cout << "positions:";

		for (int j = 0; j < remainingSites.size(); j++) {
			cout << " " << pos[remainingSites.at(j)];
		}
		cout << endl;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < remainingSites.size(); j++) {
				cout << pol[i][remainingSites.at(j)];
			}
			cout << endl;
		}
	}
}

void constructSubsample(const int n, const int s, const double* pos, bool** pol, int focalMutation) {
	vector<int> remainingSites;

	for (int j = 0; j < s; j++) {
		if (includeFocal || j != focalMutation) {
			if (relation == ANY) {
				remainingSites.push_back(j);
			} else {
				int nestedsize = 0;
				int disjointsize = 0;
				for (int i = 0; i < n; i++) {
					if (pol[i][j]) {
						if (pol[i][focalMutation]) {
							nestedsize++;
						} else {
							disjointsize++;
						}
					}
				}
				// this allows small overlapping between nested and disjoint mutations, e.g. due to recombination
				if (relation == NESTED && nestedsize * (n - k) > disjointsize * k) {
					remainingSites.push_back(j); //sharedSize/k > extSize/(n-k)
				} else if (relation == DISJOINT && nestedsize * (n - k) < disjointsize * k) {
					remainingSites.push_back(j);
				}
			}
		}
	}
	printSubsample(n, remainingSites, pos, pol);
}

void processSample(const int n, const int s, const double* pos, bool** pol) {
	vector<int> focalMutations;
	if (s >= s_min && s <= s_max) {
		if (extractFocalMutations(n, s, pol, focalMutations)) {
			for (int i = 0; i < focalMutations.size(); i++) {
				cout << endl << "// ";
				if (relation == ANY) {
					if (!includeFocal) {
						cout << "original ";
					}
					cout << "focal mutation " << (focalMutations.at(i) + 1) << endl;
				} else if (relation == NESTED) {
					if (includeFocal) {
						cout << "focal mutation " << (i + 1) << ", ";
					}
					cout << "original focal mutation " << (focalMutations.at(i) + 1) << endl;
				} else {
					cout << "original focal mutation " << (focalMutations.at(i) + 1) << endl;
				}
				constructSubsample(n, s, pos, pol, focalMutations.at(i));
			}
		}
	}
}
