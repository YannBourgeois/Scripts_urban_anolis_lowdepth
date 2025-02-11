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
 *  along with coatli.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Calculates the frequency spectrum from output of 'ms'
 *  conditional on a (possibly multiple) existence of a mutation of size k.
 *
 *  Created on: 01.02.15
 *      Author: Alexander Klassmann
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include "ntdef.h"

using namespace std;

void scan();
void compute2loci(bool** data, int chromCount, int snpCount, int** twoLoci, int relation);

int relation = ANY;

void printUsage() {
	cout << "Usage: ms2loci [-nested] [-disjoint]" << endl;
	cout << " -nested print only spectrum of nested mutations" << endl;
	cout << " -disjoint print only spectrum of disjoint mutations" << endl;
}

int main(int argc, char* argv[]) {

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "-nested") {
			relation = NESTED;
		} else if (string(argv[i]) == "-disjoint") {
			relation = DISJOINT;
		} else {
			printUsage();
			return (1);
		}
	}

	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header) {
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** s) {
	int** twoLoci = new int*[sampleSize - 1];
	for (int i = 0; i < sampleSize - 1; i++) {
		twoLoci[i] = new int[sampleSize - 1];
		for (int j = 0; j < sampleSize - 1; j++) {
			twoLoci[i][j] = 0;
		}
	}

	if (sampleSize > 0 && segsites > 0) {
		compute2loci(s, sampleSize, segsites, twoLoci, relation);
	}

	for (int i = 0; i < sampleSize - 1; i++) {
		for (int j = 0; j < sampleSize - 1; j++) {
			cout << twoLoci[i][j];
			if (j < sampleSize - 2) {
				cout << " ";
			}
		}
		cout << endl;
	}

	for (int i = 0; i < sampleSize - 1; i++) {
		delete[] twoLoci[i];
	}
	delete[] twoLoci;
}
