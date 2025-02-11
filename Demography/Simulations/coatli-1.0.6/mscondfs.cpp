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
void computeSpectra(bool** data, int chromCount, int snpCount, int k, bool inversion, int** inSpectra, int** exSpectra,
		int* inSize, int* exSize);
void printSpectra(int** inSpectra, int** exSpectra, int snpCount, int fsSize, int* inSize, int* exSize, int k,
		int relation, bool printBoth, bool printFullOuter, bool printCooc, bool printOnlySegSites, bool fold);

int k = 0;
int relation = ANY;
bool inversion = false;
bool printBoth = false;
bool printFullOuter = false;
bool printCooc = false;
bool printOnlySegSites = false;
bool fold = false;

void printUsage() {
	cout << "Usage: mscondfs k [-nested] [-disjoint] [-both] [-c] [-s] [-f] [-u] [-inv]" << endl;
	cout << " k size of focal mutation" << endl;
	cout << " -n or -a print nested spectrum" << endl;
	cout << " -d or -b print disjoint spectrum" << endl;
	cout << " -both print nested and disjoint spectra separated by '|'" << endl;
	cout << " -c include co-occurring sites" << endl;
	cout << " -s print only strictly nested and strictly disjoint" << endl;
	cout << " -f fold spectra (implies option -s)" << endl;
	cout << " -u print full outer spectrum xi_n-k+1 til xi_n-2 (all zero if no recombination)" << endl;
	cout << " -inv inversion: containing mutations are split and counted in both spectra" << endl;
}

int main(int argc, char* argv[]) {

	if (argc > 1) {
		stringstream ss(argv[1]);
		ss >> k;
	}

	for (int i = 2; i < argc; i++) {
		string s(argv[i]);
		if (s == "-inv") {
			inversion = true;
		} else if (s == "-both") {
			printBoth = true;
		} else if (s == "-n" || s == "-a") {
			relation = NESTED;
		} else if (s == "-d" || s == "-b") {
			relation = DISJOINT;
		} else if (s == "-c") {
			printCooc = true;
		} else if (s == "-f") {
			fold = true;
			printOnlySegSites = true;
		} else if (s == "-s") {
			printOnlySegSites = true;
		} else if (s == "-u") {
			printFullOuter = true;
		} else {
			printUsage();
			return (1);
		}
	}
	if (k < 1 || (printFullOuter && printOnlySegSites) || (printCooc && fold)) {
		printUsage();
		return (1);
	}
	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header) {
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** s) {
	int* inSize;
	int* exSize;
	int** inSpectra;
	int** exSpectra;

	if (sampleSize > 0 && segsites > 0) {
		inSize = new int[segsites];
		exSize = new int[segsites];
		inSpectra = new int*[segsites];
		exSpectra = new int*[segsites];
		for (int i = 0; i < segsites; i++) {
			inSpectra[i] = new int[sampleSize - 1];
			exSpectra[i] = new int[sampleSize - 1];
		}

		computeSpectra(s, sampleSize, segsites, k, inversion, inSpectra, exSpectra, inSize, exSize);
		printSpectra(inSpectra, exSpectra, segsites, sampleSize, inSize, exSize, k, relation, printBoth, printFullOuter,
				printCooc, printOnlySegSites, fold);

		for (int i = 0; i < segsites; i++) {
			delete[] inSpectra[i];
			delete[] exSpectra[i];
		}

		delete[] inSpectra;
		delete[] exSpectra;
		delete[] inSize;
		delete[] exSize;
	}
}
