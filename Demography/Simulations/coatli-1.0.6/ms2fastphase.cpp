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
 *  Created on: 09.08.2020
 *      Author: Alexander Klassmann
 *
 *  Transforms output of ms to input for fastphase
 *
 */

#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

void scan();

const int default_precision = 4;
const int default_run = 1;
const string default_individualIDbase = "HG";
bool randomize = false;
bool includeIndividualIDs = true;
bool includePositions = false;
int run = default_run;
int precision = default_precision;
int current_run = 0;
string individualIDbase = default_individualIDbase;

void randomizeGenotypes(int sampleSize, int segsites, bool **polymorphisms) {
	for (int i = 0; i < sampleSize; i += 2) {
		for (int j = 0; j < segsites; j++) {
			if (rand() % 2) {
				bool tmp = polymorphisms[i][j];
				polymorphisms[i][j] = polymorphisms[i + 1][j];
				polymorphisms[i + 1][j] = tmp;
			}
		}
	}
}

void processHeader(const int sampleSize, const int samples, const char *header) {
	if (run > samples) {
		cerr << "Error: ms file contains only " << samples << " runs." << endl;
		exit(EXIT_FAILURE);
	}

	if (sampleSize % 2 != 0) {
		cerr << "Error: ms file contains an uneven number of chromosomes." << endl;
		exit(EXIT_FAILURE);
	}
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double *positions, bool **polymorphisms) {

	current_run++;
	if (current_run == run) {

		int number_preceeding_zeros = floor(log10(sampleSize / 2)) + 1;

		if (randomize) {
			randomizeGenotypes(sampleSize, segsites, polymorphisms);
		}

		cout << sampleSize / 2 << endl;
		cout << segsites << endl;

		if (segsites > 0) {
			if (includePositions) {
				cout << "P" << std::setprecision(precision) << fixed;
				for (int i = 0; i < segsites; ++i) {
					cout << " " << positions[i];
				}
				cout << " " << endl;
			}

			for (int i = 0; i < sampleSize; ++i) {
				if (includeIndividualIDs) {
					if (i % 2 == 0) {
						cout << "#" << individualIDbase << setfill('0') << setw(number_preceeding_zeros) << (i / 2 + 1)
								<< endl;
					}
				}
				for (int j = 0; j < segsites; j++) {
					cout << polymorphisms[i][j];
				}
				cout << endl;
			}
		} else {
			cout << endl;
		}
	}
}

void printUsage() {
	cout << "Usage: ms2fastphase" << endl;
	cout << "Transforms output from ms to input for fastPHASE." << endl;
	cout << "The program assumes that each 2 consecutive chromosomes belong to the same diploid individual." << endl;
	cout << "Options: -n do not include IDs for individuals" << endl;
	cout << "  -IDbase [string] use string as base for an individual ID" << endl;
	cout << "  -M include positions (needed only for modeling options -M.. of fastPHASE)" << endl;
	cout << "  -p [integer] precision (number of decimal places) of position output (default " << default_precision
			<< ")" << endl;
	cout << "  -random randomize genotypes" << endl;
	cout << "  -run [integer] 'run' (simulation number) in ms file to be transformed (default " << default_run << ")"
			<< endl;
	cout << "Example: ms 10 1 -t 5 | ms2fastphase" << endl;
}

int main(int argc, char *argv[]) {
	for (int i = 1; i < argc; i++) {
		string s = string(argv[i]);
		if (s == "-random" || s == "-r") {
			randomize = true;
		} else if (s == "-n") {
			includeIndividualIDs = false;
		} else if (s == "-M") {
			includePositions = true;
		} else if (s == "-IDbase" && argc >= i) {
			i++;
			individualIDbase = argv[i];
		} else if ((s == "-p" || s == "-precision") && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> precision;
		} else if (s == "-run" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> run;
		} else if (s == "-h") {
			printUsage();
			return (0);
		} else {
			cout << "Error: unrecognized option '" << s << "'" << endl;
			printUsage();
			return (1);
		}
	}

	if (precision < 1) {
		cout << "Error: precision cannot be smaller than 1." << endl;
		exit(EXIT_FAILURE);
	}

	if (run < 1) {
		cout << "Error: 'run' cannot be smaller than 1." << endl;
		exit(EXIT_FAILURE);
	}

	std::ios_base::sync_with_stdio(true);

	srand(time(NULL));

	scan();

	return 0;

}
