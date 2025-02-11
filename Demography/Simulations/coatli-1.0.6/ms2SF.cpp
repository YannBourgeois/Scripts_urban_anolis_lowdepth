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
 *  Translates 'ms'-output into Sweepfinder format
 *
 *  Created on: 30.07.20
 *      Author: Alexander Klassmann
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

void scan();

const int default_length = 1000000;
const int default_run = 1;
const int default_precision = 4;
int length = default_length;
int run = default_run;
int precision = default_precision;
int current_segsite = 0;
int current_run = 0;
int number_removed = 0;
int folded = 0;

void processHeader(const int sampleSize, const int samples, const char *header) {
	if (run > samples) {
		cerr << "Error: ms file contains only " << samples << " runs." << endl;
		exit (EXIT_FAILURE);
	} else {
		cout << "position\tx\tn\tfolded" << endl;
	}
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double *positions, bool **polymorphisms) {
	double lastPosition = -1;
	current_run++;

	if (current_run == run) {
		for (int i = 0; i < segsites; i++) {
			current_segsite++;

			double currentPosition = positions[i] * length;
			if (lastPosition == currentPosition) {
				number_removed++;
				continue;
			}
			lastPosition = currentPosition;

			int x = 0;
			for (int j = 0; j < sampleSize; j++) {
				x += polymorphisms[j][i];
			}
			cout << std::setprecision(precision) << fixed << currentPosition << "\t" << x << "\t" << sampleSize << "\t"
					<< folded << endl;
		}
	}
}

void printUsage() {
	cout << "Usage: ms2SF" << endl;
	cout << "Options:" << endl;
	cout << " -length [integer] length of chromosome (positions of ms are multiplied by this number) (default "
			<< default_length << ")" << endl;
	cout << " -run [integer] 'run' (simulation number) in ms file to be transformed (default " << default_run << ")"
			<< endl;
	cout << " -folded print '1' in last column instead of '0'" << endl;
	cout << " -p [integer] precision (number of decimal places) of position output (default " << default_precision
			<< ")" << endl << endl;
	cout << " -h print options" << endl;
}

int main(int argc, char *argv[]) {
	for (int i = 1; i < argc; i++) {
		string s = string(argv[i]);
		if (s == "-length" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> length;
		} else if (s == "-run" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> run;
		} else if (s == "-folded") {
			folded = 1;
		} else if ((s == "-p" || s == "-precision") && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> precision;
		} else if (s == "-h") {
			printUsage();
			return (0);
		} else {
			cout << "Error: unrecognized option '" << s << "'" << endl;
			printUsage();
			return (1);
		}
	}

	if (length < 1) {
		cout << "Error: chromosome length cannot be smaller than 1." << endl;
		exit (EXIT_FAILURE);
	}

	if (run < 1) {
		cout << "Error: 'run' cannot be smaller than 1." << endl;
		exit (EXIT_FAILURE);
	}

	if (precision < 0) {
		cout << "Error: precision cannot be smaller than 0." << endl;
		exit (EXIT_FAILURE);
	}

	std::ios_base::sync_with_stdio(true);

	scan();

	if (number_removed > 0) {
		cerr << "Warning: removed " << number_removed << " sites yielding a position already used in SF file." << endl;
		cerr << "Please consider increasing chromosome length and/or ms output precision!" << endl;
	}

	return 0;
}
