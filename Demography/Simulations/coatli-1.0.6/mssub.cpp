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
 *  Created on: 31.07.2020
 *      Author: Alexander Klassmann
 *
 *  Subsets the chromosomes of an ms file and re-exports to ms format
 *
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

void scan();

const int default_n = 5;
const int default_precision = 4;
int n = default_n;
int precision = default_precision;
int *subsample_chromosomes;

void processHeader(const int sampleSize, const int samples, const char *header) {
	if (n > sampleSize) {
		cerr << "Error: ms file contains only " << sampleSize << " chromosomes." << endl;
		exit (EXIT_FAILURE);
	} else {
		string line(header);
		int i = line.find(' ') + 1;
		int j = line.find(' ', i);
		cout << line.substr(0, i) << n << line.substr(j) << endl;

		subsample_chromosomes = new int[n];
		for (int i = 0; i < n; i++) {
			subsample_chromosomes[i] = i;
		}
	}
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double *positions, bool** polymorphisms) {
	bool *polymorphic_in_subsample = new bool[segsites];
	int number_polymorphic = 0;

	for (int i = 0; i < segsites; i++) {
		polymorphic_in_subsample[i] = false;
		for (int j = 1; j < n; j++) {
			if (polymorphisms[subsample_chromosomes[j]][i] != polymorphisms[subsample_chromosomes[0]][i]) {
				polymorphic_in_subsample[i] = true;
				number_polymorphic++;
				break;
			}
		}
	}

	cout << endl << "//" << endl;

	cout << "segsites: " << number_polymorphic << endl;
	if (number_polymorphic > 0) {
		cout << "positions:" << std::setprecision(precision) << fixed;

		for (int i = 0; i < segsites; ++i) {
			if (polymorphic_in_subsample[i]) {
				cout << " " << positions[i];
			}
		}

		cout << " " << endl;

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < segsites; j++) {
				if (polymorphic_in_subsample[j]) {
					cout << polymorphisms[i][j];
				}
			}
			cout << endl;
		}
	} else {
		cout << endl;
	}
	free(polymorphic_in_subsample);

}

void printUsage() {
	cout << "Usage: mssub n" << endl;
	cout << "Subsets chromosomes in ms output, retaining the first n." << endl;
	cout << "Options: n subsample size (number of chromosomes to retain) (default " << default_n << ")" << endl;
	cout << "  -p [integer] precision (number of decimal places) of position output (default " << default_precision
			<< ")" << endl << endl;
	cout << "Example: ms 10 1 -t 5 | mssub" << endl;
	cout << "  retains the first " << default_n << "chromosomes out of 10" << endl;
}

int main(int argc, char *argv[]) {
	for (int i = 1; i < argc; i++) {
		string s = string(argv[i]);
		if (i == 1 && s != "-h") {
			stringstream ss(argv[i]);
			ss >> n;
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

	if (n < 2) {
		cout << "Error: sub sample size cannot be smaller than 2." << endl;
		exit (EXIT_FAILURE);
	}

	if (precision < 1) {
		cout << "Error: precision cannot be smaller than 1." << endl;
		exit (EXIT_FAILURE);
	}

	std::ios_base::sync_with_stdio(true);

	srand (time(NULL));

scan	();

	free(subsample_chromosomes);

	return 0;

}
