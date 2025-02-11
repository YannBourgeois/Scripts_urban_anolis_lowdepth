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
 *
 *  Prints the moments of the number of segregating sites S
 *  as calculated by Watterson 1975 and Klassmann & Ferretti 2016
 *
 *  Copyright 2016 Alexander Klassmann
 */
#include <iostream>
#include <cstdlib>
#include <iomanip>

/*
 * output precision
 */
#define defaultPrecision 4

using namespace std;

void printUsage() {
	cout << "Usage: Smom n [-p] [-t] [-3]" << endl;
	cout << "n sample size; with k=0, print third moments" << endl;
	cout << "-p precision (number of decimals, default 6)" << endl;
	cout << "-t theta (default 1)" << endl;
	cout << "-3 only cubic terms" << endl;
}

int main(int argc, char *argv[]) {
	bool cubicOnly = false;
	double h1 = 0, h2 = 0, h3 = 0, theta = 1.0;
	int n = 0, p = 0;

	if (argc == 1) {
		printUsage();
		exit(1);
	}

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else if (string(argv[i]) == "-t" && i + 1 < argc) {
			i++;
			theta = atof(argv[i]);
		} else if (string(argv[i]) == "-3") {
			cubicOnly = true;
		} else if (n == 0) {
			n = atoi(argv[i]);
			if (n <= 1) {
				printUsage();
				exit(1);
			}
		} else {
			printUsage();
			exit(1);
		}
	}

	if (theta < 0) {
		printUsage();
		exit(1);
	}

	if (p <= 0)
		p = defaultPrecision;

	for (int i = 1; i < n; i++) {
		if (!cubicOnly) {
			h1 += 1. / i;
			h2 += 1. / (i * i);
		}
		h3 += 1. / (i * i * i);
	}

	cout << fixed << setprecision(p);

	cout << "E[S]:\t" << h1 * theta << endl;
	cout << "E[S^2]:\t" << h1 * theta + (h1 * h1 + h2) * theta * theta << endl;
	cout << "Var[S]:\t" << h1 * theta + h2 * theta * theta << endl;
	cout << "E[S^3]:\t"
			<< h1 * theta + 3 * (h1 * h1 + h2) * theta * theta
					+ (h1 * h1 * h1 + 3 * (h1 * h2) + 2 * h3) * theta * theta * theta << endl;
	cout << "mu3[S]:\t" << h1 * theta + 3 * h2 * theta * theta + 2 * h3 * theta * theta * theta << endl;
	return (0);
}
