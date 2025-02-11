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
 *  Computes mean for each column of a table in plain ascii-format
 *
 *  Created on: 30.10.14
 *      Author: Alexander Klassmann
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

void write2dArray(int, double**&);

/*
 * output precision
 */
#define DEFAULT_PRECISION 4

using namespace std;

void printUsage() {
	cout << "Usage: colmean [-p]" << endl;
	cout << "Options: -p precision of output [default 4 decimal places]" << endl;
}

int main(int argc, char *argv[]) {
	int p = 0, nrow = 0, size = 0;
	vector<double> firstRow;
	double** m = NULL;
	string line;

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else {
			printUsage();
			return (1);
		}
	}

	if (p <= 0) {
		p = DEFAULT_PRECISION;
	}

	while (getline(cin, line)) {
		line.erase(0, line.find_first_not_of(' '));       //prefixing spaces
		line.erase(line.find_last_not_of(' ') + 1);       //surfixing spaces
		if (line.length() > 0) {
			nrow++;
			stringstream is(line);
			if (nrow == 1) {
				while (!is.eof()) {
					double tmp;
					is >> tmp;
					if (is.fail()) {
						cout << "Failed to read in number in row " << nrow << endl;
						return (1);
					}
					size++;
					firstRow.push_back(tmp);
				}
				m = new double*[size];
				for (int i = 0; i < size; i++) {
					m[i] = new double[size];
				}
				for (int i = 0; i < size; i++) {
					m[0][i] = firstRow[i];
				}
			} else {
				for (int i = 0; i < size; i++) {
					double tmp;
					is >> tmp;
					if (is.fail()) {
						cout << "Failed to read in number in row " << nrow << endl;
						return (1);
					}
					m[(nrow - 1) % size][i] += tmp;
				}
			}
		}
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			m[i][j] /= nrow / size;
		}
	}
	cout << fixed << setprecision(p);
	write2dArray(size, m);
}
