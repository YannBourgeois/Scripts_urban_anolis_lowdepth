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
 * Calculates the third moments from data in columns
 *
 *  Created on: 01.06.15
 *      Author: AlexanderKlassmann
 */
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

void write3dArray(int, double***&);

/*
 * maximal frequency spectrum size
 */
#define max_fs_size 10000

/*
 * output precision
 */
#define defaultPrecision 4

void printHelp() {
	cerr << "Usage: fs3mom [-mu][-skewness][-p precision]" << endl;
	cerr << "Options: " << endl;
	cerr << "-mu print central moments" << endl;
	cerr << "-skewness print skewness (implies -mu)" << endl;
	cerr << "-p precision (digits after decimal point) of output" << endl;
}

int main(int argc, char *argv[]) {
	double fs[max_fs_size];
	int i, j, k, fs_size = 0, rows = 0, p = 0;
	string line;
	bool central = false, skewness = false;
	double *x = NULL;
	double **xy = NULL;
	double ***xyz = NULL;

	for (i = 1; i < argc; i++) {
		if (string(argv[i]) == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else if (string(argv[i]) == "-mu") {
			central = true;
		} else if (string(argv[i]) == "-skewness") {
			skewness = true;
		} else {
			printHelp();
			return (1);
		}
	}

	if (p <= 0)
		p = defaultPrecision;

	while (getline(cin, line)) {
		rows++;
		stringstream is(line);
		fs_size = 0;
		while (!is.eof()) {
			is >> fs[fs_size++];
		}

		if (x == NULL) {
			x = new double[fs_size];

			for (i = 0; i < fs_size; i++) {
				x[i] = 0;
			}

			xy = new double*[fs_size];

			for (i = 0; i < fs_size; i++) {
				xy[i] = new double[fs_size];
				for (j = 0; j < fs_size; j++) {
					xy[i][j] = 0;
				}
			}

			xyz = new double**[fs_size];

			for (i = 0; i < fs_size; i++) {
				xyz[i] = new double*[fs_size];
				for (j = 0; j < fs_size; j++) {
					xyz[i][j] = new double[fs_size];
					for (k = 0; k < fs_size; k++) {
						xyz[i][j][k] = 0;
					}
				}
			}

		}

		for (i = 0; i < fs_size; i++) {
			x[i] += fs[i];
			for (j = 0; j < fs_size; j++) {
				xy[i][j] += fs[i] * fs[j];
				for (k = 0; k < fs_size; k++) {
					xyz[i][j][k] += fs[i] * fs[j] * fs[k];
				}
			}
		}

	}
	/*
	 * average over rows
	 */

	for (i = 0; i < fs_size; i++) {
		x[i] /= rows;
		for (j = 0; j < fs_size; j++) {
			xy[i][j] /= rows;
			for (k = 0; k < fs_size; k++) {
				xyz[i][j][k] /= rows;
			}
		}
	}

	cout << fixed << setprecision(p);

	/* calculate third moments */

	for (i = 0; i < fs_size; i++) {
		for (j = 0; j < fs_size; j++) {
			for (k = 0; k < fs_size; k++) {
				if (skewness || central) {
					xyz[i][j][k] += -x[i] * xy[j][k] - x[j] * xy[i][k] - x[k] * xy[i][j] + 2 * x[i] * x[j] * x[k];
					if (skewness) {
						xyz[i][j][k] /= sqrt(
								(xy[i][j] - x[i] * x[j]) * rows / (rows - 1) * (xy[i][k] - x[i] * x[k]) * rows
										/ (rows - 1) * (xy[j][k] - x[j] * x[k])) * rows / (rows - 1);
					}
				}
			}
		}
	}

	write3dArray(fs_size, xyz);

	for (int i = 0; i < fs_size; i++) {
		for (int j = 0; j < fs_size; j++) {
			delete[] xyz[i][j];
		}
		delete[] xyz[i];
		delete[] xy[i];
	}
	delete[] xyz;
	delete[] xy;
	delete[] x;
	return (0);
}
