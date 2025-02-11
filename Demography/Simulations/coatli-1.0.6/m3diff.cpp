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
 *  Copyright 2016 Alexander Klassmann
 */
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

#define defaultPrecision 4

int read3dArray(char*, int&, double***&);
void write3dArray(int, double***&);

void printUsage() {
	cout << "Usage: m2diff file1 file2 [-s preserve sign] [-p precision] [-max printMax] [-r relative difference]"
			<< endl;
}

void calcDiff(double*** m1, double*** m2, int size, bool sign, bool relative) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < size; k++) {
				if (relative) {
					if (m1[i][j][k] != 0 || m2[i][j][k] != 0) {
						m1[i][j][k] = (m1[i][j][k] - m2[i][j][k]) / m1[i][j][k];
					}
				} else {
					m1[i][j][k] -= m2[i][j][k];
				}
				if (!sign && m1[i][j][k] < 0)
					m1[i][j][k] *= -1;
			}
		}
	}
}

int main(int argc, char *argv[]) {
	char *fileName1 = NULL, *fileName2 = NULL;
	bool printMax = false, relative = false, sign = false;
	double ***m1, ***m2;
	int p = defaultPrecision;
	int size1 = 0, size2 = 0;

	for (int i = 1; i < argc; i++) {
		string s(argv[i]);
		if (s[0] != '-') {
			if (fileName1 == NULL) {
				fileName1 = argv[i];
			} else if (fileName2 == NULL) {
				fileName2 = argv[i];
			} else {
				printUsage();
				return 1;
			}
		} else if (s == "-max") {
			printMax = true;
		} else if (s == "-r") {
			relative = true;
		} else if (s == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else if (s == "-s") {
			sign = true;
		} else {
			printUsage();
			return 1;
		}
	}
	read3dArray(fileName1, size1, m1);
	read3dArray(fileName2, size2, m2);

	if (size1 != size2) {
		cout << "Error: first array of size " << size1 << ", second of size " << size2 << endl;
	} else {
		calcDiff(m1, m2, size1, sign, relative);

		cout << fixed << setprecision(p);

		if (printMax) {
			double max = m1[0][0][0];
			for (int i = 0; i < size1; i++) {
				for (int j = 0; j < size1; j++) {
					for (int k = 0; k < size1; k++) {
						if (abs(max) < abs(m1[i][j][k])) {
							max = m1[i][j][k];
						}
					}
				}
			}
			cout << max << endl;
		} else {
			write3dArray(size1, m1);
		}
	}

	for (int i = 0; i < size1; i++) {
		for (int j = 0; j < size1; j++) {
			delete[] m1[i][j];
			delete[] m2[i][j];
		}
		delete[] m1[i];
		delete[] m2[i];
	}
	delete m1;
	delete m2;

	return 0;
}
