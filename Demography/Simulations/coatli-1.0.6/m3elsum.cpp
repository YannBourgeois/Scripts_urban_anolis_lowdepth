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

void printUsage() {
	cout << "Usage: m3elsum file [-p precision]" << endl;
	cout << " prints the sum over all array elements" << endl;
}

int main(int argc, char *argv[]) {
	char *fileName = NULL;
	double ***m;
	int p = defaultPrecision;
	int size = 0;

	for (int i = 1; i < argc; i++) {
		string s(argv[i]);
		if (s[0] != '-') {
			if (fileName == NULL) {
				fileName = argv[i];
			} else {
				printUsage();
				return 1;
			}
		} else if (s == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else {
			printUsage();
			return 1;
		}
	}
	read3dArray(fileName, size, m);

	cout << fixed << setprecision(p);

	double sum = 0;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for( int k=0; k<size; k++){
				sum += m[i][j][k];
			}
		}
	}
	cout << sum << endl;

	for (int i = 0; i < size; i++) {
		delete[] m[i];
	}
	delete m;

	return 0;
}
