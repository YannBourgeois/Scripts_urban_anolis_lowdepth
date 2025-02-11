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
 *  prints the frequency spectrum of all samples from ms-output
 *
 *  Copyright 2014 Alexander Klassmann
 */

#include <iostream>

using namespace std;

void scan();

int main() {
	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header){
}

void processTree(const char*){
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** s) {
	int i, j, segsiteSize;

	if (sampleSize < 2)
		return;

	int* fs = new int[sampleSize - 1];

	for (i = 0; i < sampleSize - 1; i++) {
		fs[i] = 0;
	}
	if (segsites > 0) {
		for (j = 0; j < segsites; j++) {
			segsiteSize = 0;
			for (i = 0; i < sampleSize; i++) {
				if (s[i][j]) {
					segsiteSize++;
				}
			}

			if (segsiteSize)
				fs[segsiteSize - 1]++;
		}
	}

	for (i = 0; i < sampleSize - 2; i++) {
		cout << fs[i] << " ";
	}
	cout << fs[sampleSize - 2] << endl;

	delete[] fs;
}

