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
 *  Created on: 08.10.2013
 *      Author: Alexander Klassmann
 *
 *
 *  This program filters ms output.
 *  Whenever there are exactly two nested segregating sites (mutations),
 *  it prints the frequency spectrum.
 *
 */

#include <iostream>

void scan();

using namespace std;

int main() {
	std::ios_base::sync_with_stdio(true);
	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header) {
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** s) {
	unsigned int i, k;
	bool nested;

	if (segsites == 2) {
		for (k = 0; k <= 1; k++) {
			nested = false;
			for (i = 0; i < sampleSize; i++) {
				if (s[i][k] && !s[i][1 - k]) {
					nested = true;
				} else if (!s[i][k] && s[i][1 - k]) {
					nested = false;
					break;
				}
			}
			if (nested) {
				int d1 = 0, d2 = 0;
				for (i = 0; i < sampleSize; i++) {
					if (s[i][k]) {
						if (s[i][1 - k]) {
							d2++;
						} else {
							d1++;
						}
					}
				}

				for (int i = 0; i < sampleSize - 1; i++) {
					if ((i + 1 == d1 + d2) || (i + 1 == d2)) {
						cout << "1";
					} else {
						cout << "0";
					}
					if (i < sampleSize - 2) {
						cout << " ";
					} else {
						cout << endl;
					}
				}
			}
		}
	}
}
