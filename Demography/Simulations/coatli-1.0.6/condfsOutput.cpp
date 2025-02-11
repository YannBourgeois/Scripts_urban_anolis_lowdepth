/*
 *  This file is part of coatli
 *
 *  tba is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  tba is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with coatli.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: 01.02.15
 *      Author: Alexander Klassmann
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include "ntdef.h"

using namespace std;

void printSpectra(int** inSpectra, int** exSpectra, int snpCount, int fsSize, int* inSize, int* exSize, int k,
		int relation, bool printBoth, bool printFullOuter, bool printCooc, bool printOnlySegSites, bool fold) {

	for (int snp = 0; snp < snpCount; snp++) {
		if (k == inSize[snp]) {
			if (relation == ANY && !printBoth) {
				for (int j = 0; j < fsSize - 1; j++) {
					cout << inSpectra[snp][j] + exSpectra[snp][j];
					if (j < fsSize - 2)
						cout << " ";
				}
			} else {
				if (relation == NESTED || printBoth) {
					int maxj = fsSize - 1;
					if (printOnlySegSites)
						maxj = k - 1;
					if (printCooc) {
						maxj = k;
					}
					if (fold) {
						for (int j = 0; j < maxj / 2; j++) {
							cout << inSpectra[snp][j] + inSpectra[snp][k - 2 - j];
							if (j < maxj / 2 - 1 || maxj % 2 == 1)
								cout << " ";
						}
						if (maxj % 2 == 1) {
							cout << inSpectra[snp][(maxj - 1) / 2];
						}
					} else {
						for (int j = 0; j < maxj; j++) {
							cout << inSpectra[snp][j];
							if (j < maxj - 1)
								cout << " ";
						}
					}
					if (printBoth) {
						cout << " | ";
					}
				}
				if (relation == DISJOINT || printBoth) {
					int maxj = fsSize - 1;
					if (!printFullOuter)
						maxj = fsSize - k;
					if (printOnlySegSites)
						maxj = fsSize - k - 1;
					if (fold) {
						for (int j = 0; j < maxj / 2; j++) {
							cout << exSpectra[snp][j] + exSpectra[snp][fsSize - k - 2 - j];
							if (j < maxj / 2 - 1 || maxj % 2 == 1)
								cout << " ";
						}
						if (maxj % 2 == 1) {
							cout << exSpectra[snp][(maxj - 1) / 2];
						}
					} else {
						for (int j = 0; j < maxj; j++) {
							cout << exSpectra[snp][j];
							if (j < maxj - 1) {
								cout << " ";
							}
						}
					}
				}
			}
			cout << endl;
		}
	}
}
