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
 *  Created on: 01.02.2015
 *      Author: Alexander Klassmann
 */

#include <algorithm>
#include <iostream>
#include "ntdef.h"

#define AA 0
#define DA 1

using namespace std;

bool nested(bool** data, int chromCount, int snp1, int snp2) {
	bool cooc = false;
	bool inner1 = false;
	bool inner2 = false;
	for (int i = 0; i < chromCount; i++) {
		if (data[i][snp1] == DA && data[i][snp2] == DA) {
			cooc = true;
		} else if (data[i][snp1] == DA && data[i][snp2] != DA) {
			inner1 = true;
		} else if (data[i][snp1] != DA && data[i][snp2] == DA) {
			inner2 = true;
		}
	}
	return cooc && !(inner1 && inner2);
}

bool disjoint(bool** data, int chromCount, int snp1, int snp2) {
	bool cooc = false;
	for (int i = 0; i < chromCount; i++) {
		if (data[i][snp1] == DA && data[i][snp2] == DA) {
			cooc = true;
		}
	}
	return !cooc;
}

void compute2loci(bool** data, int chromCount, int snpCount, int** twoLoci, int relation) {
	int *size = new int[snpCount];

	for (int snp = 0; snp < snpCount; snp++) {
		size[snp] = 0;
		for (int chrom = 0; chrom < chromCount; chrom++) {
			if (data[chrom][snp] == DA) {
				size[snp]++;
			}
		}
	}

	for (int snp1 = 0; snp1 < snpCount; snp1++) {
		for (int snp2 = 0; snp2 < snp1; snp2++) {
			if (relation == ANY || (relation == NESTED && nested(data, chromCount, snp1, snp2))
					|| (relation == DISJOINT && disjoint(data, chromCount, snp1, snp2))) {
				twoLoci[size[snp1] - 1][size[snp2] - 1]++;
				twoLoci[size[snp2] - 1][size[snp1] - 1]++;
			}
		}
	}

	for (int i = 0; i < chromCount - 1; i++) {
		twoLoci[i][i] /= 2;
	}

	delete[] size;
}

void computeSpectra(bool** data, int chromCount, int snpCount, int k, bool inversion, int** inSpectra, int** exSpectra,
		int* inSize, int* exSize) {

	for (int snp = 0; snp < snpCount; snp++) {
		inSize[snp] = 0;
		exSize[snp] = 0;
		for (int j = 0; j < chromCount; j++) {
			inSpectra[snp][j] = 0;
			exSpectra[snp][j] = 0;
		}
	}

	for (int snp = 0; snp < snpCount; snp++) {
		for (int chrom = 0; chrom < chromCount; chrom++) {
			if (data[chrom][snp] == DA) {
				inSize[snp]++;
			} else if (data[chrom][snp] == AA) {
				exSize[snp]++;
			}
		}
	}

	for (int focalSNP = 0; focalSNP < snpCount; focalSNP++) {
		if (inSize[focalSNP] == k) {
			int inSnps = 0;
			for (int snp = 0; snp <= snpCount - 1; snp++) {
				if (snp != focalSNP) {
					int inMutationCount = 0;
					int exMutationCount = 0;
					for (int chrom = 0; chrom < chromCount; chrom++) {
						if (data[chrom][snp] == DA) {
							if (data[chrom][focalSNP] == DA) {
								inMutationCount++;
							} else if (data[chrom][focalSNP] == AA) {
								exMutationCount++;
							}
						}
					}

					if (inMutationCount / (double) inSize[focalSNP] > exMutationCount / (double) exSize[focalSNP]) {
						if (inversion) {
							if (exMutationCount > 0) {
								exSpectra[focalSNP][exMutationCount - 1]++;
							} else {
								inSpectra[focalSNP][inMutationCount - 1]++;
							}
						} else {
							inSpectra[focalSNP][inMutationCount + exMutationCount - 1]++;
						}
						if (inMutationCount < inSize[focalSNP]) {
							inSnps++;
						}
					} else if (exMutationCount / (double) exSize[focalSNP]
							> inMutationCount / (double) inSize[focalSNP]) {
						exSpectra[focalSNP][exMutationCount + inMutationCount - 1]++;

						if (exMutationCount < exSize[focalSNP]) {
							inSnps++;
						}
					}
				}
			}
		}
	}
}

