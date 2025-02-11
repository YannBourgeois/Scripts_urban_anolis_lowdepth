/**
 *  This file is part of coatli.
 *
 *  coatli is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  coatli is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTAbLITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with coatli.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  This class computes the third moments of the unfolded frequency spectrum
 *  under the assumption of constant population size.
 *
 *  Copyright 2016 Alexander Klassmann
 */

#include "ntdef.h"

double tau(int i);

void initT2(int n);
void freeT2(int n);
double tau(int n, int i, int j, int relation);
double* getCondE1(int n, int k, double theta, int relation);
double** getE2(int n, double theta, int relation, int thetaExp);

void initT3(int n);
void freeT3(int n);
double tau(int n, int h, int i, int j, int relation);

/**
 * returns E[xi_i xi_j | xi_k],
 */
double** getCondE2(int n, int k, double theta, int relation, int thetaExp, bool includeFocalMutation) {
	initT3(n);
	double** condE2 = new double*[n - 1];

	for (int i = 0; i < n - 1; i++) {
		condE2[i] = new double[n - 1];
		for (int j = 0; j < n - 1; j++) {
			condE2[i][j] = 0;
		}
	}

	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			if (thetaExp == ANY || thetaExp == QUADRATIC) {
				condE2[i][j] += tau(n, k, i + 1, j + 1, relation) * theta * theta;
			}
			if (thetaExp == ANY) {
				if (k - 1 == i && k - 1 == j && relation == NESTED_NESTED) {
					if (includeFocalMutation) {
						condE2[i][i] += 1. / k;
					}
				}
			}
			if (thetaExp == ANY || thetaExp == LINEAR) {
				if (k - 1 == i) {
					if (includeFocalMutation) {
						if (relation == NESTED_NESTED) {
							condE2[i][j] += tau(n, i + 1, j + 1, NESTED) * theta;
						} else if (relation == NESTED_DISJOINT) {
							condE2[i][j] += tau(n, i + 1, j + 1, DISJOINT) * theta;
						}
					}
				}
				if (k - 1 == j) {
					if (includeFocalMutation) {
						if (relation == NESTED_NESTED) {
							condE2[i][j] += tau(n, i + 1, j + 1, NESTED) * theta;
						} else if (relation == DISJOINT_NESTED) {
							condE2[i][j] += tau(n, i + 1, j + 1, DISJOINT) * theta;
						}
					}
				}
				if (i == j) {
					if (relation == ANY || relation == NESTED_NESTED) {
						condE2[i][j] += tau(n, k, j + 1, NESTED) * theta;
					}
					if (relation == ANY || relation == DISJOINT_DISJOINT) {
						condE2[i][j] += tau(n, k, j + 1, DISJOINT) * theta;
					}
				}
			}
			condE2[i][j] *= k;
		}
	}
	freeT3(n);
	return condE2;
}

double** getCondE2(int n, int k, double theta, int relation, int thetaExp) {
	return getCondE2(n, k, theta, relation, thetaExp, false);
}

/**
 * returns E[xi_i xi_j | xi_k] - E[xi_i | xi_k] E[xi_j | xi_k]
 */
double** getCondCov(int n, int k, double theta, int relation, int thetaExp) {
	double** condCov = getCondE2(n, k, theta, relation, thetaExp);
	double *condXi1, *condXi2;

	if (thetaExp == ANY || thetaExp == QUADRATIC) {

		int relation1 = ANY, relation2 = ANY;
		if (relation == NESTED_NESTED) {
			relation1 = relation2 = NESTED;
		} else if (relation == DISJOINT_DISJOINT) {
			relation1 = relation2 = DISJOINT;
		} else if (relation == NESTED_DISJOINT) {
			relation1 = NESTED;
			relation2 = DISJOINT;
		} else if (relation == DISJOINT_NESTED) {
			relation1 = DISJOINT;
			relation2 = NESTED;
		}

		condXi1 = getCondE1(n, k, theta, relation1);
		condXi2 = getCondE1(n, k, theta, relation2);

		for (int i = 0; i < n - 1; i++) {
			for (int j = 0; j < n - 1; j++) {
				condCov[i][j] -= condXi1[i] * condXi2[j];
			}
		}

		delete[] condXi1;
		delete[] condXi2;
	}

	return condCov;
}

double** getCondCov(int n, int k, double theta, int relation) {
	return getCondCov(n, k, theta, relation, false);
}

/**
 * returns one element of the 3D array of third moments E[xi_h xi_i xi_j]
 */
double getE3element(int n, double theta, int thetaExp, bool central, int h, int i, int j) {
	initT3(n);
	double e3 = 0;
	if (thetaExp == ANY || thetaExp == LINEAR) {
		if (h == i && i == j) {
			e3 += tau(h) * theta;
		}
	}

	if (thetaExp == ANY || thetaExp == QUADRATIC) {
		if (h == i) {
			e3 += tau(n, i, j, ANY) * theta * theta;
		}
		if (i == j) {
			e3 += tau(n, i, h, ANY) * theta * theta;
		}
		if (j == h) {
			e3 += tau(n, h, i, ANY) * theta * theta;
		}
	}

	if (thetaExp == ANY || thetaExp == CUBIC) {
		e3 += tau(n, h, i, j, ANY) * theta * theta * theta;
	}

	if (central && thetaExp != LINEAR) {
		if (thetaExp == ANY || thetaExp == CUBIC) {
			e3 += 2 * tau(h) * tau(i) * tau(j) * theta * theta * theta;
		}

		double** e2 = getE2(n, theta, ANY, thetaExp == ANY ? ANY : thetaExp - 1);

		e3 -= tau(h) * theta * e2[i - 1][j - 1] + tau(i) * theta * e2[h - 1][j - 1] + tau(j) * theta * e2[h - 1][i - 1];

		for (int i = 0; i < n - 1; i++) {
			delete[] e2[i];
		}
		delete[] e2;
	}

	freeT3(n);

	return e3;
}

/**
 * returns 3D array E[xi_h xi_i xi_j]
 */
double*** getE3(int n, double theta, int thetaExp, bool central) {
	initT3(n);
	double*** e3 = new double**[n - 1];
	for (int i = 0; i < n - 1; i++) {
		e3[i] = new double*[n - 1];
		for (int j = 0; j < n - 1; j++) {
			e3[i][j] = new double[n - 1];
			for (int k = 0; k < n - 1; k++) {
				e3[i][j][k] = 0;
			}
		}
	}

	for (int i = 0; i < n - 1; i++) {

		if (thetaExp == ANY || thetaExp == LINEAR) {
			e3[i][i][i] += tau(i + 1) * theta;
		}

		if (thetaExp == ANY || thetaExp == QUADRATIC) {
			e3[i][i][i] += 3 * tau(n, i + 1, i + 1, ANY) * theta * theta;
		}

		if (thetaExp == ANY || thetaExp == CUBIC) {
			e3[i][i][i] += tau(n, i + 1, i + 1, i + 1, ANY) * theta * theta * theta;
		}

		for (int j = i + 1; j < n - 1; j++) {

			if (thetaExp == ANY || thetaExp == QUADRATIC) {
				double t2 = tau(n, i + 1, j + 1, ANY) * theta * theta;
				e3[i][i][j] += t2;
				e3[i][j][j] += t2;
			}

			if (thetaExp == ANY || thetaExp == CUBIC) {
				e3[i][i][j] += tau(n, i + 1, i + 1, j + 1, ANY) * theta * theta * theta;
				e3[i][j][j] += tau(n, i + 1, j + 1, j + 1, ANY) * theta * theta * theta;
			}

			e3[i][j][i] = e3[j][i][i] = e3[i][i][j];
			e3[j][i][j] = e3[j][j][i] = e3[i][j][j];

			if (thetaExp == ANY || thetaExp == CUBIC) {
				for (int h = j + 1; h < n - 1; h++) {
					e3[h][i][j] = e3[i][j][h] = e3[j][h][i] = e3[h][j][i] = e3[i][h][j] = e3[j][i][h] = tau(n, h + 1,
							i + 1, j + 1,
							ANY) * theta * theta * theta;
				}
			}
		}
	}

	if (central && thetaExp != LINEAR) {
		double** e2 = getE2(n, theta, ANY, thetaExp == ANY ? ANY : thetaExp - 1);
		for (int h = 0; h < n - 1; h++) {
			for (int i = 0; i < n - 1; i++) {
				for (int j = 0; j < n - 1; j++) {
					if (thetaExp == CUBIC) {
						e3[h][i][j] += 2 * tau(h + 1) * tau(i + 1) * tau(j + 1) * theta * theta * theta;
					}
					e3[h][i][j] -= tau(h + 1) * theta * e2[i][j] + tau(i + 1) * theta * e2[h][j]
							+ tau(j + 1) * theta * e2[h][i];
				}
			}
		}
		for (int i = 0; i < n - 1; i++) {
			delete[] e2[i];
		}
		delete[] e2;
	}

	freeT3(n);

	return e3;
}

double*** getE3(int n, double theta) {
	return getE3(n, theta, ANY, false);
}

/**
 * folds the spectrum
 */
void fold(int n, double*** &unfolded) {
	int size = n / 2;
	double ***folded = new double**[size];
	for (int i = 0; i < size; i++) {
		folded[i] = new double*[size];
		for (int j = 0; j < size; j++) {
			folded[i][j] = new double[size];
		}
	}

	for (int h = 0; h < size; h++) {
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				folded[h][i][j] = unfolded[h][i][j] + unfolded[h][i][n - 2 - j] + unfolded[h][n - 2 - i][j]
						+ unfolded[h][n - 2 - i][n - 2 - j] + unfolded[n - 2 - h][i][j]
						+ unfolded[n - 2 - h][i][n - 2 - j] + unfolded[n - 2 - h][n - 2 - i][j]
						+ unfolded[n - 2 - h][n - 2 - i][n - 2 - j];
				folded[h][i][j] /= (1 + (h == n - 2 - h ? 1 : 0)) * (1 + (i == n - 2 - i ? 1 : 0))
						* (1 + (j == n - 2 - j ? 1 : 0));
			}
		}
	}

	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			delete[] unfolded[i][j];
		}
		delete[] unfolded[i];
	}
	delete[] unfolded;

	unfolded = folded;
}
