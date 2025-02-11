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
 *  This class computes the second moments of the unfolded frequency spectrum
 *  under the assumption of constant population size.
 *
 *  Copyright 2016 Alexander Klassmann
 */

#include "ntdef.h"

double tau(int i);
using namespace std;

namespace {
int lastN = 0;
double* a = 0;

/**
 * Harmonic numbers are calculated only once;
 * since the array is small, it is not necessary to delete it (although possible with freeT2)
 */
void initA(int n) {
	if (lastN < n) {
		if (a != 0) {
			delete[] a;
		}
		a = new double[n + 1];
		a[0] = 0;
		for (int i = 1; i <= n; i++) {
			a[i] = a[i - 1] + 1. / i;
		}
		lastN = n;
	}
}

double beta(int n, int i) {
	if (i > 0 && i < n) {
		return 2 * n * (a[n] - a[i - 1]) / ((n - i + 1) * (n - i)) - 2. / (n - i);
	}

	return 0.;
}

double s_2(int n, int i) {
	return 2 * (a[n - 1] - a[i - 1]) / (n - i) - beta(n, i);
}

double s(int n, int i, int j) {
	if (i + j < n) {
		return beta(n, i + j - 1) - beta(n, i + j);
	} else if (i + j == n) {
		return 1. / (n - 1);
	}
	return 0.;
}

double sa(int n, int i, int j) {
	if (i > j) {
		return (beta(n, j) - beta(n, j + 1)) / 2;
	} else if (i == j) {
		return beta(n, i) - (a[n - 1] - a[i - 1]) / (n - i);
	}
	return 0.;
}

double sb(int n, int i, int j) {
	if (i + j < n) {
		return 1. / (i * j) - 1. / (i * (i + j))
				- (beta(n, j) - beta(n, j + 1) + beta(n, i + j - 1) - beta(n, i + j)) / 2;
	} else if (i + j == n) {
		return (a[n - 1] - a[j - 1]) / (n - j) - beta(n, j) / 2 - 1. / (2 * (n - 1));
	}
	return 0.;
}

double ta(int n, int i, int j) {
	double s = sa(n, i, j);
	if (i == j) {
		s += s_2(n, i) / 2;
	}
	return s;
}

double tb(int n, int i, int j) {
	return s(n, i, j) / 2 + sb(n, i, j);
}
}

/**
 * clear Harmonic numbers (not really necessary, since the array is small)
 */
void freeT2(int n) {
	delete[] a;
	a = 0;
}

/**
 *
 */
double tau(int n, int i, int j, int relation) {
	initA(n);

	double s = 0;
	if (relation == ANY || relation == NESTED) {
		s += ta(n, i, j) + ta(n, j, i);
	}
	if (relation == ANY || relation == DISJOINT) {
		s += tb(n, i, j) + tb(n, j, i);
	}

	return s;
}

/**
 * returns the conditional frequency spectrum, given a mutation of size k
 * and the indicated relation between the two lines (nested/disjoint),
 * if the focalMutation is included, the conditional frequency of nested mutations
 * of size k is augmented by 1 (the focal mutation)
 */
double* getCondE1(int n, int k, double theta, int relation, bool includeFocalMutation) {
	double* condE1 = new double[n - 1];

	initA(n);
	for (int i = 0; i < n - 1; i++) {
		condE1[i] = tau(n, i + 1, k, relation) * theta;

		if (relation == ANY || relation == NESTED) {
			if (includeFocalMutation && i + 1 == k) {
				condE1[i] += tau(k);
			}
		}

		condE1[i] *= k;
	}
	return condE1;
}

/**
 * returns the conditional frequency spectrum, given a mutation of size k
 * and the indicated relation between the two mutations (nested/disjoint)
 */
double* getCondE1(int n, int k, double theta, int relation) {
	return getCondE1(n, k, theta, relation, false);
}

/**
 * returns E[xi_i xi_j], with the given relation between the mutations (nested or disjoint),
 * by indication only the terms quadratic in theta
 */
double** getE2(int n, double theta, int relation, int thetaExp) {
	double **e2 = new double*[n - 1];

	for (int i = 0; i < n - 1; i++) {
		e2[i] = new double[n - 1];
		for (int j = 0; j < n - 1; j++) {
			e2[i][j] = 0;
		}
	}

	initA(n);

	for (int i = 0; i < n - 1; i++) {
		if (thetaExp == ANY || thetaExp == QUADRATIC) {
			for (int j = 0; j < n - 1; j++) {
				e2[i][j] += tau(n, i + 1, j + 1, relation) * theta * theta;
			}
		}
		if (thetaExp == ANY || thetaExp == LINEAR) {
			if (relation == ANY || relation == NESTED) {
				e2[i][i] += tau(i + 1) * theta;
			}
		}
	}
	return e2;
}

/**
 * returns E[xi_i xi_j], restricted to the given relation
 * between the two mutations (nested or disjoint)
 */
double** getE2(int n, double theta, int relation) {
	return getE2(n, theta, relation, ANY);
}

/**
 * returns E[xi_i xi_j]
 */
double** getE2(int n, double theta) {
	return getE2(n, theta, ANY, ANY);
}

/**
 * returns E[xi_i xi_j]-E[xi_i]E[xi_j]
 */
double** getCentralMom(int n, double theta, int thetaExp) {

	double **mom = new double*[n - 1];

	for (int i = 0; i < n - 1; i++) {
		mom[i] = new double[n - 1];
		for (int j = 0; j < n - 1; j++) {
			mom[i][j] = 0;
		}
	}

	initA(n);

	for (int i = 0; i < n - 1; i++) {
		if (thetaExp == ANY || thetaExp == QUADRATIC) {
			for (int j = 0; j < n - 1; j++) {
				mom[i][j] += (tau(n, i + 1, j + 1, ANY) - tau(i + 1) * tau(j + 1)) * theta * theta;
			}
		}
		if (thetaExp == ANY || thetaExp == LINEAR) {
			mom[i][i] += tau(i + 1) * theta;
		}
	}

	return mom;
}

/**
 * returns quadratic part of E[xi_i xi-j]-E[xi_i]E[xi_j],
 * corresponds to the sigma defined by Fu (1995)
 */
double** getSigma(int n, double theta) {
	return getCentralMom(n, theta, QUADRATIC);
}

/**
 * returns E[xi_i xi_j]-E[xi_i]E[xi_j]
 */
double** getCov(int n, double theta) {
	return getCentralMom(n, theta, ANY);
}

/**
 * returns two-loci spectrum; the same as E2,
 * except that values on diagonal line are counted only half
 */
double** get2loci(int n, double theta, int relation) {
	double **tl = new double*[n - 1];

	for (int i = 0; i < n - 1; i++) {
		tl[i] = new double[n - 1];
	}

	initA(n);

	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			tl[i][j] = tau(n, i + 1, j + 1, relation) * theta * theta;
		}
		tl[i][i] /= 2;
	}
	return tl;
}

/**
 * calculates the above values for a folded spectrum
 */
void fold(int n, double** &unfolded) {
	int size = n / 2;
	double **folded = new double*[size];
	for (int i = 0; i < size; i++) {
		folded[i] = new double[size];
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			folded[i][j] = unfolded[i][j] + unfolded[i][n - 2 - j] + unfolded[n - 2 - i][j]
					+ unfolded[n - 2 - i][n - 2 - j];
			folded[i][j] /= (1 + (i == n - 2 - i ? 1 : 0)) * (1 + (j == n - 2 - j ? 1 : 0));
		}
	}

	for (int i = 0; i < n - 1; i++) {
		delete[] unfolded[i];
	}
	delete[] unfolded;

	unfolded = folded;
}
