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
 *  This class computes the cubic terms of the third moments of the unfolded frequency spectrum
 *  under the assumption of constant population size.
 *
 *  Copyright 2016 Alexander Klassmann
 */

#include "ntdef.h"
#include <iostream>

namespace {
double **alpha, **beta, **alpha2, **beta2;

void initAlphaBeta(int n) {
	double *a = new double[n + 1];
	a[0] = 0;
	for (int i = 1; i <= n; i++) {
		a[i] = a[i - 1] + 1. / i;
	}

	alpha = new double*[n + 1];
	beta = new double*[n + 1];
	for (int i = 0; i < n + 1; i++) {
		alpha[i] = new double[n + 1];
		beta[i] = new double[n + 1];
	}
	for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			alpha[i][j] = (a[i - 1] - a[j - 1]) / (i - j);
			beta[i][j] = 2. * i * (a[i] - a[j - 1]) / ((i - j + 1.) * (i - j)) - 2.0 / (i - j);
		}
	}

	delete[] a;
}

void initAlphaBeta2(int n) {
	alpha2 = new double*[n];
	beta2 = new double*[n];
	double s;

	//Calculate Binomial coefficients
	double** b = new double*[n];
	for (int i = 0; i < n; i++) {
		b[i] = new double[n];
	}
	b[0][0] = 1;
	for (int i = 1; i < n; i++) {
		b[i][0] = 1;
		b[0][i] = 0;
	}
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			b[i][j] = b[i - 1][j] + b[i - 1][j - 1];
		}
	}
	//end of binomial coefficient calculation

	for (int i = 0; i < n; i++) {
		alpha2[i] = new double[n];
		beta2[i] = new double[n];
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			alpha2[i][j] = beta2[i][j] = 0;
		}
	}

	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n - i + 1; j++) {
			for (int k = 2; k <= n; k++) {
				for (int t = 1; t <= k - 1; t++) {
					s = b[i - 1][t - 1] * b[n - i - j][k - t - 1] / b[n - 1][k - 1] / (k * (k - 1));
					alpha2[i][j] += s * alpha[k][t];
					beta2[i][j] += s * beta[k][t] / 2;
				}
			}
		}
	}

	//delete binomial coefficients
	for (int i = 0; i < n; i++) {
		delete[] b[i];
	}
	delete[] b;
}

double alpha3(int n, int h, int i, int j) {
	return (h + 1) * alpha2[i][j] - 2 * h * alpha2[i][j + 1] + (h - 1) * alpha2[i][j + 2];
}

double alpha4(int n, int h, int i, int j) {
	return (h + 1) * alpha2[i + 1][j] - 2 * h * alpha2[i][j + 1] + (h - 1) * alpha2[i - 1][j + 2];
}

double beta3(int n, int h, int i, int j) {
	return (h + 1) * beta2[i][j] - 2 * h * beta2[i][j + 1] + (h - 1) * beta2[i][j + 2];
}

double beta4(int n, int h, int i, int j) {
	return (h + 1) * beta2[i + 1][j] - 2 * h * beta2[i][j + 1] + (h - 1) * beta2[i - 1][j + 2];
}

double taa(int n, int h, int i, int j) {
	if (j < i) {
		if (i < h) {
			return beta4(n, i - j, i - j, j) - beta4(n, i - j, i - j + 1, j);
		} else if (h == i) {
			return beta4(n, i - j, i - j, j);
		}
	} else if (j == i) {
		if (i < h) {
			return beta2[1][j] - beta2[2][j];
		} else if (i == h) {
			return beta2[1][j];
		}
	}
	return 0.;
}

double tab(int n, int h, int i, int j) {
	if (i + j < h) {
		return beta3(n, h - i - j, i, j) - beta3(n, h - i - j, i + 1, j);
	} else if (i + j == h) {
		return beta2[i][j] - beta2[i + 1][j];
	}
	return 0.;
}

double tba(int n, int h, int i, int j) {
	if (j < i) {
		if (h + i < n) {
			return (beta[n][j] - beta[n][j + 1] - beta[n][h + j] + beta[n][h + j + 1]) / (2 * h)
					- beta4(n, i - j, i - j, j) + beta4(n, i - j, i - j + 1, j) - beta3(n, i - j, h, j)
					+ beta3(n, i - j, h + 1, j) + beta3(n, n - i - h, j, h) - beta3(n, n - i - h, j + 1, h);
		} else if (h + i == n) {
			return alpha4(n, n - h - j, n - h - j, j) - beta4(n, n - h - j, n - h - j, j) + alpha3(n, n - h - j, h, j)
					- beta3(n, n - h - j, h, j) + beta2[j][h] - beta2[j + 1][h];
		}
	} else if (j == i) {
		if (h + i < n) {
			return (beta[n][j] - beta[n][h + j]) / (2 * h) + beta3(n, n - h - j, j, h) - beta2[h][j] + beta2[h + 1][j]
					- beta2[1][j] + beta2[2][j];
		} else if (h + i == n) {
			return (alpha2[n - j][j] + alpha2[j][n - j]) / 2 + alpha2[1][j] - beta2[1][j];
		}
	}
	return 0.;
}

double tbb(int n, int h, int i, int j) {
	if (h + i + j < n) {
		return 1. / ((h + i + j) * (h + i) * h) + 1. / (i * j * (h + i)) - 1. / (i * h * (i + j))
				- (beta[n][j] - beta[n][j + 1] - beta[n][i + j] + beta[n][i + j + 1]) / (2 * i)
				- beta3(n, n - h - i - j, i, j) + beta3(n, n - h - i - j, i + 1, j);
	} else if (h + i + j == n) {
		return (alpha[n][j] - alpha[n][i + j] - (beta[n][j] - beta[n][i + j]) / 2) / i - beta2[i][j] + beta2[i + 1][j];
	}
	return 0.;
}
}

double tau(int n, int h, int i, int j, int relation) {
	double s = 0;

	if (relation == ANY || relation == NESTED_NESTED) {
		s += taa(n, h, i, j) + taa(n, h, j, i) + taa(n, i, j, h) + taa(n, j, h, i) + taa(n, j, i, h) + taa(n, i, h, j)
				+ tab(n, h, i, j) + tab(n, h, j, i);
	}
	if (relation == ANY || relation == NESTED_DISJOINT) {
		s += tab(n, i, j, h) + tab(n, i, h, j) + tba(n, j, h, i) + tba(n, j, i, h);
	}
	if (relation == ANY || relation == DISJOINT_NESTED) {
		s += tab(n, j, i, h) + tab(n, j, h, i) + tba(n, i, h, j) + tba(n, i, j, h);
	}
	if (relation == ANY || relation == DISJOINT_DISJOINT) {
		s += tba(n, h, i, j) + tba(n, h, j, i) + tbb(n, h, i, j) + tbb(n, i, j, h) + tbb(n, j, h, i) + tbb(n, h, j, i)
				+ tbb(n, j, i, h) + tbb(n, i, h, j);
	}

	return s;
}

void initT3(int n) {
	initAlphaBeta(n);
	initAlphaBeta2(n);
}

void freeT3(int n) {
	for (int i = 0; i < n + 1; i++) {
		delete[] alpha[i];
		delete[] beta[i];
	}
	delete[] alpha;
	delete[] beta;

	for (int i = 0; i < n; i++) {
		delete[] alpha2[i];
		delete[] beta2[i];
	}
	delete[] alpha2;
	delete[] beta2;
}
