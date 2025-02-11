/**
 *  This file is part of ntx.
 *
 *  ntx is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ntx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ntx.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * This class computes first and secondary moments of the (un)folded frequency spectrum
 * under the assumption of constant population size.
 * The formulas are derived in Y. Fu: "Statistical properties of segregating sites",
 * Theoretical population biology, Vol. 48, 1995
 *
 * Copyright 2014 Alexander Klassmann
 */

#define KRONECKER( A, B )  (((A)==(B))?1:0)

/**
 * Harmonic numbers: 1+1/2+1/3+...+1/i
 */
double* getHarmonicNumbers(int n) {
	double* HarmonicNumbers = new double[n];

	if (!HarmonicNumbers)
		return 0;

	HarmonicNumbers[0] = 0;

	for (int i = 1; i < n; i++) {
		HarmonicNumbers[i] = HarmonicNumbers[i - 1] + 1.0 / i;
	}

	return HarmonicNumbers;
}

/**
 * The expected Xi spectrum under constant population size
 * cf. equation (1) of Fu 1995.
 */
double* getXi0(int n) {
	double* xi0 = new double[n - 1];
	for (int i = 1; i < n; i++) {
		xi0[i - 1] = 1. / i;
	}
	return xi0;
}

/**
 * the expected Eta spectrum under constant population size
 * cf. equations (6) and (7) of Fu 1995.
 */
double* getEta0(int n) {
	double* eta0 = new double[n / 2];
	for (int i = 1; i <= n / 2; i++) {
		eta0[i - 1] = (1. / i + 1. / (n - i)) / (1. + KRONECKER(i, n - i));
	}
	return eta0;
}

/*
 * cf. equation (6) of Fu 1995.
 */
double getBeta(int i, double *HarmonicNumbers, int n) {
	double ai = HarmonicNumbers[i - 1], an = HarmonicNumbers[n - 1];
	double beta = 0;

	beta = 2.0 * n * (an + (1.0 / n) - ai) / ((n - i + 1.0) * (n - i))
			- 2.0 / (n - i);

	return beta;
}

/*
 * cf. equation (2) of Fu 1995.
 */
double getSigma_ii(int i, double *HarmonicNumbers, int n) {
	double sigma_ii = 0;
	double ai = HarmonicNumbers[i - 1], an = HarmonicNumbers[n - 1];

	if (2 * i < n) {
		sigma_ii = getBeta(i + 1, HarmonicNumbers, n);
	} else {
		if (2 * i == n) {
			sigma_ii = 2.0 * (an - ai) / (n - i) - 1.0 / (i * i);
		} else {
			sigma_ii = getBeta(i, HarmonicNumbers, n) - 1.0 / (i * i);
		}

	}

	return sigma_ii;
}

/*
 * cf. equation (3) of Fu 1995.
 */
double getSigma_ij(int i, int j, double *HarmonicNumbers, int n) {
	double sigma_ij = 0;

	if (i == j) {
		return getSigma_ii(i, HarmonicNumbers, n);
	}

	if (i < j) {
		int tmp = i;
		i = j;
		j = tmp;
	}

	double ai = HarmonicNumbers[i - 1], aj = HarmonicNumbers[j - 1], an =
			HarmonicNumbers[n - 1];

	if (i + j < n) {
		sigma_ij = (getBeta(i + 1, HarmonicNumbers, n)
				- getBeta(i, HarmonicNumbers, n)) / 2.0;
	} else {
		if (i + j == n) {
			sigma_ij = ((an - ai) / (n - i) + (an - aj) / (n - j))
					- ((getBeta(i, HarmonicNumbers, n)
							+ getBeta(j + 1, HarmonicNumbers, n)) / 2.0)
					- (1.0 / (i * j));
		} else {
			sigma_ij = ((getBeta(j, HarmonicNumbers, n)
					- getBeta(j + 1, HarmonicNumbers, n)) / 2.0)
					- (1.0 / (i * j));
		}
	}
	return sigma_ij;
}

/*
 * cf. equation (9) of Fu 1995
 */
double getRho_ii(int i, double *HarmonicSums, int n) {
	double rho_ii;

	rho_ii = getSigma_ii(i, HarmonicSums, n)
			+ getSigma_ii(n - i, HarmonicSums, n)
			+ 2 * getSigma_ij(i, n - i, HarmonicSums, n);
	rho_ii /= (1.0 + KRONECKER(i, n - i)) * (1.0 + KRONECKER(i, n - i));

	return rho_ii;
}

/*
 * cf. equation (9) of Fu 1995
 */
double getRho_ij(int i, int j, double *HarmonicSums, int n) {
	double rho_ij;

	rho_ij = getSigma_ij(i, j, HarmonicSums, n)
			+ getSigma_ij(i, n - j, HarmonicSums, n)
			+ getSigma_ij(n - i, j, HarmonicSums, n)
			+ getSigma_ij(n - i, n - j, HarmonicSums, n);
	rho_ij /= ((1.0 + KRONECKER(i, n - i)) * (1.0 + KRONECKER(j, n - j)));

	return rho_ij;
}

/**
 * sigma matrix = (sigma)ij
 */
double** getSigma(int n) {
	int i, j;
	double **sigma = new double*[n - 1];
	double *HarmonicNumbers = getHarmonicNumbers(n);

	for (i = 0; i < n - 1; i++) {
		sigma[i] = new double[n - 1];
	}

	for (i = 0; i < n - 1; i++) {
		for (j = i; j < n - 1; j++) {
			if (i == j) {
				sigma[i][i] = getSigma_ii(i + 1, HarmonicNumbers, n);
			} else {
				sigma[j][i] = sigma[i][j] = getSigma_ij(i + 1, j + 1,
						HarmonicNumbers, n);
			}
		}
	}
	delete[] HarmonicNumbers;
	return sigma;
}

double** getCov(int n, double theta) {
	double *xi0 = getXi0(n);
	double **cov = getSigma(n);
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			cov[i][j] *= theta * theta;
		}
		cov[i][i] += xi0[i] * theta;
	}
	delete[] xi0;
	return cov;
}

/**
 * rho matrix = (rho)ij
 */
double** getRho(int n) {
	int i, j;
	double **rho = new double*[n / 2];
	double *HarmonicNumbers = getHarmonicNumbers(n);

	for (i = 0; i < n / 2; i++) {
		rho[i] = new double[n / 2];
	}

	for (i = 0; i < n / 2; i++) {
		for (j = i; j < n / 2; j++) {
			if (i == j) {
				rho[i][i] = getRho_ii(i + 1, HarmonicNumbers, n);
			} else {
				rho[j][i] = rho[i][j] = getRho_ij(i + 1, j + 1, HarmonicNumbers,
						n);
			}
		}
	}
	return rho;
}

/*
 * set weighting for theta estimator in nominator
 */
double* getW(double* fs0, int fs_size) {
	int i;
	double sum = 0;
	double* w = new double[fs_size];

	for (i = 0; i < fs_size; i++) {
		w[i] = fs0[i];
		sum += w[i];
	}
	for (i = 0; i < fs_size; i++) {
		w[i] /= sum;
	}
	return w;
}

/*
 * first coefficient in nominator of test statistic
 */
double getWeightedFirstMoment(double* weight, double* fs0, int fs_size) {
	int i;
	double alpha = 0;

	for (i = 0; i < fs_size; i++) {
		alpha += weight[i] * weight[i] / fs0[i];
	}

	return alpha;
}

/*
 * second coefficient in nominator of test statistic
 */
double getWeightedSecondMoment(double* weight, double* fs0, double** secMom0,
		int fs_size) {
	int i, j;
	double beta = 0;

	for (i = 0; i < fs_size; i++) {
		for (j = 0; j < fs_size; j++) {
			beta += (weight[i] / fs0[i]) * secMom0[i][j] * (weight[j] / fs0[j]);
		}
	}

	return beta;
}
