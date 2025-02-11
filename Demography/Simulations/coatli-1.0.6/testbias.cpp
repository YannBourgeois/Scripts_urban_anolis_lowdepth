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
 *
 *  Prints the moments of the number of segregating sites S
 *  as calculated by Watterson 1975 and Klassmann & Ferretti 2016
 *
 *  Copyright 2016 Alexander Klassmann
 */
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>

double*** getE3(int, double);
double** getE2(int, double);
double** getSigma(int, double);
double* getE1(int, double);

/*
 * output precision
 */
#define DEFAULT_PRECISION 4

using namespace std;

void printUsage() {
	cerr << "Usage: skewness n [-p] [-t]" << endl;
	cerr << "n sample size; with k=0, print third moments" << endl;
	cerr << "-p precision (number of decimals, default " << DEFAULT_PRECISION << ")" << endl;
	cerr << "-t theta (default 1)" << endl;
}

/*
 * compute Omega for unfolded (full)  tests
 */
double* getOmega(int n, string test) {
	double *w1, *w2, *omega, sum_w1 = 0, sum_w2 = 0;
	int i;

	w1 = new double[n - 1];
	w2 = new double[n - 1];

	omega = new double[n - 1];

	if (test == "D") {
		for (i = 1; i < n; i++) {
			w1[i - 1] = n - i;
			w2[i - 1] = 1.0 / i;
		}
	} else if (test == "E") {
		for (i = 1; i < n; i++) {
			w1[i - 1] = 1.0 / (n - 1);
			w2[i - 1] = 1.0 / i;
		}
	} else if (test == "D2") {
		for (i = 1; i < n; i++) {
			w1[i - 1] = 1.0 / i;
			w2[i - 1] = (i == 1) ? 1.0 : 0.0;
		}
	} else if (test == "F") {
		for (i = 1; i < n; i++) {
			w1[i - 1] = n - i;
			w2[i - 1] = (i == 1) ? 1.0 : 0.0;
		}
	} else if (test == "H") {
		for (i = 1; i < n; i++) {
			w1[i - 1] = n - i;
			w2[i - 1] = i;
		}
	}

	for (i = 0; i < n - 1; i++) {
		sum_w1 += w1[i];
		sum_w2 += w2[i];
	}

	for (i = 0; i < n - 1; i++) {
		omega[i] = (w1[i] / sum_w1 - w2[i] / sum_w2);
	}

	delete[] w1;
	delete[] w2;

	return omega;
}

double getSkewness(int n, double theta, double* omega, double*** e3) {
	double** sigma = getSigma(n, 1.);
	double* xi0 = getE1(n, 1.);
	double s = 0, s1 = 0, s2 = 0;
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			for (int k = 0; k < n - 1; k++) {
				s += omega[i] * omega[j] * omega[k] * e3[i][j][k] / (xi0[i] * xi0[j] * xi0[k]);
			}
		}
	}
	for (int i = 0; i < n - 1; i++) {
		s1 += omega[i] * omega[i] / xi0[i];
		for (int j = 0; j < n - 1; j++) {
			s2 += omega[i] / xi0[i] * omega[j] / xi0[j] * sigma[i][j];
		}
	}
	double var = sqrt(theta * s1 + theta * theta * s2);

	delete[] omega;
	for (int i = 0; i < n - 1; i++) {
		delete[] sigma[i];
	}
	delete[] sigma;

	return s / (var * var * var);
}

double getBias(int n, double theta, double* omega, double** e2, double*** e3) {
	double** sigma = getSigma(n, 1.);
	double* xi0 = getE1(n, 1.);
	double* xiS = new double[n - 1];
	double* xiS2 = new double[n - 1];
	double s = 0, s1 = 0, s2 = 0;

	double an = 0, bn = 0;
	for (int i = 1; i <= n - 1; i++) {
		an += 1. / i;
		bn += 1. / (i * i);
	}

	for (int k = 0; k < n - 1; k++) {
		xiS[k] = 0;
		xiS2[k] = 0;
		for (int i = 0; i < n - 1; i++) {
			xiS[k] += e2[k][i];
			for (int j = 0; j < n - 1; j++) {
				xiS2[k] += e3[k][i][j];
			}
		}
	}

	for (int i = 0; i < n - 1; i++) {
		s1 += omega[i] * omega[i] / xi0[i];
		for (int j = 0; j < n - 1; j++) {
			s2 += omega[i] / xi0[i] * omega[j] / xi0[j] * sigma[i][j];
		}
	}

	for (int k = 0; k < n - 1; k++) {
		s += omega[k] / xi0[k] * (xiS[k] / an * s1 + (xiS2[k] - xiS[k]) / (an * an + bn) * s2);
	}

	double var = sqrt(theta * s1 + theta * theta * s2);

	delete[] omega;
	for (int i = 0; i < n - 1; i++) {
		delete[] sigma[i];
	}
	delete[] sigma;
	delete[] xi0;
	delete[] xiS2;
	delete[] xiS;

	return -s / (2 * var * var * var);
}

int main(int argc, char *argv[]) {
	double theta = 1.0;
	int n = 0, p = 0;

	if (argc == 1) {
		printUsage();
		exit(1);
	}

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else if (string(argv[i]) == "-t" && i + 1 < argc) {
			i++;
			theta = atof(argv[i]);
		} else if (n == 0) {
			n = atoi(argv[i]);
			if (n <= 1) {
				printUsage();
				exit(1);
			}
		} else {
			printUsage();
			exit(1);
		}
	}

	if (theta < 0) {
		printUsage();
		exit(1);
	}

	if (p <= 0)
		p = DEFAULT_PRECISION;

	cout << fixed << setprecision(p);

	double*** e3 = getE3(n, theta);
	double** e2 = getE2(n, theta);

	cout << "Test\tSkewness\tBias" << endl;
	cout << "Tajima's D:\t" << setw(p + 3) << getSkewness(n, theta, getOmega(n, "D"), e3) << "\t" << setw(p + 3)
			<< getBias(n, theta, getOmega(n, "D"), e2, e3) << endl;
	cout << "Fu&Li's D:\t" << setw(p + 3) << getSkewness(n, theta, getOmega(n, "D2"), e3) << "\t" << setw(p + 3)
			<< getBias(n, theta, getOmega(n, "D2"), e2, e3) << endl;
	cout << "Fu&Li's F:\t" << setw(p + 3) << getSkewness(n, theta, getOmega(n, "F"), e3) << "\t" << setw(p + 3)
			<< getBias(n, theta, getOmega(n, "F"), e2, e3) << endl;
	cout << "Fay&Wu's H:\t" << setw(p + 3) << getSkewness(n, theta, getOmega(n, "H"), e3) << "\t" << setw(p + 3)
			<< getBias(n, theta, getOmega(n, "H"), e2, e3) << endl;
	cout << "Zeng's E:\t" << setw(p + 3) << getSkewness(n, theta, getOmega(n, "E"), e3) << "\t" << setw(p + 3)
			<< getBias(n, theta, getOmega(n, "E"), e2, e3) << endl;

	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < n - 1; j++) {
			delete[] e3[i][j];
		}
		delete[] e2[i];
	}
	delete[] e3;
	delete[] e2;

	return (0);
}
