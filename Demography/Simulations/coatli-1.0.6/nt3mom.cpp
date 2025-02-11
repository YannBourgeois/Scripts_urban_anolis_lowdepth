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
 *  Prints the third moments of the frequency spectrum with
 *  constant population.
 *
 *  Copyright 2016 Alexander Klassmann
 */
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "ntdef.h"

double** getCondCov(int, int, double, int, int);
double** getCondE2(int, int, double, int, int);
double*** getE3(int, double, int, bool);
double getE3element(int, double, int, bool, int, int, int);
void fold(int, double*** &);
void write3dArray(int, double***&);

using namespace std;

void printUsage() {
	cerr << "Usage: nt3mom n [k] [-p] [-t] [-mu] ..." << endl;
	cerr << "n sample size; print third moments E[xi_i xi_j xi_k]" << endl;
	cerr << "k print expected values conditional on mutations of size k (E[xi_i xi_j|xi_k])" << endl;
	cerr << "-p precision (number of decimals, default " << DEFAULT_PRECISION << ")" << endl;
	cerr << "-t theta (default 1)" << endl;
	cerr << "-ijk i j k calculate only the single element E[xi_i xi_j xi_k]" << endl;
	cerr << "-cov print covariance cov[xi_i,xi_j] conditional on mutations of size k" << endl;
	cerr << "-nn or -aa print conditional moments for two nested mutations" << endl;
	cerr << "-nd or -ab print conditional moments with first nested, second disjoint" << endl;
	cerr << "-dn or -ba print conditional moments with first disjoint, second nested" << endl;
	cerr << "-dd or -bb print conditional moments for two disjoint mutations" << endl;
	cerr << "-mu print central moments mu3[xi_i xi_j xi_k]" << endl;
	cerr << "-f print folded (nested, disjoint, etc.) conditional moments" << endl;
	cerr << "-1 use only terms linear in theta" << endl;
	cerr << "-2 use only terms quadratic in theta" << endl;
	cerr << "-3 use only terms cubic in theta" << endl;
}

int main(int argc, char *argv[]) {
	double theta = 1.0;
	bool cov = false;
	bool folded = false;
	bool central = false;
	int thetaExp = ANY;
	int k = 0, n = 0, p = 0, fs_size, relation = 0, ii = 0, jj = 0, kk = 0;

	if (argc == 1) {
		printUsage();
		exit(1);
	}

	for (int i = 1; i < argc; i++) {
		string s(argv[i]);
		if (s == "-p" && i + 1 < argc) {
			i++;
			p = atoi(argv[i]);
		} else if (s == "-t" && i + 1 < argc) {
			i++;
			theta = atof(argv[i]);
		} else if (s == "-cov") {
			cov = true;
		} else if (s == "-f") {
			folded = true;
		} else if (s == "-nn" || s == "-aa") {
			relation = NESTED_NESTED;
		} else if (s == "-dd" || s == "-bb") {
			relation = DISJOINT_DISJOINT;
		} else if (s == "-nd" || s == "-ab") {
			relation = NESTED_DISJOINT;
		} else if (s == "-dn" || s == "-ba") {
			relation = DISJOINT_NESTED;
		} else if (s == "-1") {
			thetaExp = LINEAR;
		} else if (s == "-2") {
			thetaExp = QUADRATIC;
		} else if (s == "-3") {
			thetaExp = CUBIC;
		} else if (((s == "-ijk") || (s == "-hij")) && i + 3 < argc) {
			i++;
			ii = atoi(argv[i]);
			if (ii < 1 || ii >= n) {
				printUsage();
				exit(1);
			}
			i++;
			jj = atoi(argv[i]);
			if (jj < 1 || jj >= n) {
				printUsage();
				exit(1);
			}
			i++;
			kk = atoi(argv[i]);
			if (kk < 1 || kk >= n) {
				printUsage();
				exit(1);
			}
		} else if (s == "-mu" || s == "-mu3") {
			central = true;
		} else if (n == 0) {
			n = atoi(argv[i]);
			if (n <= 1) {
				printUsage();
				exit(1);
			}
		} else if (k == 0) {
			k = atoi(argv[i]);
			if (k < 1 || k >= n) {
				printUsage();
				exit(1);
			}
		} else {
			printUsage();
			exit(1);
		}
	}

	if (theta < 0 || (relation != ANY && k == 0) || (k == 0 && cov) || (k > 0 && ii > 0)) {
		printUsage();
		exit(1);
	}

	if (p <= 0)
		p = DEFAULT_PRECISION;

	fs_size = n - 1;

	cout << fixed << setprecision(p);

	if (ii > 0) {
		cout << getE3element(n, theta, thetaExp, central, ii, jj, kk) << endl;
	} else if (k > 0) {
		double** condCov =
				cov ? getCondCov(n, k, theta, relation, thetaExp) : getCondE2(n, k, theta, relation, thetaExp);
		if (folded) {
			int nrow, ncol;
			if (relation == 1 || relation == 3) {
				nrow = k;
			} else {
				nrow = n - k;
			}
			if (relation == 1 || relation == 4) {
				ncol = k;
			} else {
				ncol = n - k;
			}

			double** foldedCondCov = new double*[nrow / 2];
			for (int i = 0; i < nrow / 2; i++) {
				foldedCondCov[i] = new double[ncol / 2];
			}

			for (int i = 0; i < nrow / 2; i++) {
				for (int j = 0; j < ncol / 2; j++) {
					foldedCondCov[i][j] = condCov[i][j];
					if (i != nrow - 2 - i) {
						foldedCondCov[i][j] += condCov[nrow - 2 - i][j];
					}
					if (j != ncol - 2 - j) {
						foldedCondCov[i][j] += condCov[i][ncol - 2 - j];
						if (i != nrow - 2 - i) {
							foldedCondCov[i][j] += condCov[nrow - 2 - i][ncol - 2 - j];
						}
					}
				}
			}

			for (int i = 0; i < nrow / 2; i++) {
				for (int j = 0; j < ncol / 2; j++) {
					cout << foldedCondCov[i][j];
					if (j < ncol / 2 - 1)
						cout << "\t";
				}
				cout << endl;
			}

			for (int i = 0; i < nrow / 2; i++) {
				delete[] foldedCondCov[i];
			}
			delete[] foldedCondCov;

		} else {
			for (int i = 0; i < fs_size; i++) {
				for (int j = 0; j < fs_size; j++) {
					cout << condCov[i][j];
					if (j < fs_size - 1)
						cout << "\t";
				}
				cout << endl;
			}
		}
		for (int i = 0; i < fs_size; i++) {
			delete[] condCov[i];
		}
		delete[] condCov;
	} else {
		double*** e3 = getE3(n, theta, thetaExp, central);

		if (folded) {
			fold(n, e3);
			fs_size /= 2;
		}

		write3dArray(fs_size, e3);

		for (int i = 0; i < fs_size; i++) {
			for (int j = 0; j < fs_size; j++) {
				delete[] e3[i][j];
			}
			delete[] e3[i];
		}
		delete[] e3;
	}
	return (0);
}
