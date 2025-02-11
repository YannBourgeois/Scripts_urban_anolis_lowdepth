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
 *  Prints the second moments resp. conditional first moment of the frequency spectrum with
 *  constant population size using the formulas given by Fu 1995, Ferretti 2012.
 *
 *  Copyright 2016 Alexander Klassmann
 */
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "ntdef.h"

using namespace std;

void initT2(int);
void freeT2(int);
double* getCondE1(int, int, double, int);
double** getE2(int, double, int, int);
double** get2loci(int, double, int);
double** getCentralMom(int, double, int);
void fold(int, double**&);
void write2dArray(int, double**&);

void printUsage() {
	cout << "Usage: nt2mom n [k] [-p] [-t] ..." << endl;
	cout << "n sample size; print second moments E[xi_i xi_j]" << endl;
	cout << "k print first moments conditional on mutations of size k" << endl;
	cout << "-p precision (number of decimals, default " << DEFAULT_PRECISION << ")" << endl;
	cout << "-t theta (default 1)" << endl;
	cout << "-n or -a print moments for a nested mutations" << endl;
	cout << "-d or -b print moments for a disjoint mutations" << endl;
	cout << "-cov print covariance matrix Cov[xi_i,xi_j]" << endl;
	cout << "-sigma print the coefficients sigma matrix of Fu 95 (quadratic terms of cov)" << endl;
	cout << "-f print folded (shared, extern, etc.) second moments" << endl;
	cout << "-1 use only terms linear in theta" << endl;
	cout << "-2 use only terms quadratic in theta" << endl;
	cout << "-enclosing print contribution of enclosing mutations to the nested spectrum" << endl;
	cout << "-cooccurring print contribution of co-occurring mutations to the nested spectrum" << endl;
	cout << "-strictlynested print contribution of strictly nested mutations to the nested spectrum" << endl;
	cout << "-complementary print contribution of complementary mutations to the disjoint spectrum" << endl;
	cout << "-strictlydisjoint print contribution of strictly disjoint mutations to the disjoint spectrum" << endl;
}

int main(int argc, char *argv[]) {
	double **e2 = NULL, theta = 1.0;
	double *e1_nested = NULL, *e1_disjoint = NULL;
	bool sigma = false;
	bool twoLoci = false;
	bool folded = false;
	bool cov = false;
	int thetaExp = ANY;
	int k = 0, n = 0, p = 0, fs_size, relation = 0;

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
		} else if (s == "-f") {
			folded = true;
		} else if (s == "-n" || s == "-a") {
			relation = NESTED;
		} else if (s == "-d" || s == "-b") {
			relation = DISJOINT;
		} else if (s == "-enclosing") {
			relation = ENCLOSING;
		} else if (s == "-cooccurring") {
			relation = COOCCURRING;
		} else if (s == "-strictlynested") {
			relation = STRICTLY_NESTED;
		} else if (s == "-complementary") {
			relation = COMPLEMENTARY;
		} else if (s == "-strictlydisjoint") {
			relation = STRICTLY_DISJOINT;
		} else if (s == "-1") {
			thetaExp = LINEAR;
		} else if (s == "-2") {
			thetaExp = QUADRATIC;
		} else if (s == "-cov") {
			cov = true;
		} else if (s == "-sigma") {
			sigma = true;
		} else if (s == "-2loci") {
			twoLoci = true;
		} else {
			if (n == 0) {
				n = atoi(argv[i]);
				if (n <= 1) {
					printUsage();
					exit(1);
				}
			} else {
				k = atoi(argv[i]);
				if (k <= 0 || k >= n) {
					printUsage();
					exit(1);
				}
			}
		}
	}

	if (n < 2 || theta < 0 || (folded && (relation != ANY || k > 0)) || (sigma && cov) || (sigma && twoLoci)
			|| (cov && twoLoci) || (k > 0 && (sigma || cov || twoLoci))) {
		printUsage();
		return (0);
	}

	if (p <= 0)
		p = DEFAULT_PRECISION;

	cout << fixed << setprecision(p);

	if (k > 0) {
		if (thetaExp != QUADRATIC) {
			e1_nested = getCondE1(n, k, theta, NESTED);
			e1_disjoint = getCondE1(n, k, theta, DISJOINT);

			if (relation == NESTED || relation == ENCLOSING || relation == COOCCURRING || relation == STRICTLY_NESTED) {
				for (int l = 1; l < n; l++) {
					if ((relation == ENCLOSING && l <= k) || (relation == COOCCURRING && l != k)
							|| (relation == STRICTLY_NESTED && l >= k)) {
						cout << 0;
					} else {
						cout << e1_nested[l - 1];
					}
					if (l < n - 1) {
						cout << " ";
					}
				}
			} else if (relation == DISJOINT || relation == COMPLEMENTARY || relation == STRICTLY_DISJOINT) {
				for (int l = 1; l <= n - k; l++) {
					if ((relation == COMPLEMENTARY && l != n - k) || (relation == STRICTLY_DISJOINT && l == n - k)) {
						cout << 0;
					} else {
						cout << e1_disjoint[l - 1];
					}
					if (l < n - k) {
						cout << " ";
					}
				}
			} else {
				for (int l = 1; l <= n - k; l++) {
					cout << e1_nested[l - 1] + e1_disjoint[l - 1];
					if (l < n - 1) {
						cout << " ";
					}
				}
				for (int l = n - k + 1; l < n; l++) {
					cout << e1_nested[l - 1];
					if (l < n - 1) {
						cout << " ";
					}
				}
			}
			cout << endl;
		}
	} else {
		if (cov) {
			e2 = getCentralMom(n, theta, thetaExp);
		} else if (sigma) {
			e2 = getCentralMom(n, theta, QUADRATIC);
		} else if (twoLoci) {
			e2 = get2loci(n, theta, relation);
		} else {
			e2 = getE2(n, theta, relation, thetaExp);
		}

		if (folded) {
			fold(n, e2);
			fs_size = n / 2;
		} else {
			fs_size = n - 1;
		}
		write2dArray(fs_size, e2);

		for (int i = 0; i < fs_size; i++) {
			delete[] e2[i];
		}

		delete[] e2;
	}

	return (0);
}
