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
 * The program computes the frequency spectra using formulas of (1).
 * The spectra are derived there under the condition of exactly two
 * nested mutations within a coalescent and vanishing theta.
 *
 * The program prints the spectra
 * - d2: the number of individuals with both mutations
 * - d1: the number of individuals containing only one (the older) mutation
 * - the standard frequency spectrum.
 *
 * The spectrum d0 of indivuals with none of the two mutations is equal to d1.
 *
 * Reference:
 * (1) A. Hobolth & C. Wiuf:
 * The genealogy, site frequency spectrum and ages of two nested mutant alleles.
 * Theoretical Population Biology 2009.
 *
 *
 *  Created on: 08.10.2013
 *      Author: Alexander Klassmann
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

double binomial(int n, int k) {
	double b = 1.0;

	// binomial(n, k) = binomial(n, n-k)
	if (k > n - k)
		k = n - k;

	// Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
	for (int i = 0; i < k; ++i) {
		b *= (n - i);
		b /= (i + 1);
	}

	return b;
}

int main(int argc, char *argv[]) {
	unsigned int n, d0, d1, d2 = 0, l;
	double c, h;
	double *p1, *p2;

	if (argc != 2) {
		cout << "Usage: M2fs [sample size]";
		exit(1);
	}

	n = atoi(argv[1]);

	p1 = new double[n - 2];
	p2 = new double[n - 2];

	h = 0;
	for (l = 1; l < n; l++) {
		h += 1. / l;
	}
	c = 1. / ((n - 1) * h - 2. * (n - 1.) * (n - 1.) / n);

	cout << fixed << setprecision(4);

	cout << "d2: ";

	// eq (23) of (1)
	for (d2 = 1; d2 < n - 1; d2++) {
		double p = 0;
		for (l = 3; l <= n - d2 + 1; l++) {
			p += ((l - 2.) / l) * binomial(n - d2 - 1, l - 2)
					/ binomial(n - 2, l - 2);
		}
		p2[d2 - 1] = p * c;
		cout << p2[d2 - 1];
		if (d2 < n - 2)
			cout << " ";
		else
			cout << endl;
	}

	cout << "d1: ";

	// corrected version of eq (24) of (1)
	for (d1 = 1; d1 < n - 1; d1++) {
		double p = 0;
		for (d0 = 1; d0 < n - d1; d0++) {
			for (l = 3; l <= d0 + d1 + 1; l++) {
				p +=
						(1. / l)
								* (binomial(d0 + d1 - 2, l - 3)
										/ binomial(n - 2, l - 2));
			}
		}
		p1[d1 - 1] = p * c;
		cout << p1[d1 - 1];
		if (d1 < n - 2)
			cout << " ";
		else
			cout << endl;
	}

	cout << "fs: ";
	for (int i = 0; i < n - 1; i++) {
		double p = p2[i];
		for (int j = 0; j < i; j++) {
			p += p2[j] * 1. / (n - j - 2);
		}
		cout << p;
		if (i < n - 2)
			cout << " ";
		else
			cout << endl;
	}
}
