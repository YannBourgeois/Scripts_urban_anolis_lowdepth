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
 *  This class computes the first moment of the unfolded frequency spectrum
 *  under the assumption of constant population size.
 *
 *  Copyright 2016 Alexander Klassmann
 */


double tau(int i) {
	return 1. / i;
}

double* getE1(int n, double theta) {
	double* xi = new double[n - 1];
	for (int i = 1; i < n; i++) {
		xi[i - 1] = tau(i) * theta;
	}
	return xi;
}

