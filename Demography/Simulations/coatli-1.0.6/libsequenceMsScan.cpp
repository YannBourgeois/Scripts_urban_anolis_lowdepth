/*
 *  This file is part of ntx.
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
 *  along with ntx.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Parses output of 'ms' using C library "libsequence".
 * It relies on the library "libsequence" (Thornton, K. (2003) "libsequence:
 * a C++ class library for evolutionary genetic analysis." Bioinformatics 19(17): 2325-2327),
 * downloadable from https://molpopgen.github.io/libsequence/
 *
 * Libsequence allows to handle different input than 'ms', however the ms command line parameters seem
 * not to be parsed and using option -T of 'ms' (output of the coalescent tree) leads to an error;
 * furthermore, it seems that only the first 'run' is parsed.
 *
 * In order to use this parser out-comment the lines in the Makefile header accordingly
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <Sequence/VariantMatrix.hpp>
#include <Sequence/VariantMatrixViews.hpp>
#include <Sequence/variant_matrix/msformat.hpp>

using namespace std;
using namespace Sequence;

void processHeader(const int sampleSize, const int samples, const char *header);
void processSample(const int sampleSize, const int segsites, const double *positions, bool **polymorphisms);

void scan() {

	std::ios_base::sync_with_stdio(true);

	VariantMatrix vm = Sequence::from_msformat(std::cin);

	/* libsequence seems not to parse command line parameters and only data of the first 'run' */
	char header[] = "ms dummy 1";

	processHeader(vm.nsam(), 1, header);

	if (vm.nsam() > 0) {
		/* translate numbers into booleans */
		bool **polymorphisms = new bool*[vm.nsam()];

		for (int i = 0; i < vm.nsam(); i++) {
			polymorphisms[i] = new bool[vm.nsites()];
		}

		for (int i = 0; i < vm.nsam(); i++) {
			for (int j = 0; j < vm.nsites(); j++) {
				polymorphisms[i][j] = (vm.cat(j, i) == 1);
			}
		}

		processSample(vm.nsam(), vm.nsites(), vm.cpbegin(), polymorphisms);

		for (int i = 0; i < vm.nsam(); i++) {
			delete[] polymorphisms[i];
		}
		delete[] polymorphisms;
	}
}
