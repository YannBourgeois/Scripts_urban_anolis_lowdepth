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
 *  Created on: 18.05.2012
 *      Author: klasal
 *
 *  Filters output of ms; retains only segregating sites between
 *  specified positions
 */

#include <iostream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>

using namespace std;

void scan();

double lb=0, rb=1;

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cout << "Usage: mscut [lower boundary] [upper boundary]\n";
		return (1);
	}
	lb = atof(argv[1]);
	rb = atof(argv[2]);

	std::ios_base::sync_with_stdio(true);
	cout << fixed << setprecision(4);

	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header){
	cout << header << "Window: [" << lb << ',' << rb << "[\n\n";
}

void processTree(const char*){
}

void processSample(const int sampleSize, const int segsites, const double* positions,
		bool** polymorphisms){
	cout << "\n//\n";
	int lbIndex = 0;
	int rbIndex = 0;

	for (int i = 0; i < segsites; ++i) {
		if (positions[i] < lb) {
			lbIndex++;
		}
		if (positions[i] < rb) {
			rbIndex++;
		}
	}

	cout << "segsites: " << rbIndex - lbIndex << '\n';
	if (rbIndex - lbIndex > 0) {
		cout << "positions:";
		for (int i = lbIndex; i < rbIndex; ++i) {
			cout << " " << positions[i];
		}

		cout << '\n';

		for (int i = 0; i < sampleSize; ++i) {
			for (int j = lbIndex; j < rbIndex; j++) {
				cout << polymorphisms[i][j];
			}
			cout << '\n';
		}
	} else {
		cout << '\n';
	}
}
