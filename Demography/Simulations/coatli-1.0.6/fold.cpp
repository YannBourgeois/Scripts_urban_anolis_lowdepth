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
 * Folds an unfolded spectrum
 *
 *  Created on: 10.09.16
 *      Author: AlexanderKlassmann
 */
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

string sep = "|";

void printHelp() {
	cerr << "Usage: cat unfolded.fs | fold" << endl;
}

int main(int argc, char *argv[]) {
	vector<double> nestedFs;
	vector<double> disjointFs;
	string tmp;
	string line;

	while (getline(cin, line)) {
		nestedFs.clear();
		disjointFs.clear();
		stringstream is(line);
		double m;
		while (is >> m) {
			nestedFs.push_back(m);
		}
		if (!is.eof()) {
			is.clear();
			is >> tmp;
			if (tmp != sep) {
				cout << "Error: invalid separator '" << tmp << "'." << endl;
				return (1);
			}
			while (is >> m) {
				disjointFs.push_back(m);
			}
		}

		if (nestedFs.size() > 0) {
			for (int i = 0; i < nestedFs.size() / 2; i++) {
				cout << nestedFs[i] + nestedFs[nestedFs.size() - 1 - i];
				if (i < (nestedFs.size() - 1) / 2) {
					cout << " ";
				}
			}
			if (nestedFs.size() % 2 == 1) {
				cout << nestedFs[nestedFs.size() / 2];
			}
			if (disjointFs.size() == 0) {
				cout << endl;
			} else {
				cout << " " << sep << " ";
				for (int i = 0; i < disjointFs.size() / 2; i++) {
					cout << disjointFs[i] + disjointFs[disjointFs.size() - 1 - i];
					if (i < (disjointFs.size() - 1) / 2) {
						cout << " ";
					}
				}
				if (disjointFs.size() % 2 == 1) {
					cout << disjointFs[disjointFs.size() / 2];
				}
				cout << endl;
			}
		}
	}
	return (0);
}
