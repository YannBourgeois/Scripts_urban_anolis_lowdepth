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
 *  Copyright 2016 Alexander Klassmann
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

using namespace std;

#define max_size 100000

int read2dArray(char *fileName, int &size, double** &m) {
	ifstream file;
	if (fileName != NULL) {
		file.open(fileName, ios::in);
	}
	istream& in = fileName != NULL ? file : cin;

	double* firstRow;
	string line;

	if (in) {
		getline(in, line);

		if (in.eof()) {
			cout << "Error: array misspecified" << endl;
			exit(1);
		}

		while (std::isspace(*line.rbegin()))
			line.erase(line.length() - 1);

		stringstream is(line);
		size = 0;
		firstRow = new double[max_size];
		while (!is.eof()) {
			if (size == max_size) {
				cout << "Error: array bigger than " << max_size << endl;
				exit(1);
			}
			is >> firstRow[size++];
		}

		m = new double*[size];
		for (int i = 0; i < size; i++) {
			m[i] = new double[size];
		}

		for (int j = 0; j < size; j++) {
			m[0][j] = firstRow[j];
		}

		delete[] firstRow;

		for (int i = 1; i < size; i++) {
			getline(in, line);
			if (in.eof()) {
				cout << "Error: array misspecified" << endl;
				exit(1);
			}
			int j = 0;
			while (std::isspace(*line.rbegin()))
				line.erase(line.length() - 1);

			if (line.length() > 0) {
				stringstream is(line);
				while (!is.eof()) {
					is >> m[i][j++];
				}
			} else {
				i--;
			}
		}
	} else {
		cout << "Error: could not read file";
		if (fileName != NULL) {
			cout << " " << fileName;
		}
		cout << endl;
		exit(1);
	}
	return (0);
}

void write2dArray(int size, double** &m) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << m[i][j];
			if (j < size - 1) {
				cout << " ";
			}
		}
		cout << endl;
	}
}
