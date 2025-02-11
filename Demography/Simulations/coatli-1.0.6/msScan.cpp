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
 *  Copyright 2014 Alexander Klassmann
 */

#include <fstream>
#include <cstdlib>
#include <cstring>

void processHeader(const int sampleSize, const int samples, const char *header);

void processTree(const char*);

void processSample(const int sampleSize, const int segsites, const double *positions, bool** polymorphisms);

void scan() {
	int maxsegsites = 1000;
	int i, j, sampleSize, samples, segsites;
	bool **polymorphisms;
	double *positions;
	char dummy[20], line[2048], *polyline;
	FILE *file;

	file = stdin;
	fgets(line, 1000, file);

	sscanf(line, "%s %d %d", dummy, &sampleSize, &samples);

	line[strlen(line) - 1] = '\0';
	processHeader(sampleSize, samples, line);

	polyline = new char[maxsegsites + 2];
	positions = new double[maxsegsites];
	polymorphisms = new bool*[sampleSize];
	for (i = 0; i < sampleSize; i++) {
		polymorphisms[i] = new bool[maxsegsites];
	}

	while (fgets(line, 1000, file) != NULL) {
		if (line[0] == '(') {
			processTree(line);
		} else if (line[0] == 's') {

			sscanf(line, "segsites: %d", &segsites);
			if (segsites > maxsegsites) {
				delete[] polyline;
				delete[] positions;
				maxsegsites = 2 * segsites;
				polyline = new char[maxsegsites + 2];
				positions = new double[maxsegsites];
				for (i = 0; i < sampleSize; i++) {
					delete[] polymorphisms[i];
					polymorphisms[i] = new bool[maxsegsites];
				}
			}

			if (segsites > 0) {
				fscanf(file, "%s", dummy);

				for (i = 0; i < segsites; i++) {
					fscanf(file, "%lf", &positions[i]);
				}
				for (i = 0; i < sampleSize; i++) {
					fscanf(file, "%s", polyline);
					for (j = 0; j < segsites; j++) {
						polymorphisms[i][j] = (polyline[j] == '1');
					}
				}
			}

			processSample(sampleSize, segsites, positions, polymorphisms);

		}
	}
	for (i = 0; i < sampleSize; i++) {
		delete[] polymorphisms[i];
	}
	delete[] polyline;
	delete[] polymorphisms;
	delete[] positions;
}
