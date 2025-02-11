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
 *  Draws 'ms'-output as aligned mutations in latex-tikz
 *
 *  Created on: 21.05.14
 *      Author: Alexander Klassmann
 */
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "tikzdef.h"

using namespace std;

void scan();

const double default_scale = 1.0;
const double relative_radius = 0.45;
const char* colors[] = { "blue", "brown", "cyan", "green", "lime", "magenta", "olive", "orange", "pink", "purple",
		"red", "teal", "violet", "yellow" };
const int ncolors = 14;
bool usePositions = false;
bool useDifferentColors = false;
bool useMutationSizeColors = false;
bool printXaxis = false;
double scale = default_scale;
double rscale = default_scale;

void processHeader(const int sampleSize, const int samples, const char* header) {
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** polymorphisms) {
	double radius;
	if (segsites > sampleSize) {
		radius = relative_radius / segsites;
	} else {
		radius = relative_radius / sampleSize;
	}
	int* colorScheme = new int[segsites];
	if (useDifferentColors) {
		for (int i = 0; i < segsites; i++) {
			colorScheme[i] = i + 1;
		}
	} else if (useMutationSizeColors) {
		for (int j = 0; j < segsites; j++) {
			int s = 0;
			for (int i = 0; i < sampleSize; i++) {
				if (polymorphisms[i][j]) {
					s++;
				}
			}
			colorScheme[j] = s;
		}
	}

	cout << "\\begin{tikzpicture}[color=" << linecolor;
	if (!useDifferentColors) {
		cout << ", ball color=" << ballcolor;
	}
	cout << "]" << endl;
	cout << "\\def\\radius{" << radius * 10 * scale << "}" << endl;
	cout << "\\def\\radiusscale{" << rscale << "}" << endl;

	if (useDifferentColors || useMutationSizeColors) {
		int max = useDifferentColors ? segsites : sampleSize;
		for (int i = 0; i < max; i++) {
			cout << "\\colorlet{m" << (i + 1) << "}{" << colors[i % ncolors] << "}" << endl;
		}
	}

	cout << "\\foreach \\i in {1,...," << sampleSize << "}{" << endl;
	cout << "\\draw (0,{((\\i-" << 0.25 << ")*\\radius*2.2)})--(" << 10 * scale << ",{((\\i-" << 0.25
			<< ")*\\radius*2.2)});" << endl;
	cout << "}" << endl;

	if (printXaxis) {
		cout << "\\draw" << endl;
		cout << " (0,0)--(" << 10 * scale << ",0)" << endl;
		int range = usePositions ? 10 : segsites;
		for (int j = 0; j < range; j++) {
			double x;
			if (usePositions) {
				x = double(j) / range * 10 * scale;
			} else {
				x = ((j + 0.5) / range) * 10 * scale;
			}
			cout << " (" << x << ", " << 0.15 * scale << ")--(" << x << ", " << -0.15 * scale << ") node[below]{";
			if (usePositions) {
				if (j == 0) {
					cout << "0";
				} else if (j < 10) {
					cout << "0." << j;
				}
			} else {
				cout << j + 1;
			}
			cout << "}" << endl;
		}
		if (usePositions) {
			cout << " (10, " << 0.15 * scale << ")--(10, " << -0.15 * scale << ") node[below]{1}";
		}
		cout << ";" << endl;
	}

	for (int i = 0; i < sampleSize; i++) {
		for (int j = 0; j < segsites; j++) {
			if (polymorphisms[i][j]) {
				double x;

				if (usePositions) {
					x = positions[j] * 10 * scale;
				} else {
					x = ((j + 0.5) * 1.0 / segsites) * 10 * scale;
				}
				cout << "\\shade";
				if (useDifferentColors || useMutationSizeColors) {
					cout << "[ball color=m" << colorScheme[j] << "] ";
				}
				cout << "(" << x << ",{(" << sampleSize - i - 0.25
						<< "*\\radius*2.2)}) circle({\\radius*\\radiusscale});" << endl;
			}
		}
	}
	cout << "\\end{tikzpicture}" << endl;
	delete[] colorScheme;
}

void printUsage() {
	cout << "Usage: ms2tikzhaplo -c -colors -pos -r [scale radius] -s [scale]" << endl;
	cout << "Options:" << endl;
	cout << " -c complete LaTeX document (for standalone compiling)" << endl;
	cout << " -colors use different colors for each mutation" << endl;
	cout << " -colors2 color mutations with respect to their frequency" << endl;
	cout << " -r scale radius of mutation balls (default 1)" << endl;
	cout << " -s scale size of whole picture (default 1)" << endl;
	cout << " -xaxis draw x-axis below sequences" << endl;
	cout << " -pos use mutation positions from 'ms' (instead of equally spaced)" << endl;
	cout << " -hap input consists of haplotype data only (without 'ms' header)" << endl;
}

int main(int argc, char* argv[]) {
	bool hap = false, complete = false;
	for (int i = 1; i < argc; i++) {
		string s = string(argv[i]);
		if (s == "-hap") {
			hap = true;
		} else if (s == "-pos") {
			usePositions = true;
		} else if (s == "-c") {
			complete = true;
		} else if (s == "-xaxis") {
			printXaxis = true;
		} else if (s == "-colors") {
			useDifferentColors = true;
		} else if (s == "-colors2") {
			useMutationSizeColors = true;
		} else if (s == "-r" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> rscale;
		} else if (s == "-s" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> scale;
		} else {
			printUsage();
			return (0);
		}
	}
	if ((hap && usePositions) || rscale <= 0 || scale <= 0) {
		printUsage();
	}

	std::ios_base::sync_with_stdio(true);

	if (complete) {
		cout << "\\documentclass{standalone}" << endl;
		cout << "\\usepackage{tikz}" << endl;
		cout << "\\begin{document}" << endl;
	}
	cout << fixed << setprecision(precision);

	if (hap) {
		int a;
		string line;
		vector < vector<bool> > haplotypes;
		while (getline(cin, line)) {
			stringstream is(line);
			vector<bool> v;
			while (!is.eof()) {
				is >> a;
				v.push_back(a != 0);
			}
			haplotypes.push_back(v);
		}
		if (haplotypes.size() > 0 && haplotypes[0].size() > 0) {
			double* positions = NULL;
			bool** polymorphisms = new bool*[haplotypes.size()];
			for (int i = 0; i < haplotypes.size(); i++) {
				polymorphisms[i] = new bool[haplotypes[i].size()];
				for (int j = 0; j < haplotypes[i].size(); j++) {
					polymorphisms[i][j] = haplotypes[i][j];
				}
			}

			processSample(haplotypes.size(), haplotypes[0].size(), positions, polymorphisms);

			for (int i = 0; i < haplotypes.size(); i++) {
				delete[] polymorphisms[i];
			}
			delete[] polymorphisms;
			delete[] positions;
		}
	} else {
		scan();
	}
	if (complete) {
		cout << "\\end{document}" << endl;
	}
	return 0;
}
