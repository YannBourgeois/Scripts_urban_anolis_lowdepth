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
 *  Created on: 20.10.2013
 *      Author: Alexander Klassmann
 *
 *  Filters out subtrees of a specified size by using the newick-string from ms output with option -T
 *
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

bool getSubtrees(string, unsigned int, vector<string>&);
void getLeavesInTree(string, vector<unsigned int>&);
void scan();

string newick;
unsigned int k = 0, l = 0, n = 0;
bool negation = false, complementTree = false, extractSubtree = false, printAll = false, reprintNewick = false,
		subtreeFound = false, ignoreMultipleSubtrees = false;
vector<string> subtrees;
vector<unsigned int> subtreeLeaves;
vector<unsigned int> specifiedLeaves;

/**
 * reads one row of integers from file
 */
bool readline(char* fileName, vector<unsigned int>& leaves) {
	int leave;
	ifstream file;
	file.open(fileName, ios::in);

	if (file.is_open()) {
		while (file >> leave) {
			leaves.push_back(leave);
		}
	} else {
		cout << "Error: Unable to open file " << fileName;
		return (false);
	}
	return (true);
}

void printUsage() {
	cout << "Usage: msTsub k [-n] [-d] [-a][-c][-i]" << endl;
	cout << "Options: -all print out all (including empty trees, if no subtrees of size k)" << endl;
	cout << "         -n or -a print nested trees of size k" << endl;
	cout << "         -d or -b print disjoint subtrees" << endl;
	cout << "         -T reprint newick string" << endl;
	cout << "         -i ignore multiple subtrees of size k (select one randomly)" << endl;
}

int main(int argc, char *argv[]) {

	for (int i = 1; i < argc; i++) {
		string s(argv[i]);
		if (s == "-all") {
			printAll = true;
		} else if (s == "-d" || s == "-b") {
			complementTree = true;
		} else if (s == "-f") {
			i++;
			if (readline(argv[i], specifiedLeaves)) {
				k = specifiedLeaves.size();
				sort(specifiedLeaves.begin(), specifiedLeaves.end());
			}
		} else if (s == "-i") {
			ignoreMultipleSubtrees = true;
		} else if (s == "-n" or s == "-a") {
			extractSubtree = true;
		} else if (s == "-T") {
			reprintNewick = true;
		} else if (s == "-v") {
			negation = true;
		} else {
			k = atoi(argv[i]);
		}
	}

	if (k < 1) {
		printUsage();
		return (1);
	}

	scan();
}

void processHeader(const int sampleSize, const int samples, const char* header) {
	string line(header);
	int i = line.find(' ') + 1;
	int j = line.find(' ', i);
	if (extractSubtree) {
		cout << line.substr(0, i) << (complementTree ? sampleSize - k : k) << line.substr(j);
	} else {
		cout << line << endl;
	}
}

void processTree(const char* c) {
	string line(c);
	subtreeFound = getSubtrees(line, k, subtrees);
	if (subtreeFound && ignoreMultipleSubtrees) {
		string rs = subtrees.at(rand() % subtrees.size());
		subtrees.clear();
		subtrees.push_back(rs);
	}
	if (subtreeFound && specifiedLeaves.size() > 0) {
		for (int i = subtrees.size() - 1; i >= 0; i--) {
			getLeavesInTree(subtrees.at(i), subtreeLeaves);
			sort(subtreeLeaves.begin(), subtreeLeaves.end());
			if (subtreeLeaves != specifiedLeaves) {
				subtrees.erase(subtrees.begin() + i);
			}
		}
		subtreeFound = subtrees.size() > 0;
	}
	newick = line;
	if (reprintNewick)
		cout << newick << endl;
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** polymorphisms) {
	if (extractSubtree) {
		if (subtreeFound) {
			for (unsigned int s = 0; s < subtrees.size(); s++) {
				if (!printAll || s > 0) {
					if (reprintNewick)
						cout << newick << endl;
					cout << segsites << endl;
					cout << "positions: ";
					for (int i = 0; i < segsites; i++) {
						cout << " " << positions[i];
					}
					cout << endl;
				}
				for (int j = 0; j < sampleSize; j++) {
					getLeavesInTree(subtrees.at(s), subtreeLeaves);
					bool leafInSubtree = false;
					for (int i = 0; i < subtreeLeaves.size(); i++) {
						if (subtreeLeaves.at(i) == j + 1) {
							leafInSubtree = true;
							break;
						}
					}
					if (complementTree != leafInSubtree) {
						for (int i = 0; i < segsites; i++) {
							cout << polymorphisms[j][i];
						}
						cout << endl;
					}
				}
			}
		}
	} else if (subtreeFound != negation) {
		if (!printAll) {
			//cout << endl << comment << endl;
			if (reprintNewick)
				cout << newick << endl;
			cout << segsites << endl;
			cout << "positions: ";
			for (int i = 0; i < segsites; i++) {
				cout << " " << positions[i];
			}
			cout << endl;
		}
		for (int i = 0; i < sampleSize; ++i) {
			for (int j = 0; j < segsites; j++) {
				cout << polymorphisms[i][j];
			}
			cout << '\n';
		}
	}

}
