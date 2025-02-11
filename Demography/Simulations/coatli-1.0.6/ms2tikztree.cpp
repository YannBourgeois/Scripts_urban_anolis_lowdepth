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
 *  Created on: 02.01.2015
 *      Author: Alexander Klassmann
 *
 *  Draws a latex-tikz tree from the output of 'ms', if run with option '-T' (include tree structure)
 *
 */

#include "Node.h"
#include <map>
#include <list>
#include <sstream>
#include <iostream>

void scan();
Node* getTree(const string line);
void tree2tikz(Node* node, bool complete, bool labels, double, double);

Node* root = NULL;
bool complete = false;
bool labels = false;
double scale = 1.;
double rscale = 1.;

void printUsage() {
	cout << "Usage: cat ms.output | ms2tikztree [-n][-r][-s][-l] " << endl;
	cout << "Options:" << endl;
	cout << "-c complete (standalone) latex document" << endl;
	cout << "-r scale radius of mutation balls (default 1)" << endl;
	cout << "-s scale size of tree (default 1; yields picture of about 10 cm width)" << endl;
	cout << "-l print node labels" << endl;
	cout << "Note: 'ms' has to be run with option -T (print tree structure) and recombination rate 0!" << endl;
}

void getDescendentLeaves(Node* node, list<int>* leaves) {
	if (node->children.size() == 0) {
		int number;
		stringstream ss(node->name);
		ss >> number;
		leaves->push_back(number);
	} else {
		for (int i = 0; i < node->children.size(); i++) {
			getDescendentLeaves(node->children.at(i), leaves);
		}
	}
}

void constructMap(Node* node, map<Node*, list<int>*>* nodeMap) {
	list<int>* leaves;
	for (int i = 0; i < node->children.size(); i++) {
		leaves = new list<int>;
		getDescendentLeaves(node->children.at(i), leaves);
		leaves->sort();
		nodeMap->insert(pair<Node*, list<int>*>(node->children.at(i), leaves));
		constructMap(node->children.at(i), nodeMap);
	}
}

void addMutationsToTree(Node* node, int sampleSize, int segsites, bool** s) {
	list<int>* leaves;
	map<Node*, list<int>*> nodeMap;
	constructMap(node, &nodeMap);
	for (int j = 0; j < segsites; j++) {
		list<int> mutatedLeaves;
		for (int i = 0; i < sampleSize; i++) {
			if (s[i][j]) {
				mutatedLeaves.push_back(i + 1);
			}
		}
		for (std::map<Node*, list<int>*>::iterator it = nodeMap.begin(); it != nodeMap.end(); ++it) {
			leaves = it->second;
			if (*leaves == mutatedLeaves) {
				it->first->mutations++;
			}
		}
	}
}

void processHeader(const int sampleSize, const int samples, const char* header) {
}

void processTree(const char* newick) {
	string s(newick);
	if (s.length() == 0) {
		printUsage();
		return;
	}
	root = getTree(s);
}

void processSample(const int sampleSize, const int segsites, const double* positions, bool** s) {
	if (root != NULL) {
		addMutationsToTree(root, sampleSize, segsites, s);
	} else {
		printUsage();
	}
}

int main(int argc, char *argv[]) {
	for (int i = 1; i < argc; i++) {
		string s(argv[i]);
		if (s == "-l") {
			labels = true;
		} else if (s == "-c") {
			complete = true;
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
	if (rscale <= 0 || scale <= 0) {
		printUsage();

	}

	scan();

	if (root != NULL) {
		tree2tikz(root, complete, labels, scale, rscale);
	} else {
		cout << "Error: could not read tree." << endl;
	}

	return 0;
}

