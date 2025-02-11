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
 *      Author: klasal
 *
 *  Parses newick Tree output from ms
 *
 */

#include <sstream>
#include "Node.h"

const char lp = '(', rp = ')', comma = ',', dp = ':';

/**
 * parse newick string
 */
Node* getTree(const string newick) {
	int j = 0, lastPos = 0, nodeNumber = 0;
	char lastChar = '^';
	stringstream ss;

	Node* rootNode = NULL;
	Node* currentParentNode = NULL;
	Node* currentNode = NULL;
	Node* newNode;

	do {
		switch (newick.at(j)) {
		case lp:
			newNode = new Node();
			nodeNumber++;
			ss << nodeNumber;
			ss >> newNode->id;
			ss.clear();
			if (rootNode == NULL) {
				rootNode = newNode;
				rootNode->name = "";
			} else {
				newNode->parent = currentNode;
				currentNode->children.push_back(newNode);
				currentParentNode = currentNode;
			}
			currentNode = newNode;
			lastChar = lp;
			lastPos = j;
			break;
		case rp:
			ss << newick.substr(lastPos + 1, j - lastPos - 1);
			ss >> currentNode->length;
			ss.clear();
			if (currentNode->parent == rootNode) {
				return rootNode;
			}
			currentNode = currentParentNode;
			currentParentNode = currentParentNode->parent;
			lastChar = rp;
			lastPos = j;
			break;
		case dp:
			currentNode->name = newick.substr(lastPos + 1, j - lastPos - 1);
			lastChar = dp;
			lastPos = j;
			break;
		case comma:
			ss << newick.substr(lastPos + 1, j - lastPos - 1);
			ss >> currentNode->length;
			ss.clear();
			currentNode = currentParentNode;
			currentParentNode = currentParentNode->parent;
			lastChar = comma;
			lastPos = j;
			break;
		default:
			if (lastPos == j - 1) {
				if (lastChar == lp || lastChar == comma) {
					newNode = new Node();
					nodeNumber++;
					ss << nodeNumber;
					ss >> newNode->id;
					ss.clear();
					newNode->parent = currentNode;
					currentNode->children.push_back(newNode);
					currentParentNode = currentNode;
					currentNode = newNode;
				}
			}
		}
		j++;
	} while (j < newick.length());
	return NULL;
}

