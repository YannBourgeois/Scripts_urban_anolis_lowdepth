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
 *  Draws tikz tree from output of 'ms'.
 *
 */

#include <iostream>
#include <iomanip>
#include "Node.h"
#include "tikzdef.h"

int getSize(Node* node) {
	int s = 0;
	if (node->children.size() == 0) {
		return 1;
	} else {
		for (int i = 0; i < node->children.size(); i++) {
			s += getSize(node->children.at(i));
		}
	}
	return s;
}

double getHeight(Node* node) {
	double height = 0;
	while (node->parent != NULL) {
		height += node->length;
		node = node->parent;
	}
	return height;
}

double getMaxHeight(Node* node) {
	if (node->children.size() == 0) {
		return getHeight(node);
	} else {
		double max = 0;
		for (int i = 0; i < node->children.size(); i++) {
			double h = getMaxHeight(node->children.at(i));
			if (h > max) {
				max = h;
			}
		}
		return max;
	}
}

void printTree(Node* node, int level) {
	cout << level << " : " << node->id << " " << node->name << "," << node->length << " " << getSize(node) << " "
			<< getHeight(node) << endl;
	for (int i = 0; i < node->children.size(); i++) {
		printTree(node->children.at(i), level + 1);
	}
}

void printNodes(Node* node, int position, int leaves, double treeHeight) {
	cout << "\\coordinate (" << node->id << ") at (" << (position + (getSize(node) - 1) / 2.0) / (double) (leaves - 1)
			<< "," << treeHeight - getHeight(node) << ");" << endl;
	int s = position;
	for (int i = 0; i < node->children.size(); i++) {
		printNodes(node->children.at(i), s, leaves, treeHeight);
		s += getSize(node->children.at(i));
	}
}

void printLabels(Node* node, int position, int leaves, double treeHeight) {
	if (node->name.length() > 0) {
		cout << "\\node[below] at (" << node->id << ") {$" << node->name << "$};" << endl;
	}
	int s = position;
	for (int i = 0; i < node->children.size(); i++) {
		printLabels(node->children.at(i), s, leaves, treeHeight);
		s += getSize(node->children.at(i));
	}
}

void printPath(Node* node) {
	if (node->children.size() > 0) {
		cout << "(" << node->id << ")";
		for (int i = 0; i < node->children.size(); i++) {
			cout << " edge ";
			for (int j = 0; j < node->children.at(i)->mutations; j++) {
				cout << "node [pos=" << (j + 1) / (double) ((node->children.at(i)->mutations) + 1)
						<< ",style=mutation] {} ";
			}
			cout << "(" << node->children.at(i)->id << ")" << endl;
		}
		for (int i = 0; i < node->children.size(); i++) {
			printPath(node->children.at(i));
		}
	}
}

void printPreamble(bool header, double scale, double rscale) {
	cout << fixed << setprecision(precision);
	if (header) {
		cout << "\\documentclass{standalone}" << endl;
		cout << "\\usepackage{tikz}" << endl;
		cout << "\\begin{document}" << endl;
	}
	cout << "\\tikzstyle{mutation}=[circle,shading=ball,ball color=" << ballcolor << ",scale=" << rscale << "]" << endl;
	cout << "\\begin{tikzpicture}[scale=" << 10 * scale << ",color=" << linecolor << ",thick]" << endl;
}

void printEnd(bool header) {
	cout << "\\end{tikzpicture}" << endl;
	if (header) {
		cout << "\\end{document}" << endl;
	}
}

void tree2tikz(Node* node, bool header, bool labels, double scale, double rscale) {
	printPreamble(header, scale, rscale);
	printNodes(node, 0, getSize(node), getMaxHeight(node));
	cout << "\\path" << endl;
	printPath(node);
	cout << ";" << endl;
	if (labels) {
		printLabels(node, 0, getSize(node), getMaxHeight(node));
	}
	printEnd(header);
}

