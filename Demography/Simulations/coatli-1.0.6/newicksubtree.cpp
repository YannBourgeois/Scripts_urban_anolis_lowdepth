/*
 *  newicksubtree.cpp
 *
 *  Created on: 20.10.2013
 *      Author: Alexander Klassmann
 *
 *  Isolates subtree-strings of a specified size from newick string
 *
 *
 *  remove out-comment of main function if running as standalone tool!
 *
 *  compile: g++ newicksubtree.cpp -o newicksubtree
 *
 *  example: echo "(A,(B,C))" | ./newicksubtree 2
 *
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace std;

const char lp = '(', rp = ')', sep = ',', dp = ':';

/**
 * fills vector with all substrees of size k, represented by substrings of the original string
 */
bool getSubtrees(const string newick, const unsigned int k, vector<string>& subtrees) {
	unsigned int i, j, lpc, rpc, nc;
	char last = '^';
	string ss;
	subtrees.clear();
	for (i = 0; i < newick.length(); i++) {
		ss = newick.substr(i);
		lpc = 0;
		rpc = 0;
		nc = 0;
		j = 0;
		do {
			if (ss.at(j) == lp) {
				lpc++;
				last = lp;
			} else if (ss.at(j) == rp) {
				rpc++;
				if (last == sep)
					nc++;
				last = rp;
			} else if (ss.at(j) == sep) {
				if (last == lp)
					nc++;
				last = sep;
			}
			j++;
		} while ((lpc > rpc) && j < ss.length());
		if (lpc == rpc && nc == k) {
			subtrees.push_back(ss.substr(0, j));
		}
	}
	return subtrees.size() > 0;
}

/**
 * returns random subtree; empty string, if no subtree is found
 */
string getRandomSubtree(const string newick, const unsigned int k) {
	vector<string> subtrees;
	getSubtrees(newick, k, subtrees);
	return subtrees.size() > 0 ? subtrees.at(rand() % subtrees.size()) : "";
}

/**
 * get terminal nodes out of string;
 * Assumes that nodes are coded as integers!
 */
void getLeavesInTree(const string newick, vector<unsigned int> &leaves) {
	leaves.clear();
	unsigned int i = 1, j = 0;
	while (i < newick.length()) {
		char c = newick.at(i);
		if (c != lp && c != rp && c != sep && c != dp
				&& (newick.at(i - 1) == lp || newick.at(i - 1) == sep)) {
			j = newick.find(dp, i);
			leaves.push_back((unsigned int) atoi(newick.substr(i, j - i).c_str()));
			if (j > 0)
				i = j + 6;
		} else
			i++;
	}
}

/*
 int main(int argv, char *argc[]){
 string line;
 vector<string> subtrees;
 unsigned int k = atoi(argc[1]);
 getline(cin, line);
 getSubtrees(line,k,subtrees);
 for(int i=0;i<subtrees.size();i++){
 cout << subtrees.at(i)<<endl;
 }
 }
 */
