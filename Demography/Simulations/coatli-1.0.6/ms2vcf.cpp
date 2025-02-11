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
 *  Translates 'ms'-output into vcf
 *
 *  Created on: 03.07.20
 *      Author: Alexander Klassmann
 */
#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

void scan();

const int default_ploidy = 2;
const int default_length = 1000000;
const string default_variantbase = ".";
const string default_chrombase = "";
const string default_samplebase = "HG";
const char *alphabet[] = { "A", "C", "G", "T" };
string variantbase = default_variantbase;
string chrombase = default_chrombase;
string samplebase = default_samplebase;
int ploidy = default_ploidy;
int length = default_length;
int alphabet_size = 4;
int current_segsite = 0;
int current_chrom = 0;
int number_removed = 0;

void processHeader(const int sampleSize, const int samples, const char *header) {
	if (sampleSize % ploidy != 0) {
		cerr << "Warning: sample size is not a multiple of ploidy!" << endl;
	}

	cout << "##fileformat=VCFv4.2" << endl;
	cout << "##source=ms2vcf" << endl;
	cout << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << endl;
	for (int i = 1; i <= samples; i++) {
		cout << "##contig=<ID=" << chrombase << i << ",length=" << length << ">" << endl;
	}
	cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;

	cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

	int number_preceeding_zeros = floor(log10(sampleSize / ploidy)) + 1;

	for (int j = 0; j < sampleSize / ploidy; j++) {
		cout << "\t" << samplebase << setfill('0') << setw(number_preceeding_zeros) << (j + 1);
	}
	cout << endl;
}

void processTree(const char*) {
}

void processSample(const int sampleSize, const int segsites, const double *positions, bool **polymorphisms) {
	int lastPosition = -1;
	current_chrom++;

	for (int i = 0; i < segsites; i++) {
		current_segsite++;

		int currentPosition = (int) (positions[i] * length + 0.5);
		if (lastPosition == currentPosition) {
			number_removed++;
			continue;
		}
		lastPosition = currentPosition;

		/* random indices for REF and ALT alleles */
		int REF = rand() % alphabet_size;
		int ALT = (REF + rand() % (alphabet_size - 1) + 1) % alphabet_size;

		cout << chrombase << current_chrom << "\t" << currentPosition << "\t" << variantbase;

		/* if there is a non-trivial id base, append number */
		if (variantbase != ".") {
			cout << current_segsite;
		}
		cout << "\t" << alphabet[REF] << "\t" << alphabet[ALT] << "\t.\tPASS\t.\tGT";

		for (int j = 0; j < sampleSize / ploidy; j++) {
			cout << "\t";
			for (int k = 0; k < ploidy; k++) {
				cout << polymorphisms[j * ploidy + k][i];
				if (k < ploidy - 1) {
					cout << "|";
				}
			}
		}
		cout << endl;
	}
}

void printUsage() {
	cout << "Usage: ms2vcf" << endl;
	cout << "Options:" << endl;
	cout << " -length [integer] length of chromosome (positions of ms are multiplied by this number) (default "
			<< default_length << ")" << endl;
	cout << " -ploidy [integer] ploidy of individuals (default " << default_ploidy << ")" << endl;
	cout << " -chrombase [string] base name for chromosomes (default \"" << default_chrombase << "\")" << endl;
	cout << " -samplebase [string] base name for sample individuals (default \"" << default_samplebase << "\")" << endl;
	cout << " -variantbase [string] base name for variant ids, (default \"" << default_variantbase << "\")" << endl;
	cout << " -h print options" << endl;
}

int main(int argc, char *argv[]) {
	for (int i = 1; i < argc; i++) {
		string s = string(argv[i]);
		if (s == "-ploidy" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> ploidy;
		} else if (s == "-length" && argc >= i) {
			i++;
			stringstream ss(argv[i]);
			ss >> length;
		} else if (s == "-chrombase" && argc >= i) {
			i++;
			chrombase = argv[i];
		} else if (s == "-samplebase" && argc >= i) {
			i++;
			samplebase = argv[i];
		} else if (s == "-variantbase" && argc >= i) {
			i++;
			variantbase = argv[i];
		} else if (s == "-h") {
			printUsage();
			return (0);
		} else {
			cout << "Error: unrecognized option '" << s << "'" << endl;
			printUsage();
			return (1);
		}
	}

	if (length < 1) {
		cout << "Error: chromosome length cannot be smaller than 1." << endl;
		return (1);
	}
	if (ploidy < 1) {
		cout << "Error: ploidy cannot be smaller than 1." << endl;
		return (1);
	}

	std::ios_base::sync_with_stdio(true);

	srand (time(NULL));

	scan ();

	if (number_removed > 0) {
		cerr << "Warning: removed " << number_removed << " sites yielding a position already used in vcf." << endl;
		cerr << "Please consider increasing chromosome length in vcf and/or ms output precision!" << endl;
	}

	return 0;
}
