/**
 *  This file is part of coatli.
 *
 *  coatli is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  coatli is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTAbLITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with coatli.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  This class computes the cubic terms of the third moments of the unfolded frequency spectrum
 *  under the assumption of constant population size.
 *  The computation uses the unsimplified equations (Supplement A, Part 2) and is very slow!
 *
 *  Copyright 2016 Alexander Klassmann
 */

#include "ntdef.h"

namespace {
double b(int n, int k) {
	if (n == -1 && k == -1) {
		return 1;
	}
	if (k < 0 || n < 0 || k > n) {
		return 0;
	}

	long long b = 1;

	if (k > n - k) {
		k = n - k;
	}

	for (int i = 0; i < k; ++i) {
		b *= n - i;
		b /= i + 1;
	}

	return b;
}

double taa(int n, int h, int i, int j) {
	double s = 0;
	for (int kpp = 2; kpp <= n; kpp++) {
		for (int kp = 2; kp <= kpp; kp++) {
			for (int k = 2; k <= kp; k++) {
				double tmp = 0;
				for (int t = 1; t <= kp - 1; t++) {
					for (int t1 = 0; t1 <= kpp - 2; t1++) {
						for (int t2 = 1; t2 <= kpp - 1; t2++) {
							tmp += b(kp - t - 1, k - 2) / b(kp - 1, k - 1) * t / kp * b(t1 - 1, t - 2)
									* b(kpp - t1 - t2 - 1, kp - t - 1) / b(kpp - 1, kp - 1) * t2 / kpp
									* b(h - i - 1, t1 - 1) * b(i - j - 1, t2 - 2) * b(n - h - 1, kpp - t1 - t2 - 1)
									/ b(n - 1, kpp - 1);
						}
					}
				}
				s += tmp / ((k - 1) * (kp - 1) * (kpp - 1));
			}
		}
	}
	return s;
}

double tab(int n, int h, int i, int j) {
	double s = 0;
	for (int kpp = 2; kpp <= n; kpp++) {
		for (int kp = 2; kp <= kpp; kp++) {
			for (int k = 2; k <= kp; k++) {
				double tmp = 0;
				for (int t = 2; t <= kp - 1; t++) {
					for (int t1 = 1; t1 <= kpp - 2; t1++) {
						for (int t2 = 1; t2 <= kpp - t1; t2++) {
							tmp += b(kp - t - 1, k - 2) / b(kp - 1, k - 1) * t / kp * b(t1 - 1, t - 2)
									* b(kpp - t1 - t2 - 1, kp - t - 1) / b(kpp - 1, kp - 1) * t1 / kpp
									* b(h - i - j - 1, t1 - 2) * b(i - 1, t2 - 1) * b(n - h - 1, kpp - t1 - t2 - 1)
									/ b(n - 1, kpp - 1);
						}
					}
				}
				s += tmp / ((k - 1) * (kp - 1) * (kpp - 1));
			}
		}
	}
	return s;
}

double tba_3(int n, int h, int i, int j) {
	double s = 0;
	for (int kpp = 2; kpp <= n; kpp++) {
		for (int kp = 2; kp <= kpp; kp++) {
			for (int k = 2; k <= kp; k++) {
				double tmp = 0;
				for (int t = 1; t <= kp - 1; t++) {
					for (int t1 = 0; t1 <= kpp - 2; t1++) {
						for (int t2 = 1; t2 <= kpp - t1; t2++) {
							tmp += b(kp - t - 1, k - 2) / b(kp - 1, k - 1) * t / kp * b(t1 - 1, t - 2)
									* b(kpp - t1 - t2 - 1, kp - t - 1) / b(kpp - 1, kp - 1) * (kpp - t1 - t2) / kpp
									* b(h - i - 1, t1 - 1) * b(i - 1, t2 - 1) * b(n - h - j - 1, kpp - t1 - t2 - 2)
									/ b(n - 1, kpp - 1);
						}
					}
				}
				s += tmp / ((k - 1) * (kp - 1) * (kpp - 1));
			}
		}
	}
	return s;
}

double tba_2(int n, int h, int i, int j) {
	double s = 0;
	for (int kpp = 2; kpp <= n; kpp++) {
		for (int kp = 2; kp <= kpp; kp++) {
			for (int k = 2; k <= kp; k++) {
				double tmp = 0;
				for (int t = 1; t <= kp - 1; t++) {
					for (int t1 = 1; t1 <= kpp - 1; t1++) {
						for (int t2 = 1; t2 <= kpp - t1; t2++) {
							tmp += b(kp - t - 1, k - 2) / b(kp - 1, k - 1) * (kp - t) / kp * b(t1 - 1, t - 1)
									* b(kpp - t1 - t2 - 1, kp - t - 2) / b(kpp - 1, kp - 1) * t1 / kpp
									* b(h - j - 1, t1 - 2) * b(i - 1, t2 - 1) * b(n - h - i - 1, kpp - t1 - t2 - 1)
									/ b(n - 1, kpp - 1);
						}
					}
				}
				s += tmp / ((k - 1) * (kp - 1) * (kpp - 1));
			}
		}
	}
	return s;
}

double tba_1(int n, int h, int i, int j) {
	double s = 0;
	for (int kpp = 2; kpp <= n; kpp++) {
		for (int kp = 2; kp <= kpp; kp++) {
			for (int k = 2; k <= kp; k++) {
				double tmp = 0;
				for (int t = 1; t <= kp - 1; t++) {
					for (int t1 = 1; t1 <= kpp - 1; t1++) {
						for (int t2 = 1; t2 <= kpp - t1; t2++) {
							tmp += b(kp - t - 1, k - 2) / b(kp - 1, k - 1) * (kp - t) / kp * b(t1 - 1, t - 1)
									* b(kpp - t1 - t2 - 1, kp - t - 2) / b(kpp - 1, kp - 1) * t2 / kpp
									* b(h - 1, t1 - 1) * b(i - j - 1, t2 - 2) * b(n - h - i - 1, kpp - t1 - t2 - 1)
									/ b(n - 1, kpp - 1);
						}
					}
				}
				s += tmp / ((k - 1) * (kp - 1) * (kpp - 1));
			}
		}
	}
	return s;
}

double tba(int n, int h, int i, int j) {
	return tba_3(n, i, j, h) + tba_2(n, i, h, j) + tba_1(n, h, i, j);
}

double tbb(int n, int h, int i, int j) {
	double s = 0;
	for (int kpp = 2; kpp <= n; kpp++) {
		for (int kp = 2; kp <= kpp; kp++) {
			for (int k = 2; k <= kp; k++) {
				double tmp = 0;
				for (int t = 1; t <= kp - 1; t++) {
					for (int t1 = 1; t1 <= kpp - 2; t1++) {
						for (int t2 = 1; t2 <= kpp - t1; t2++) {
							tmp += b(kp - t - 1, k - 2) / b(kp - 1, k - 1) * (kp - t) / kp * b(t1 - 1, t - 1)
									* b(kpp - t1 - t2 - 1, kp - t - 2) / b(kpp - 1, kp - 1) * (kpp - t1 - t2) / kpp
									* b(h - 1, t1 - 1) * b(i - 1, t2 - 1) * b(n - h - i - j - 1, kpp - t1 - t2 - 2)
									/ b(n - 1, kpp - 1);
						}
					}
				}
				s += tmp / ((k - 1) * (kp - 1) * (kpp - 1));
			}
		}
	}
	return s;
}
}

double tau(int n, int h, int i, int j, int relation) {
	double s = 0;

	if (relation == ANY || relation == NESTED_NESTED) {
		s += taa(n, h, i, j) + taa(n, h, j, i) + taa(n, i, j, h) + taa(n, j, h, i) + taa(n, j, i, h) + taa(n, i, h, j)
				+ tab(n, h, i, j) + tab(n, h, j, i);
	}
	if (relation == ANY || relation == NESTED_DISJOINT) {
		s += tab(n, i, j, h) + tab(n, i, h, j) + tba(n, j, i, h) + tba(n, j, h, i);
	}
	if (relation == ANY || relation == DISJOINT_NESTED) {
		s += tab(n, j, h, i) + tab(n, j, i, h) + tba(n, i, j, h) + tba(n, i, h, j);
	}
	if (relation == ANY || relation == DISJOINT_DISJOINT) {
		s += tba(n, h, i, j) + tba(n, h, j, i) + tbb(n, h, i, j) + tbb(n, i, j, h) + tbb(n, j, h, i) + tbb(n, h, j, i)
				+ tbb(n, j, i, h) + tbb(n, i, h, j);
	}

	return s;
}

void initT3(int n) {
}

void freeT3(int n) {
}
