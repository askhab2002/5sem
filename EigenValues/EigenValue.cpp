#include <iostream>
#include <math.h>
#include "SignChange.hpp"

using namespace std;

double EigenValue(const int &n, double *M, const int &k, const double &eps) {
	double max = 0;
	double current = 0;

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			current += fabs(M[i * n + j]);
		}

		if(max < current) {
			max = current;
		}
		current = 0;
	}
 /*
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(fabs(M[i * n + j]) > max) {
				max = fabs(M[i * n + j]);
			}
		}
	}  */

	double a, b, c;

	a = -max;
	b = max;
//        cout << "   max = " << max << endl;

//	int Sign = 0;

	for(int i = 0; b - a > eps; i++) {
		c = (b + a)/2;
		if( i < 10) {
//			cout << " c = " << c << endl;
//                        Sign = SignChange(n, M, c);
//			cout << " SignChange = " << Sign << endl;
		}

		if(SignChange(n, M, c) < k) {
			a = c;
		}
		else {
			b = c;
		}
//		cout << "   a = " << a << "   b = " << b << endl;
	}

	return (b + a)/2;
}
