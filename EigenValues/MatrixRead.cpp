#include <iostream>
#include <fstream>
#include <math.h>
#include "MatrixRead.hpp"

using namespace std;

double *MatrixRead(const int &n, const int &m, const int &k, const string &File, int *correctness) {

        if((n <= 0) || (n < m) || (m <= 0)) {
		*correctness = -1;
		return NULL;
	}
	
	if(k > 0 && k < 5) {
		double *Matrix = new double[n * n];
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				Matrix[j + i * n] = f(k, n, i + 1, j + 1);
			}
		}

		*correctness = 1;
		
                  
		return Matrix;
	}
        
	if(k != 0) {
		*correctness = -1;

		return NULL;
	}

        ifstream input(File);
	if(!(input.is_open())) {
		*correctness = -1;
		return NULL;
	}

	double *Matrix = new double[n * n];
	double value = 0;
        int size = 0;

	for(; input >> value; size++) {
		if(size > n * n) {
			delete [] Matrix;
                        *correctness = -2;
			input.close();
			cout << " ____________________" << endl;
			return NULL;
		}
		Matrix[size] = value;
	}

	if(size != n * n) { 
		delete [] Matrix;
		*correctness = -2;
		input.close();
		return NULL;
	}

	input.close();
/*
	for(int i = 0; i < m; i++) {
		for(int j = 0; j < m; j++) {
			cout << std::uppercase << std::scientific << Matrix[j + n * i] << " ";
		}
		cout << endl;
	}
*/        
	*correctness = 1;
	return Matrix;
}

double f(const int &k, const int &n,const int &i, const int &j) {
	switch(k) {
		case 1:
			return (n - max(i, j) + 1);
		case 2:
			if(i == j) {
				return 2;
			}
			if(fabs(fabs(i - j) - 1) < 1e-4) {
				return -1;
			}
			return 0;
		case 3:
			if(i == j && i < n) {
				return 1;
			}
		        if(i == n) {
			        return j;
			}
		        if(j == n) {
			        return i;
			}
		        return 0;	
		case 4:
			return 1/((double)i + (double)j - 1);
		default:
			return -1;
	}
}
