
#include <iostream>
#include <math.h>
#include "MatrixRead.hpp"
#include "MatrixJordan.hpp"

using namespace std;

void RowsAddition(const int &size_, const int &i, const int &j, const double &h, double *matrix_, double *inverse_) {
	for(int k = 0; k < size_; k++) {
		matrix_[size_ * i + k] = h * matrix_[size_ * j + k] + matrix_[size_ * i + k];
		inverse_[size_ * i + k] = h * inverse_[size_ * j + k] + inverse_[size_ * i + k];
	}
	
}
void RowMultiply(const int &size_, const int &i, const double &h, double *matrix_, double *inverse_) {
	for(int k = 0; k < size_; k++) {
		matrix_[size_ * i + k] = h * matrix_[size_ * i + k];
		inverse_[size_ * i + k] = h * inverse_[size_ * i + k];
	}
}
void RowsSwap(const int &size_, const int &i, const int &j, double *matrix_, double *inverse_) {
	for(int k = 0; k < size_; k++) {
		swap(matrix_[size_ * i + k], matrix_[size_ * j + k]);
		swap(inverse_[size_ * i + k], inverse_[size_ * j + k]);
	}
}
int NonZero(const int &size_, const int &j, double *matrix_) {
	for(int k = j; k < size_; k++) {
		if(fabs(matrix_[size_ * k + j] - 0) >  1e-5) {
			return k + 1;
		}
	}
	return -1;
}

void Print(const int &n, const int &l, const int &m, double *matrix_) {
	cout << endl;
	for(int i = 0; i < min(n, m); i++) {
		for(int j = 0; j < min(l, m); j++) {
			cout << std::uppercase << std::scientific << matrix_[i * n + j] << " ";
		}
		cout << endl;
	}
}

void JordanInverse(const int &n, const int &m,  double *M, double *I) {
        
        int ZeroRow = 0;
        for(int t = 0; t < n; t++) {
                ZeroRow = NonZero(n, t, M); 
                if(ZeroRow == -1) {
                        I = NULL;
                        cout << "    Матрица вырождена" << endl;
                        return;
                }
                if(t != ZeroRow - 1) {
                        RowsSwap(n, t, ZeroRow - 1, M, I); 
                }
                RowMultiply(n, t, (1/M[t * n + t]), M, I);
                for(int v = 0; v < n; v++) {
                        if(v == t) {
                                continue;
                        }
                        RowsAddition(n, v, t, -M[v * n + t], M, I); 
                }
        }
        
	cout << "-----------Обратная матрица:------------" << endl;
        Print(n, n, m, I);

       

        cout << "-----------Метод Жордана закончился успешно------------------" << endl;
        
}

double Norm(double *M, double *I, int n) {
	double *K = new double[n * n];

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			K[i * n + j] = 0;
			for(int k = 0; k < n; k++) {
				K[i * n + j] = K[i * n + j] + M[i * n + k] * I[k * n + j];
			}
		}
	}
//        Print(n, n, n, M);
//	Print(n, n, n, I);
//	Print(n, n, n, K);

	for(int i = 0; i < n; i++) {
		K[i * n + i]--;
	}
        
	double sum = 0;
	double max = 0;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			sum += fabs(K[i * n + j]);
		}
                if(i == 0) {
			max = sum;
		}

		if(sum > max) {
			max = sum;
		}
		sum = 0;
	}

	delete [] K;
	cout << "    Norm = " << std::uppercase << std::scientific << max << endl;
	return max;
}
