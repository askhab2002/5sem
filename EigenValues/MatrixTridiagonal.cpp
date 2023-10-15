
#include <math.h>
#include <iostream>

using namespace std;

void UMultiply(double *A, const int &n, const int &k, double *x);

void Print(const int &n, const int &l, const int &m, double *matrix_) { //Печать 
        cout << endl;
        for(int i = 0; i < min(n, m); i++) {
                for(int j = 0; j < min(l, m); j++) {
                        cout << std::uppercase << std::scientific << matrix_[i * n + j] << " ";
                }
                cout << endl;
        }
}

double *Multiply(double *A, double *B, const int &n, const int &m, const int &k, const int &l) { // Ненужная функция, неиспользуется
	if(m != k) {
		return NULL;
	}

	double *C = new double[n * l];
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < l; j++) {
			C[i * n + j] = 0;
			for(int t = 0; t < m; t++) {
				C[i * n + j] += A[i * m + t] * B[t * l + j];
			}
		}
	} 
	return C;
}

double *X(const int &n, double *a) {// Преобразуем столбец матрицы в вектор X, который породит матрицу отражений
	double *b = new double[n];
	for(int i = 0; i < n; i++) {
		b[i] = a[i];
	}

	double a_ = 0;
	double sum = 0;
        for(int i = 0; i < n; i++) {
                sum += a[i] * a[i];
	}

	a_ = sqrt(sum);
	 
	a[0] = a[0] - a_;
        
	sum = 0;
	for(int i = 0; i < n; i++) {
		sum += a[i] * a[i];
	}
	a_ = sqrt(sum);
	

	for(int i = 0; i < n; i++) {
		a[i] = - (a[i])/(a_);
	}  

        delete [] b;

	return a;
}

void Tridiagonal(const int &n, double *B) { //Тридиагонализация
        double *A = new double[n * n];
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			A[i * n + j] = B[i * n + j];
		}
	}

	for(int k = 1; k <= n - 2; k++) { 

                double *a = new double[n - k];
		for(int i = 0; i < n - k; i++) {
			a[i] = A[(i + k) * n + k - 1];
		}
                
		double *x;
		x = X(n - k, a);
		if(isnan(x[0])) {
			continue;
		}

                UMultiply(A, n, k, x);            
                
		delete [] a;
	}
        
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) { 
			if(fabs(A[i * n + j]) < 1e-5) {
				B[i * n + j] = 0;
				continue;
			}
			B[i * n + j] = A[i * n + j];
		}
	}

	delete [] A;
}

void UMultiply(double *A, const int &n, const int &k, double *x) { //Умножение матрицы на матрицу отражений
      
	double *y = new double[n - k];
	for(int i = 0; i < n - k; i++) {
		y[i] = 0;
		for(int j = 0; j < n - k; j++) {
                        y[i] += A[(i + k) * n + (j + k)] * x[j];
		}
	}

        double *z = new double[n - k];

	double sum = 0;
	for(int j = 0; j < n - k; j++) {
                sum += x[j] * y[j];
        }

	for(int i = 0; i < n - k; i++) {
		z[i] = 2 * y[i] - 2 * sum * x[i];
        }
        
	double row = 0;
	for(int i = k; i < n; i++) {
		row += fabs(A[(k - 1) * n + i]) * fabs(A[(k - 1) * n + i]);
	}
	row = sqrt(row);
        
        A[k * n + (k - 1)] = row;
	A[(k - 1) * n + k] = row;

	for(int i = k + 1; i < n; i++) {
		A[i * n + (k - 1)] = 0;
		A[(k - 1) * n + i] = 0;
	}
        
	for(int i = k; i < n; i++) {
		for(int j = k ; j < n; j++)  {
	                A[i * n + j] -= x[i - k] * z[j - k] + x[j - k] * z[i - k];
	        }
        }

	delete [] y;
	delete [] z;

}





