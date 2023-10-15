#include <iostream>
#include <math.h>

using namespace std;

class Matrix {
	public:
		int size_ = 0;
		double *matrix_ = NULL;
		double *inverse_ = NULL;
	public:
		Matrix() {
			size_ = 0;
			matrix_ = NULL;
			inverse_ = NULL;
		}
		Matrix(const int &n, const int &m, const int &k, const string &File) {
			int correctness = 0;
			matrix_ = MatrixRead(n, m, k, File, &correctness);

			inverse_ = new double[n * n];
			for(int i = 0; i < n * n; i++) {
				inverse_[i] = 0;
			}
			for(int i = 0; i < n; i++) {
				inverse_[i * n + i] = 1;
			}

			size_ = n;

			if(correctness != 1) {
				size_ = 0;
				matrix_ = NULL;
			}
		}
		Matrix(const Matrix& M) {
			size_ = M.size_;
			matrix_ = new double[size_ * size_];
			for(int k = 0; k < size_ * size_; k++) {
				matrix_[k] = M.matrix_[k];
			}

			inverse_ = new double[size_ * size_];
                        for(int i = 0; i < size_ * size_; i++) {
                                inverse_[i] = M.inverse_[i];
                        }
                        
	        }	
		~Matrix() {
			delete [] matrix_;
			delete [] inverse_;
		}
                
		void RowsAddition(const int &i, const int &j, const double &h) {
		        for(int k = 0; k < size_; k++) {
				matrix_[size_ * i + k] = h * matrix_[size_ * j + k] + matrix_[size_ * i + k];
				inverse_[size_ * i + k] = h * inverse_[size_ * j + k] + inverse_[size_ * i + k];
			//	cout << h << " * " << matrix_[size_ * j + k] << " + " << matrix_[size_ * i + k] << " ";
			}
			cout << endl;
		}
		void RowMultiply(const int &i, const double &h) {
		        for(int k = 0; k < size_; k++) {
		                matrix_[size_ * i + k] = h * matrix_[size_ * i + k];
				inverse_[size_ * i + k] = h * inverse_[size_ * i + k];
		        }
	        }	       
		void RowsSwap(const int &i, const int &j) { 
			for(int k = 0; k < size_; k++) {
				swap(matrix_[size_ * i + k], matrix_[size_ * j + k]);
				swap(inverse_[size_ * i + k], inverse_[size_ * j + k]);
			}
		}
		int NonZero(const int &j) {
			for(int k = j; k < size_; k++) {
				if(matrix_[size_ * k + j] != 0) {
					return k + 1;
				}
			}
			return -1;
		}
                
		double Get(const int &i, const int &j) {
			return matrix_[(i - 1) * size_ + (j - 1)];
		}
		void Add(const int &i, const int &j, double h) {
			matrix_[(i - 1) * size_ + (j - 1)] = h;
		}
		void Print() {
			for(int i = 0; i < size_; i++) {
				for(int j = 0; j < size_; j++) {
					cout << matrix_[i * size_ + j] << " ";
				}
				cout << endl;
			}
		}
		void PrintInverse() {
                        for(int i = 0; i < size_; i++) {
                                for(int j = 0; j < size_; j++) {
                                        cout << inverse_[i * size_ + j] << " ";
                                }
                                cout << endl;
                        }
                }
		 
};
