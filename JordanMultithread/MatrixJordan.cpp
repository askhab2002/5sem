
#include <iostream>
#include <math.h>
#include <thread>
#include <mutex>
#include "MatrixRead.hpp"
#include "MatrixJordan.hpp"

using namespace std;

std::mutex MyMutex;

void Addition(const int &n, const int &threads, const int &t, double *M, double *I);
void OneInvert(const int &n, const int &t, double *M, double *I, const int &threads);
void RowsAdditionThread(const int &size_, const int &start, const int &stop, const int &i, const int &j, const double &h, double *matrix_, double *inverse_);

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

void JordanInverse(const int &n, const int &m,  double *M, double *I, const int &threads) {

	//	int ZeroRow = 0;
	cout << "   Cores = " << std::thread::hardware_concurrency() << endl;
	for(int t = 0; t < n; t++) { 
		OneInvert(n, t, M, I, threads);
		/*
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

		   thread thr1(Addition, 0, n/2, t, M, I);
		   thread thr2(Addition, n/2, n, t, M, I);

		   thr1.join();
		   thr2.join(); */

		/*                for(int v = 0; v < n; v++) {
				  if(v == t) {
				  continue;
				  }
				  RowsAddition(n, v, t, -M[v * n + t], M, I); 
				  }   */


	}

	cout << "-----------Обратная матрица:------------" << endl;
	Print(n, n, m, I);



	cout << "-----------Метод Жордана закончился успешно------------------" << endl;

}

void Addition(const int &n, const int &threads, const int &t, double *M, double *I) {
	//	lock_guard<mutex> guard(MyMutex);
	double h = 0;
	int start = 0;
	int stop = 0;

	for(int i = 0; i < n; i++) {
		if(i == t) {
			continue;
		}
//		cout << " От строки " << i << " отнимается строка " << t << " thread: " << std::this_thread::get_id() << endl;

		//		RowsAddition(n, i, t, -M[i * n + t], M, I);
		h = M[i * n + t];
		/*		for(int k = 0; k < n; k++) {

				M[n * i + k] -= h * M[n * t + k];
				I[n * i + k] -= h * I[n * t + k];
				} */

		thread *Threads_ = new thread[threads];
		for(int j = 0; j < threads - 1; j++) {

			start = j * (n/threads);
			stop = (j + 1) * (n/threads);
			cout << "  start = " << start << "  stop = " << stop << endl;

			Threads_[j] = thread(RowsAdditionThread, cref(n), cref(start), cref(stop), cref(i), cref(t), cref(h), M, I);
		}
		cout << "  start = " << stop << "  stop = " << n << endl;
		Threads_[threads - 1] = thread(RowsAdditionThread, cref(n), cref(stop), cref(n), cref(i), cref(t), cref(h), M, I);

		for(int j = 0; j < threads; j++)  {
			Threads_[j].join();
		}

		delete [] Threads_;

	} 

}

void RowsAdditionThread(const int &size_, const int &start, const int &stop, const int &i, const int &j, const double &h, double *matrix_, double *inverse_) {
	for(int k = start; k < stop; k++) {
		matrix_[size_ * i + k] = h * matrix_[size_ * j + k] + matrix_[size_ * i + k];
		inverse_[size_ * i + k] = h * inverse_[size_ * j + k] + inverse_[size_ * i + k];
	}

}

void OneInvert(const int &n, const int &t, double *M, double *I, const int &threads) {
	int ZeroRow = threads;
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

	Addition(n, threads, t, M, I);

/*
	thread *Threads_ = new thread[threads];

	int start = 0;
	int stop = 0;

	for(int i = 0; i < threads - 1; i++) {

		start = i * (n/threads);
		stop = (i + 1) * (n/threads);
		cout << "  start = " << start << "  stop = " << stop << endl;

		Threads_[i] = thread(Addition, cref(n), cref(start), cref(stop), cref(t), M, I);
	}
	cout << "  start = " << stop << "  stop = " << n << endl;
	Threads_[threads - 1] = thread(Addition, cref(n), cref(stop), cref(n), cref(t), M, I);

	for(int i = 0; i < threads; i++)  {
		Threads_[i].join();
	}

	delete [] Threads_;  */
	/*
	   thread Threads[3];
	   int start, stop;

	   for(int i  = 0; i < 2; i++) {
	   start = i * (n/3);
	   stop = (i + 1) * (n/3);
	   Threads[i] = thread(Addition, cref(n), cref(start), cref(stop), cref(t), M, I);
	   }
	   Threads[2] = thread(Addition, cref(n), cref(stop), cref(n), cref(t), M, I);
	   for(int i = 0; i < 3; i++) {
	   Threads[i].join();
	   }  */
}

double Norm(double *M, double *I, const int &n) {
	double *K = new double[n * n];

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			K[i * n + j] = 0;
			for(int k = 0; k < n; k++) {
				K[i * n + j] = K[i * n + j] + M[i * n + k] * I[k * n + j];
			}
		}
	}

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
