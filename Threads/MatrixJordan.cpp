
#include <iostream>
// #include <barrier>
#include <atomic>
#include <pthread.h>
#include <chrono>
#include <math.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include "MatrixRead.hpp"
#include "MatrixJordan.hpp"

using namespace std;

//std::atomic_bool stop_threads = false;
///int stop_threads = false;

void RowsAddition(const int &size_, const int &i, const int &j, const double &h, double *matrix_, double *inverse_) {
	for(int k = 0; k < size_; k++) {
		matrix_[size_ * i + k] = h * matrix_[size_ * j + k] + matrix_[size_ * i + k];
		inverse_[size_ * i + k] = h * inverse_[size_ * j + k] + inverse_[size_ * i + k];
	}

}

void RowsAddition_(const int &size_, const int &start, const int &stop, const int &i, const int &j, const double &h, double *matrix_, double *inverse_) {

	for(int k = start; k < stop; k++) {
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
 

void synchronize(int total_threads)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;

	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);
}



void JordanInverse(const int n, const int m,  double *M, double *I, const int start, const int stop, barrier *my_barrier) {

	int ZeroRow = m;
	ZeroRow = 0;

	for(int t = 0; t < n; t++) {

		if(start == 0) {
			ZeroRow = NonZero(n, t, M); 
			if(ZeroRow == -1) {
				I = NULL;
				cout << "    Матрица вырождена" << endl;
			        //myBarrier->threadsWaiting.store(0);
                                //myBarrier->waitVariable.notify_all();
				//flag
///				my_barrier->waitVariable.notify_all();
///				stop_threads = true;
				return;
			}
			if(t != ZeroRow - 1) {
				RowsSwap(n, t, ZeroRow - 1, M, I); 
			}
			RowMultiply(n, t, (1/M[t * n + t]), M, I);
		}
//		synchronize(4);
//		flag
///                if(stop_threads) {
///			return;
///		}

		my_barrier->wait();
		for(int v = start; v < stop; v++) {
			if(v == t) {
				continue;
			}
			RowsAddition(n, v, t, -M[v * n + t], M, I); 
		}
//		synchronize(4);
///                if(stop_threads) {
///                        return;
///                }
		my_barrier->wait();

	}       

////	cout << "-----------Метод Жордана закончился успешно------------------" << endl;

}

double Norm(double *M, double *I, int n) {
///	if(stop_threads) {
///	        return 0;
///	}
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
