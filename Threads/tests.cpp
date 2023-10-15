#include <iostream>
#include <chrono>
#include <math.h>
#include <string>
#include <thread>
#include <vector>
#include "MatrixJordan.hpp"
#include "MatrixRead.hpp"

using namespace std;



int main(int argc, char *argv[]) { 
	if(argc != 5 && argc != 6) {
		cout << "    Вы ввели некорректную командную строку" << endl;
		return 0;
	} 
//        cout << " ------------------------------ " << endl;
//        cout << "    hardware_concurrency() threads = " << std::thread::hardware_concurrency() << endl;
	string n_(argv[1]);
	string m_(argv[2]);
	string k_(argv[3]);
	string threads_(argv[4]);
	int n = stoi(n_);
	int m = stoi(m_);
	int k = stoi(k_);
	int threads = stoi(threads_);
        
	if(threads > n) {
		threads = n;
	}
	

	string Filename = "";
	
	if(argc == 6) {
	        Filename = string(argv[5]);
	}

	int correctness = 0;
        double *M = MatrixRead(n, m, k, Filename, &correctness);

        if(correctness  != 1) {
                cout << "    Matrix hasnt read" << endl;
                return 0;
        }

///	cout << " -----------Введенная матрица:-------------" << endl;
///        Print(n, n, m, M);

	double *Inverse = new double[n * n];
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			Inverse[i * n + j] = 0;
		}
	}
	for(int i = 0; i < n; i++) {
		Inverse[i * n + i] = 1;
	}
        
        barrier *my_barrier = new barrier(threads);

        auto now = std::chrono::system_clock::now();
      
	thread *Threads = new thread[threads];
        int start = 0;
        int stop = 0;

        for(int k = 0; k < threads - 1; k++) {
                start = k * (n/threads);
                stop = (k + 1) * (n/threads);
                Threads[k] = thread(JordanInverse, n, m, M, Inverse, start, stop, my_barrier);
        }
        Threads[threads - 1] = thread(JordanInverse, n, m, M, Inverse, stop, n, my_barrier);
        for(int k = 0; k < threads; k++) {
                Threads[k].join();
        }

        delete [] Threads;

        auto then = std::chrono::system_clock::now();
	std::chrono::duration<double> ttime = then - now;

///	Print(n, n, m, Inverse);

        cout << "    Время = " << ttime.count() << endl;
	delete [] M;

        M = MatrixRead(n, m, k, Filename, &correctness);
	Norm(M, Inverse, n);
	
	delete [] Inverse;
	delete [] M;
	return 0;
}
