#include <iostream>
#include <math.h>
#include <string>
#include "MatrixJordan.hpp"
#include "MatrixRead.hpp"

using namespace std;

int main(int argc, char *argv[]) { 
	if(argc != 5 && argc != 6) {
		cout << "    Вы ввели некорректную командную строку" << endl;
		return 0;
	} 
        cout << " ------------------------------ " << endl;

	string n_(argv[1]);
	string m_(argv[2]);
	string k_(argv[3]);
	string threads_(argv[4]);
	int n = stoi(n_);
	int m = stoi(m_);
	int k = stoi(k_);
	int threads = stoi(threads_);

	clock_t start, stop;

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

	cout << " -----------Введенная матрица:-------------" << endl;
        Print(n, n, m, M);

	double *Inverse = new double[n * n];
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			Inverse[i * n + j] = 0;
		}
	}
	for(int i = 0; i < n; i++) {
		Inverse[i * n + i] = 1;
	}
        
	start = clock();
	JordanInverse(n, m, M, Inverse, threads);
	stop = clock();

	cout << "    Затраченное время на нахождение обратной матрицы - " << std::uppercase << std::scientific << ((double)stop - (double)start)/((double)CLOCKS_PER_SEC) << endl;

	delete [] M;

        M = MatrixRead(n, m, k, Filename, &correctness);

	start = clock();
	Norm(M, Inverse, n);
	stop = clock();

	cout << "    Затраченное время на нахождение нормы - " << std::uppercase << std::scientific << ((double)stop - (double)start)/((double)CLOCKS_PER_SEC) << endl;
	
	delete [] Inverse;
	delete [] M;
	return 0;
}
