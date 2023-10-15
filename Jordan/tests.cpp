#include <iostream>
#include <math.h>
#include <string>
#include "MatrixJordan.hpp"
#include "MatrixRead.hpp"

using namespace std;

int main(int argc, char *argv[]) { 
	if(argc != 4 && argc != 5) {
		cout << "    Вы ввели некорректную командную строку" << endl;
		return 0;
	} 
        cout << " ------------------------------ " << endl;

	string n_(argv[1]);
	string m_(argv[2]);
	string k_(argv[3]);
	int n = stoi(n_);
	int m = stoi(m_);
	int k = stoi(k_);

	clock_t start, stop;

	string Filename = "";
	
	if(argc == 5) {
	        Filename = string(argv[4]);
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
	JordanInverse(n, m, M, Inverse);
	stop = clock();

	cout << "    Затраченное время на нахождение обратной матрицы - " << std::uppercase << std::scientific << ((double)stop - (double)start)/((double)CLOCKS_PER_SEC) << endl;

	delete [] M;

        M = MatrixRead(n, m, k, Filename, &correctness);
	Norm(M, Inverse, n);
	
	delete [] Inverse;
	delete [] M;
	return 0;
}
