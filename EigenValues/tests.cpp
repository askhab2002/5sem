#include <iostream>
#include <math.h>
#include <string>
#include "MatrixRead.hpp"
#include "MatrixTridiagonal.hpp"
#include "SignChange.hpp"
#include "EigenValue.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	if(argc != 6 && argc != 7) {
		cout << "    Вы ввели некорректную командную строку   " << endl;
		return 0;
	}

	string n_(argv[1]);
	string m_(argv[2]);
        string e_(argv[3]);
	string k_(argv[4]);
	string l_(argv[5]);
	int n = stoi(n_);
	int m = stoi(m_);
	int k = stoi(k_);
	int l = stoi(l_);
	double e = stod(e_);
        cout << " e = " << e << endl;

	clock_t start, stop;

	string Filename = "";

	if(argc == 7) {
		Filename = string(argv[6]);
	}

	int correctness = 0;
       
        	
	
	double *M = MatrixRead(n, m, k, Filename, &correctness); 

	if(correctness != 1) {
		cout << "    Matrix hasnt read" << endl;
		return 0;
	}

	cout << "-----------Введенная матрица:--------------" << endl;

	Print(n, n, m, M);
        
	start = clock();
	Tridiagonal(n, M);
	stop = clock();
        cout << "    Time to tridiagonolize:    " << std::uppercase << std::scientific << ((double)stop - (double)start)/((double)CLOCKS_PER_SEC) << endl;
	

	cout << "-----------Тридиагональная матрица:-----------" << endl;

	Print(n, n, m, M);
        
	start = clock();
	double Eigen = EigenValue(n, M, l, e);
	stop = clock();

	cout << " -------------Собственное значение = " << Eigen << endl;
        
	cout << "    Time to find eigenvalue:    " << std::uppercase << std::scientific << ((double)stop - (double)start)/((double)CLOCKS_PER_SEC) << endl;

	delete [] M;
	
	return 0;
}
