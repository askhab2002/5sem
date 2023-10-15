#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

using namespace std;

//Функция получает адреса переменных x, i, j и название файла, откуда будет считываться последовательность
//Возвращает 0 при успешном завершении, -1 - ошибка открытия файла, -2 - ошибка данных в файлах(нет цифр).
//x - максимальное значение последрвательности, i - номер первого максимума, j - последний максимум.

int SequenceMax(double *x, int *i, int *j, const string file) {
	ifstream input(file);
	if(!(input.is_open())) {
		return -1;
	}

	double value;
        vector<double> sequence;
        
	while(input >> value) {
                sequence.push_back(value);
	}
        
	if(sequence.size() == 0) {
		return -2;
	}

	input.close();
        
	value = sequence[0];
	*i = 1;
        
	cout << "  Элементы последовательности: . Сделал для удобства проверки, хотя оно не требуется в задаче." << endl;

        for(int k = 0; k < sequence.size(); k++) {
		if(fabs(value - sequence[k]) < 1e-5) {
			*j = k + 1;
		} 
		else { 
		        if(sequence[k] > value) {
			        value = sequence[k];
			        *i = k + 1;
				*j = *i;
			}
		}
	        cout << sequence[k] << endl;
        }	
        
	*x = value;
	return 0;
}
