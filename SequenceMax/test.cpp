#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "Sequence.hpp"

using namespace std;

int main() {
	double x = -1;
	int i = -1;
	int j = -1;
        int error = -1;

	error = SequenceMax(&x, &i, &j, "text.txt");
        
	cout << " x = " << x << " i = " << i << " j = " << j << " error code = " << error << endl;
	return 0;
}

