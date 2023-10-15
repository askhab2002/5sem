#include <iostream>
#include <math.h>
#include <limits>

using namespace std;

double Multiplier(const int &n, double *M) {
	double max1 = 0;
	double max2 = 0;

	for(int i = 0; i < n; i++) {
		if(fabs(M[i * n + i]) > max1) {
			max1 = fabs(M[i * n + i]);
		}
	}

	for(int i = 0; i < n - 1; i++) {
                if(fabs(M[(i + 1) * n + i]) > max2) {
                        max2 = fabs(M[(i + 1) * n + i]);
                }
        }
//        cout << "   max1 = " << max1 << "    max2 = " << max2 << endl;

	if(max2 > max1) {
		max1 = max2;
	}

	return 1/(4 * max1);
}


int Sign(const int &n, double *M) {
	double Alpha = Multiplier(n, M);
//	cout << "   Alpha = " << std::uppercase << std::scientific << Alpha << endl;

	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			M[i * n + j] *= Alpha;
		}
	}

	double x = M[0];
	double y = 1;
	int m = 0;

	if(x * y < 0) {
		m++;
	}

	double a, b, c, max;
	double u, v;

	for(int i = 1; i < n; i++) {
		a = M[i * n + i];
		b = M[i * n + i - 1];

		if(fabs(b) < std::numeric_limits<double>::epsilon()) {
			if(b < 0) {
				b = -1;
			}
			else {
				b = 1;
			}
		}

		if(fabs(x) < fabs(b * b * y)) {
			max = fabs(b * b * y);
		}
		else { max = fabs(x); }

		c = (1/std::numeric_limits<double>::epsilon())/max;

		u = c * (a * x - b * b * y);
		v = c * x;
	        
//	        cout << "  u = " << u << "  x = " << x << endl;	
		if(u * x < 0) {
			m++;
		}

		x = u;
		y = v;
	}

        Alpha = 1/Alpha;
	for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                        M[i * n + j] *= Alpha;
                }
        }  

	return m;
}

int SignChange(const int &n, double *M, const double &k) {
	for(int i = 0; i < n; i++) {
		M[i * n + i] = M[i * n + i] - k;
	}

	int m = Sign(n, M);

	for(int i = 0; i < n; i++) {
		M[i * n + i] = M[i * n + i] + k;
	}

	return m;
}



