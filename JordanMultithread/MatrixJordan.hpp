
using namespace std;

void RowsAddition(const int &size_, const int &i, const int &j, const double &h, double *matrix_, double *inverse_);
void RowMultiply(const int &size_, const int &i, const double &h, double *matrix_, double *inverse_);
void RowsSwap(const int &size_, const int &i, const int &j, double *matrix_, double *inverse_);
int NonZero(const int &size_, const int &j, double *matrix_);
void Print(const int &n, const int &l, const int &m, double *matrix_);
void JordanInverse(const int &n, const int &m, double *M, double *I, const int &threads);
double Norm(double *M, double *I, const int &n);
