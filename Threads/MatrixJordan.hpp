#pragma once
#include <mutex>
#include <condition_variable>
#include <atomic>

using namespace std;

void RowsAddition(const int &size_, const int &i, const int &j, const double &h, double *matrix_, double *inverse_);
void RowMultiply(const int &size_, const int &i, const double &h, double *matrix_, double *inverse_);
void RowsSwap(const int &size_, const int &i, const int &j, double *matrix_, double *inverse_);
int NonZero(const int &size_, const int &j, double *matrix_);
void Print(const int &n, const int &l, const int &m, double *matrix_);
double Norm(double *M, double *I, int n);

class barrier {
        public:
        const unsigned int threadCount;
        std::atomic<unsigned int>threadsWaiting;
//      std::atomic<unsigned int>threadsOut;
        std::condition_variable waitVariable;
//      std::condition_variable waitOut;
        std::mutex mutex;
        public:
        barrier(unsigned int n) : threadCount(n) {
                threadsWaiting = 0;
//              threadsOut = 0;
        }
        barrier(const barrier &) = delete;
        void wait() {
                std::unique_lock<std::mutex> lock(mutex);
                if (threadsWaiting.fetch_add(1) >= threadCount - 1) {
                        threadsWaiting.store(0);
                        waitVariable.notify_all();
                }
                else {
                        waitVariable.wait(lock);
                }
/*
                if(threadsOut.fetch_add(1) >= threadCount - 1) {
                        threadsOut.store(0);
                        waitOut.notify_all();
                }
                else {
                        waitOut.wait(lock);
                }    */
        }
};

void JordanInverse(const int n, const int m,  double *M, double *I, const int start, const int stop, barrier *my_barrier);


