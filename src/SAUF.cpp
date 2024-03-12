#include "SAUF.h"

using namespace std;
using namespace arma;


int NrThreads;

extern const double NaN=numeric_limits<double>::quiet_NaN();
extern const double Infinity=numeric_limits<double>::infinity();



void ParallelFor(const bint _beg, const bint _end, function<void(bint)> parfor, const bint NrThread) {
	for (bint i = _beg; i < _end; i += NrThread) {
		vector<thread> threads(NrThread);
		for (bint j = i; j < i + NrThread && j < _end; j++) {
			threads[j % NrThread] = thread(parfor, j);
		}
		for (bint j = i; j < i + NrThread && j < _end; j++) {
			threads[j % NrThread].join();
		}
	}
}

void ParallelFor(const bint _end, function<void(bint)> parfor, const bint NrThread) {
	ParallelFor(0, _end, parfor, NrThread);
}



