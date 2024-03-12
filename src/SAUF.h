#ifndef _SAUF_H_
#define _SAUF_H_
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG  //should I use this?
#define ARMA_BLAS_UNDERSCORE

// faster with, during debuggin comment out
//#define ARMA_NO_DEBUG
//#define ARMA_EXTRA_DEBUG
//#define ARMA_PRINT_ERRORS

#include<vector>
#include<functional>
#include<algorithm>	
#include<cstdlib>
#include<sstream>
#include<set>
#include<list>
#include<set>
#include<tuple>
#include<ctime>
#include<cassert>
#include<cmath>
#include<iterator>
#include<thread>
#include<random>
#include<numeric>
#include <iomanip>  


// typedef long long bint;
using bint = long long;

// undef possible stupid macros
#undef NAN
#undef min
#undef max

#include "armadillo"

using namespace std;
using namespace arma;

extern const unsigned int n;
extern const double NaN;
extern const double Infinity;

void ParallelFor(const bint _beg, const bint _end, function<void(bint)> parfor, const bint NrThread=1000);
void ParallelFor(const bint _end, function<void(bint)> parfor, const bint NrThread=1000);
#endif