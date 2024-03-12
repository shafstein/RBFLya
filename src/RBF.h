#ifndef _RBF_H_
#define _RBF_H_
#include "SAUF.h"

 class RBFLya {
	public:
    bint N;
    double c;
    function<vec(vec)> f;
	function<double(vec)> pq;
	function<double(double,double)> psi1, psi2;
    vec alpha, beta;
    vector<vec> X, fX;
	vec A;
	mat Am;
   public:
    RBFLya(function<vec(vec)> _f, function<double(vec)> _pq, function<double(double, double)> _psi1, function<double(double, double)> _psi2, double _c):
        f(_f), pq(_pq), psi1(_psi1), psi2(_psi2), c(_c), N(0) {  
	};
	void FixVertices(const vector<vec> &_X); 
	void WriteA(void);
	void WriteAm(void);
	void SolveRBF(void); 
	void SolveRBFm(void); 
	double V(const vec &x) const;
	double OrbDerV(const vec &x) const; 
 };

#endif