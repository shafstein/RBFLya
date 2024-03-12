#include "RBF.h"
using namespace std;
using namespace arma;


extern "C" {
	// LAPACK	
	void dpptrf_(char *UPLO, bint *N, double *A, int *INFO);
	void dpptrs_(char *UPLO, bint *N, bint *NRHS, double *AP, double *B, bint *LDA, int *INFO);
	void dpprfs_(char *UPLO, bint *NN, bint *NRHS, double *APorig, double *AP, double *Borig, bint *LDB, double *B, bint *LDX, double *FERR, double *BERR, double *WORK, bint *IWORK, int *INFO);
	// BLAS
	void dspmv_(char *UPLO, bint *N, double *ALPHA, double *AP, double *B, bint *INCX, double *BETA, double *Y, bint *INCY);
	void dtpmv_(char *UPLO, char *TRANS, char *DIAG, bint *N, double *AP, double *X, bint *INCX);
}


void RBFLya::FixVertices(const vector<vec> &_X) {
    X = _X;
    N = X.size();
	fX.resize(N);
    for (bint i=0;i<N;i++){
      fX[i] = f(X[i]);
    }
}

void RBFLya::WriteA(void) {
	A.set_size(N * (N + 1) / 2);
	beta.set_size(N);
	ParallelFor(N, [&](bint k) {
		bint offset = k * (k + 1) / 2;
		beta(k) = -pq(X[k]);
		for (bint j = 0; j <= k; j++) {
			vec x_j_m_x_k = X[j] - X[k];
			double dist = norm(x_j_m_x_k, 2);
			if (1.0 - c * dist > 0.0) {   // if not psi1(n,c)=psi2(n,c)=0 and thus A(j,k)=0
				A[offset + j] = -psi2(dist, c) * dot(x_j_m_x_k, fX[j]) * dot(x_j_m_x_k, fX[k]) - psi1(dist, c) * dot(fX[j], fX[k]);
			}
			else {
				A[offset + j] = 0.0;
			}
		}
	});
}

void RBFLya::WriteAm(void) {
	Am.set_size(N,N);
	beta.set_size(N);
	ParallelFor(N, [&](bint k) {
		beta(k) = -pq(X[k]);
		for (bint j = 0; j < N; j++) {
			vec x_j_m_x_k = X[j] - X[k];
			double dist = norm(x_j_m_x_k, 2);
			if (1.0 - c * dist > 0.0) {   // if not psi1(n,c)=psi2(n,c)=0 and thus A(j,k)=0
				Am(j,k) = -psi2(dist, c) * dot(x_j_m_x_k, fX[j]) * dot(x_j_m_x_k, fX[k]) - psi1(dist, c) * dot(fX[j], fX[k]);
			}
			else {
				Am(j,k)=0.0;
			}
		}
	});
}

void RBFLya::SolveRBF(void) {
	cout << "HERE WE GO WITH LAPACK CHOLESKY DECOMPOSITION" << endl;
	char UPLO = 'U';
	int INFO; 
    // Cholesky factorize in-place
	dpptrf_(&UPLO, &N, A.memptr(), &INFO);
	bint NRHS = 1;
	alpha = beta;
    // solve A*alpha = beta (solution overwrites alpha)
	dpptrs_(&UPLO, &N, &NRHS, A.memptr(), alpha.memptr(), &N, &INFO);
}

void RBFLya::SolveRBFm(void) {
	Am=chol(Am);
	alpha = solve(trimatl(Am.t()),beta);
	alpha=solve(trimatu(Am),alpha);
}

double RBFLya::V(const vec &x) const {
	double ret = 0.0;
	for (bint k = 0; k < N; k++) {
		vec x_k_m_x = X[k] - x;
		double dist = norm(x_k_m_x, 2);
		if (1.0 - c * dist > 0) {   // if not psi1(n,c)=0
			ret += alpha(k) * dot(x_k_m_x, fX[k]) * psi1(dist, c);
		}
	}
	return ret;
}

double RBFLya::OrbDerV(const vec &x) const {
	double ret = 0.0;
	vec fx = f(x);
	for (bint k = 0; k < N; k++) {
		vec x_k_m_x = X[k] - x;
		double dist = norm(x_k_m_x, 2);
		if (1 - c * dist > 0) {   // if not psi1(n,c)=0
			ret += -alpha(k) * (psi1(dist, c) * dot(fx, fX[k]) + psi2(dist, c) * dot(x_k_m_x, fx) * dot(x_k_m_x, fX[k]));
		}
	}
	return ret;
}