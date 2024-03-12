#include "RBF.h"
using namespace arma;


const unsigned int n = 2;
int main(int argc, char **argv) {
	bool packed = true;
	function<vec(const vec &)> f = [](const vec &x)->vec {
		vec fx(n);
		fx(0) = x(1);
		fx(1) = -x(0) + 1.0 / 3.0 * pow(x(0), 3) - x(1);
		return fx;
	};
	function<double(const vec &)> pq = [&f](const vec &x)->double {
		return pow(norm(x, 2), 2) * (1 + pow(norm(f(x), 2), 2));
	};
	function<double(double,double)> psi1 = [](double r, double c)->double {
		return 1 -c * r > 0 ? -22 * pow(c, 2) * pow(1 - c * r, 7) * (1 + 7 * c * r + 16 * pow(c * r, 2)) : 0;
	};
	function<double(double,double)> psi2 = [](double r, double c)->double {
		return 1 - c * r > 0 ? 528 * pow(c, 4) * pow(1 - c * r, 6) * (1 + 6 * c * r) : 0;
	};

	/*
	function<double(double,double)> psi1 = [](double r, double c)->double {
		return 1 - c * r > 0 ? -130 * pow(c, 2) * pow(1 - c * r, 9) * (5 + 45 * c * r + 159 * pow(c * r, 2) + 231 * pow(c * r, 3)) : 0;
	};
	function<double(double,double)> psi2 = [](double r, double c)->double {
		return 1 - c * r > 0 ? 17160 * pow(c, 4) * pow(1 - c * r, 8) * (1 + 8 * c * r + 21 * pow(c * r, 2)) : 0;
	};
	*/
	double c = 1.5;
	RBFLya R(f,pq,psi1,psi2,c);
	bint xN=200,yN=250;  // must be even numbers to avoid (0,0)
	//bint xN=300,yN=250;
	//bint xN=250,yN=400;
	//bint xN=250,yN=500;
	//bint xN=500,yN=500;
	double xMin=-2.6, yMin=-2.6;
	double xMax=2.6, yMax=2.6;
	vector<vec> X(xN*yN);
	for(bint i=0;i<xN;i++){
		for(bint j=0;j<yN;j++){
		    X[j+yN*i] = vec {xMin+i*(xMax-xMin)/(xN-1), yMin+j*(yMax-yMin)/(yN-1)};		
		}
	}
	R.FixVertices(X);
	wall_clock timer;
	if(packed == true){
		timer.tic();
		R.WriteA();
		cout<<"write A packed "<<timer.toc()<<endl;
		timer.tic();
		R.SolveRBF();
		cout<<"solve packed "<<timer.toc()<<endl;
		R.A.set_size(0);
	}
	else{ 
		timer.tic();
		R.WriteAm();
		cout<<"write A unpacked "<<timer.toc()<<endl;
		timer.tic();
		R.SolveRBFm();
	}
	cout<<"solved in "<<timer.toc()<<" sec"<<endl;


	vec xMatlab(xN), yMatlab(yN);
	for(bint i=0;i<xN;i++){
		xMatlab(i)=X[i*yN](0);
	}
	for(bint j=0;j<yN;j++){
		yMatlab(j)=X[j](1);
	}
	mat VMatlab(yN,xN), OrbDerVMatlab(yN,xN);
	ParallelFor(yN,[&](bint j)->void {
		for(bint i=0;i<xN;i++){
			VMatlab(j,i)=R.V(X[j+yN*i]);
			OrbDerVMatlab(j,i)=R.OrbDerV(X[j+yN*i]);
		}
	});
	xMatlab.save("x.txt", raw_ascii);
	yMatlab.save("y.txt", raw_ascii);
	VMatlab.save("V.txt", raw_ascii);
	OrbDerVMatlab.save("OrbDerV.txt", raw_ascii);
}

/* to plot with Matlab / Octave
	// plot V
	load -ascii 'x.txt'
	load -ascii 'y.txt'
	load -ascii 'V.txt'
	[X,Y]=meshgrid(x,y);
	surf(X,Y,V)
	xlabel('X')
	ylabel('Y')
	zlabel('V(X,Y)')

	// plot V'
 	load -ascii 'x.txt'
 	load -ascii 'y.txt'
	load -ascii 'OrbDerV.txt'
	[X,Y]=meshgrid(x,y);
	surf(X,Y,OrbDerV)
	xlabel('X')
	ylabel('Y')
	zlabel("V'(X,Y)")
*/
