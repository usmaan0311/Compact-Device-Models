#ifndef POISSONHEADERDEF
#define POISSONHEADERDEF
#include<string>

class Poisson{

const double q=1.609e-19;
const double ep0=8.4e-14;
const double vt=25.9e-3;
double eps_s=11.9*ep0;
double eps_ox=3.9*ep0;
double ni=1e+10;
double Na;
double tox;
double tsi;
double Cox;
double phi_ms;
//double phif;
double Vgs;
double tdep;
double tol;
double dx;
int maxiter;
int nx;

	public:
		double phif;
		Poisson(double,double,double,double,double, double, int, int);
		void PrintVal();
		double* XVec(double*);
		double* RVec(double*);
		int* IVec(int*);
		double* Init(double*,double);
		double l2_norm(double*, double*);
		double* FinDiff(double*,double*,int*, int*);
		void WriteSol(double*, double*,double*,int*,std::string, std::string);


};




#endif
