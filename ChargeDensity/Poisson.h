#ifndef POISSONHEADERDEF
#define POISSONHEADERDEF
#include<string>

class Poisson{
	protected:
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
	double phif;
	double tdep;
	double tol;
	double dx;
	int maxiter;
	int nx;

	public:
		Poisson(double,double,double,double, int, int);
		void PrintVal();
		double* XVec(double*);
		double* RVec(double*);
		int* IVec(int*);
		double* Init(double*,double);
		double l2_norm(double*, double*);
		double* FinDiff(double*,double*,int*, int*,double);
		void WriteSol(double*, double*,double*,int*,std::string, std::string);


};

class Charge : public Poisson
{
	double n0 = ni*ni/Na;
	public:
		Charge(double,double,double,double,int,int);
		double* N_inv(double*, int);
		double Simpson(double*,int, double);
		void WriteCharge(double*,double*,double*,double*,std::string,std::string,int);


};


#endif
