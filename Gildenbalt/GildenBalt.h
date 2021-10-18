#ifndef GILDENBALTHEADERDEF
#define GILDENBALTHEADERDEF

class Gildenbalt
{
const double ep0=8.84e-14;
const double e=1.609e-19;
const double Vt=25.9e-3;
double Vfb=0.21;
double eps_s=11.4*ep0;
double eps_ox=3.9*ep0;
double ni=1e+10;
double mu=1000;
double Vgb;
double tox;
double Na;
double Cox;
double Gamma;
double phif;
double phin;
double G;
double xg;
double xn;
double xsub;
	public:
	Gildenbalt(double Vgb,double tox, double Na, double phin);

	double s(double a, double b, double c);
	double sigma(double a, double c, double tau);
	double Phis();
	void ShowPar();

};

#endif
