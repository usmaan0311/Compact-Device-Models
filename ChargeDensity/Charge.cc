#include<string>
#include<cmath>
#include<cassert>
#include<fstream>
#include"Poisson.h"

Charge::Charge(double tsi, double tox, double Na, double tol, int maxiter, int nx):Poisson::Poisson(tsi, tox, Na, tol, maxiter, nx)
{
}

double* Charge::N_inv(double* A,int n)
{
double* qi=new double[n];
for(int i=0; i<n; i++)
{

qi[i]=0;

}

for(int i=0; i<n; i++)
{
qi[i]=q*n0*exp(A[i]/vt);

}
return qi;

}

double Charge::Simpson(double* A, int n,double delx)
{
double s=0,h=delx;
double* Nin=N_inv(A,n);
for(int i=0; i<n; i++)
{
if(i==0)
{

s+=(3*h/8)*Nin[i];
}
else if(i==n-1)
{

s+=(3*h/8)*Nin[i];

}
else
{
if(i%3 !=0)
{
s+= (9*h/8)*Nin[i];

}
else
{

s+= (6*h/8)*Nin[i];

}

}

}

return s;
}

/*

double Charge::Trapz(double* A, int n,double delx)
{
double s=0,h=delx;
double* Nin=N_inv(A,n);
for(int i=0; i<n; i++)
{
if(i==0)
{

s+=(h/2)*Nin[i];
}
else if(i==n-1)
{

s+=(h/2)*Nin[i];

}
else
{
s+=h*Nin[i];

}



}

return s;
}

*/

void Charge::WriteCharge(double* Ninv,double* X,double* Qin, double* Vgs,std::string strn,std::string strq, int ng)
{

	std::ofstream WriteN(strn);
	WriteN.precision(16);
	WriteN.setf(std::ios::scientific);
	assert(WriteN.is_open());
	for(int i=0; i<nx;i++)
	{

		WriteN<<X[i]<<" "<<Ninv[i]<<"\n";


	}
	WriteN.close();

	std::ofstream Writeq(strq);
	Writeq.precision(16);
	Writeq.setf(std::ios::scientific);
	assert(Writeq.is_open());
	for(int i=0; i<ng; i++)
	{

		Writeq<<Vgs[i]<<" "<<Qin[i]<<"\n";


	}
	Writeq.close();




}

