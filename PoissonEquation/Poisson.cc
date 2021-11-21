#include<iostream>
#include<cmath>
#include<fstream>
#include<cassert>
#include<string>
#include"Poisson.h"


Poisson::Poisson(double phi_ms, double tsi, double tox,double Na, double Vgs, double tol,int maxiter, int nx)
{

this->tsi=tsi;
this->tox=tox;
this->Na=Na;
this->Vgs=Vgs;
this->phi_ms=phi_ms;
Cox=eps_ox/this->tox;
phif=vt*log(this->Na/ni);
tdep=sqrt(2*eps_s*phif/(q*this->Na));
this->tol=tol;
this->maxiter=maxiter;
this->nx=nx;
dx=tsi/(this->nx -1);
}

void Poisson::PrintVal()
{

	std::cout<<" Na = "<<Na<<std::endl;
	std::cout<<" tox = "<<tox<<std::endl;
	std::cout<<" tsi = "<<tsi<<std::endl;
	std::cout<<" Cox = "<<Cox<<std::endl;
	std::cout<<" phif = "<<phif<<std::endl;
	std::cout<<" Vgs = "<<Vgs<<std::endl;
	std::cout<<" tdep = "<<tdep<<std::endl;
	std::cout<<" tol = "<<tol<<std::endl;
	std::cout<<" dx = "<<dx<<std::endl;
	std::cout<<" nx = "<<nx<<std::endl;
	std::cout<<" maxiter = "<<maxiter<<std::endl;


}


double* Poisson::XVec(double* X)
{

	for(int i=0; i<nx; i++)
	{

		X[i]=i*dx;


	}


return X;

}


double* Poisson::RVec(double* R)
{

	for(int i=0; i<maxiter; i++)
	{

		R[i]=0;


	}


return R;

}


int* Poisson::IVec(int* I)
{

	for(int i=0; i<maxiter; i++)
	{

		I[i]=0;


	}


return I;

}





double* Poisson::Init(double* A, double a)
{

	for(int i=0; i<nx; i++)
	{

		A[i]=a;


	}

return A;

}


double Poisson::l2_norm(double* A, double* B)
{

	double norm, Bsum=0.0, Dsum=0.0;

	for(int i=0; i<nx; i++)
	{

		Bsum+=pow(B[i],2);
		Dsum+=pow((A[i] - B[i]),2);

	}

	norm = sqrt(Dsum/Bsum);

return norm;

}



double* Poisson::FinDiff(double* A, double* resid,int* iteration, int* iter)
{
	double* mat=new double [nx];
	*iter=0;
	double diff=1.0 + tol;

	while(diff>tol && *iter<maxiter)
	{

		for(int i=0; i<nx; i++)
		{

			mat[i]=A[i];


		}

		for(int i=1; i<nx-1; i++)
		{
// implicit scheme
			
		A[i] = (mat[i+1] + mat[i-1] + (pow(dx,2)*q*Na/eps_s)*(  exp(-mat[i+1]/vt) -1 \
				-exp(-2*phif/vt)*( exp(mat[i+1]/vt) - 1 ) ) )/2 ;
// Boundary condition
			
		A[nx-1]=0;
	
		A[0] = ( A[1] + ((Vgs - phi_ms)*dx*Cox/eps_s) )/(1 + (dx*Cox/eps_s));
	

		}

		diff = l2_norm(A,mat);
		resid[*iter]=diff;
		iteration[*iter]=*iter;
		*iter+=1;
		
	}


return A;

}

void Poisson::WriteSol(double* A, double* X,double* R, int* Itr,  std::string strg, std::string strgr)
{
	std::ofstream WriteA(strg);
	WriteA.precision(16);
	WriteA.setf(std::ios::scientific);
	assert(WriteA.is_open());
	for(int i=0; i<nx; i++)
	{
		WriteA<<X[i]<<" "<<A[i]<<"\n";


	}

	WriteA.close();

	std::ofstream WriteR(strgr);
	WriteR.precision(16);
	WriteR.setf(std::ios::scientific);
	assert(WriteR.is_open());
	for(int i=0; i<nx; i++)
	{
		WriteR<<Itr[i]<<" "<<R[i]<<"\n";


	}

	WriteR.close();



}

