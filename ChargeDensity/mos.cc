#include<iostream>
#include"Poisson.h"

int main(int argc, char* argv[])
{
int iter, maxiter=200000000, ng=50, nx=atoi(argv[1]);
double init=1.0, tsi=200e-7, tox=2e-7, Na=1e+17, Vgm=4, tol=atof(argv[2]);
Poisson pos(tsi=tsi, tox=tox, Na=Na, tol=tol, maxiter=maxiter, nx=nx );
Charge qinv(tsi=tsi, tox=tox, Na=Na, tol=tol, maxiter=maxiter, nx=nx );


double dVg=2*Vgm/(ng -1);
double* Vgs=new double [ng];
for(int i=0; i<ng;i++)
{

Vgs[i]= -Vgm + i*dVg;

}

double* Ninv=new double [nx];
double Qin[ng];

double* X=new double [nx];
double* A=new double [nx];
double* R=new double [maxiter];
int* Itr=new int [maxiter];
R=pos.RVec(R);
Itr=pos.IVec(Itr);

X=pos.XVec(X);
A=pos.Init(A,init);
pos.PrintVal();
A=pos.FinDiff(A,R,Itr,&iter,Vgs[ng-1]);

Ninv=qinv.N_inv(A,nx);

double delx=(X[nx-1] - X[0])/(nx -1);
for(int i=0; i<ng-1; i++)
{
	
A=pos.Init(A,init);
A=pos.FinDiff(A,R,Itr,&iter,Vgs[i]);
Qin[i]=qinv.Simpson(A,nx,delx);

}


std::cout<<"converged in Iterations:\t"<<iter<<std::endl;

pos.WriteSol(A,X,R,Itr,"solution.dat","Residual.dat");
qinv.WriteCharge(Ninv,X,Qin,Vgs,"InversionDensity.dat","InversionCharge.dat",ng);

delete [] Vgs;
delete [] Ninv;
delete [] X;
delete [] A;
delete [] R;
delete [] Itr;

return 0;

}
