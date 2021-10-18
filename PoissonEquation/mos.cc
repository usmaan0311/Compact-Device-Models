#include<iostream>
#include"Poisson.h"

int main(int argc, char* argv[])
{
int iter, maxiter=2000000, nx=atoi(argv[1]);
double init=1.0, tsi=200e-7, tox=2e-7, Na=1e+17, Vgs=1.5, tol=atof(argv[2]);
Poisson pos(tsi=tsi, tox=tox, Na=Na, Vgs=Vgs, tol=tol, maxiter=maxiter, nx=nx );

double* X=new double [nx];
double* A=new double [nx];
double* R=new double [maxiter];
int* Itr=new int [maxiter];
R=pos.RVec(R);
Itr=pos.IVec(Itr);

X=pos.XVec(X);
A=pos.Init(A,init);
pos.PrintVal();
A=pos.FinDiff(A,R,Itr,&iter);

std::cout<<"converged in Iterations:\t"<<iter<<std::endl;

pos.WriteSol(A,X,R,Itr,"solution.dat","Residual.dat");

delete [] X;
delete [] A;
delete [] R;
delete [] Itr;

return 0;

}
