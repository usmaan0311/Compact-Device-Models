#include<iostream>
#include<fstream>
#include<cassert>
#include"GildenBalt.h"

double* allot(double*,double,int); 
void write(double*,double*,int);

int main(int argc, char* argv[])
{

int n=atoi(argv[1]);
double tox=1e-7, Na=1e+16, phin=0.0, Vgm=atof(argv[2]);
double* Vg=new double [n];
double* Ps=new double [n];
allot(Vg,Vgm,n);
for(int i=0;i<n;i++)
{
Gildenbalt mos1(Vg[i], tox, Na, phin);
Ps[i]=mos1.Phis();
mos1.ShowPar();
std::cout<<"Phis =\t"<<Ps[i]<<std::endl;
}

write(Ps,Vg,n);

delete [] Ps;
delete [] Vg;

return 0;
}

double* allot(double* A, double Vgm, int n)
{
double dv = Vgm/(n -1);
for(int i=0; i<n; i++)
{

	A[i] = 0.25 + (i+1)*dv;

}

return A;
}

void write(double* Ps, double* Vg, int n)
{
	std::ofstream writePs("Phis.dat");
	writePs.precision(10);
	writePs.setf(std::ios::scientific);
	assert(writePs.is_open());
	
	for(int i=0; i<n; i++)
	{

		writePs<<Vg[i]<<" "<<Ps[i]<<"\n";

	}
	writePs.close();


}
