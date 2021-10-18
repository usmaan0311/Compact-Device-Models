#include<iostream>
#include"GildenBalt.h"


int main(int argc, char* argv[])
{

double Vg=0.5, tox=1e-7, Na=1e+16, phin=1.0;
Gildenbalt mos1(Vg, tox, Na, phin);
double phis=mos1.Phis();
mos1.ShowPar();
std::cout<<"Phis =\t"<<phis<<std::endl;


return 0;
}

