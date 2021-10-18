#include<iostream>
#include<cmath>
#include"GildenBalt.h"

Gildenbalt::Gildenbalt(double Vgb, double tox, double Na, double phin)
{
this->Vgb=Vgb;
this->tox=tox;
this->Na=Na;
this->phin=phin;
Cox=(eps_ox)/(this->tox);
Gamma=sqrt(2*e*(eps_s)*(this->Na))/Cox;
phif=Vt*log((this->Na)/(ni));
G=Gamma/sqrt(Vt);
xg=(this->Vgb - Vfb)/Vt;
xn=(2*phif + phin)/Vt;
xsub=xg + pow(G,2)/2 - G*pow((xg - 1 + pow(G,2)/4 ),1/2);
}

double Gildenbalt::s(double a, double b, double c)
{

return 0.5*(a + b - sqrt( pow( (a - b),2 ) + c ));


}

double Gildenbalt::sigma(double a, double c, double tau)
{

double v = a + c;
double mu = pow(v,2)/tau + pow(c,2)/2 - a;
return (a*v)/(mu + (v/(mu))*c*( pow(c,2)/3 - a ) );

}

double Gildenbalt::Phis()
{
double eta = s(xsub,xn+3,5);
double a=pow( (xg - eta), 2 ) - pow(G,2)*eta + pow(G,2);
double c=2*(xg - eta) + pow(G,2);
double tau=xn - eta + log(abs(a)/pow(G,2));
double x0=eta + sigma(a,c,tau);

double Delta0 = exp(x0 - xn);
double p = 2*(xg - x0) + pow(G,2)*(1 + Delta0);
double q = pow((xg - x0),2) - pow(G,2)*(x0 + Delta0 -1);
double x = x0 + 2*q/(p + sqrt( pow(p,2) - 2*(2 - pow(G,2)*Delta0)*q ));
return x*Vt;
}
void Gildenbalt::ShowPar()
{

std::cout<<"Cox\t:\t"<<Cox<<std::endl;
std::cout<<"Gamma\t:\t"<<Gamma<<std::endl;
std::cout<<"phif\t:\t"<<phif<<std::endl;
std::cout<<"G\t:\t"<<G<<std::endl;
std::cout<<"xg\t:\t"<<xg<<std::endl;
std::cout<<"xn\t:\t"<<xn<<std::endl;
std::cout<<"xsub\t:\t"<<xsub<<std::endl;

}
