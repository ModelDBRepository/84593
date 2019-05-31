#include "header.h"

double alpha_n(double V)
{
  return -0.01*(V+34.0)/(exp(-0.1*(V+34.0))-1.0);
}

double beta_n(double V)
{
  return 0.125*exp(-(V+44.0)/25.0);
}

double n_inf(double V)
{
  return alpha_n(V)/(alpha_n(V)+beta_n(V));
}

double alpha_m(double V)
{
  return -0.1*(V+33.0)/(exp(-0.1*(V+33.0)) - 1.0);
}

double beta_m(double V)
{
  return 4.0*exp(-(V+58.0)/12.0);
}

double m_inf(double V)
{
  return alpha_m(V)/(alpha_m(V)+beta_m(V));
}

double alpha_h(double V)
{
  return 0.07*exp(-(V+50.0)/10.0);
}

double beta_h(double V)
{
  return 1.0/(exp(-0.1*(V+20.0))+1.0);
}

double h_inf(double V)
{
  return alpha_h(V)/(alpha_h(V)+beta_h(V));
}

double power_2(double x)
{
  return x*x;
}

double power_3(double x)
{
  return x*x*x;
}

double power_4(double x)
{
  return x*x*x*x;
}

double m_Ca_inf(double V)
{
  return 1.0/(1.0+exp(-(V+20.0)/9.0));
}

double m_KNa_inf(double x)
{
   return m_inf_max/(1.0+pow(EC_50/x,Hill_order));
}

double phi_Na(double x)
{
  return(x*x*x/(x*x*x+Kp*Kp*Kp));
}

int within_proxy(double x,double target,double tolerance)
{
  while(x>500.0) x-=500.0;
  
  if(x>=(target-tolerance) && x<=(target+tolerance))
    return 1;
  else return 0;
}
