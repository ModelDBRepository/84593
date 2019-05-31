#include <stdio.h>
#include <math.h>
#include "header.h"
#include "myutils.c"
#include "rk4.c"
#include "nrutil.c"

/***********************************
Variables:
1  Vs 
   m is replaced by m_inf
2  h
3  n
4  Ca_s
5  Vd
6  Ca_d
   m_Ca is replaced by m_Ca_inf
7  Na_s
************************************/

#define NUM_var 7  

void derivatives(double, double [], double []);
void copy_step(int, double [], double []);

/*global variables */
double I_soma;

int main(void)
{
  FILE *outfp,*f_inst,*cycles;
  double t,y[1+NUM_var],y_next[1+NUM_var],dydx[1+NUM_var]; /*y[0] not used*/
  extern double I_soma;
  int c=0,first_passing;
  double spk_old,spk_new;
  int cycle_spk;
  
  if((outfp=fopen("out.dat.a3","w"))==NULL
     || (f_inst=fopen("f_inst.dat.a3","w"))==NULL
     || (cycles=fopen("cycles.dat.a3","w"))==NULL){
    printf("can't open files\n");
    exit(0);
  }
  
  t=0.0;
  I_soma=0.0;
  first_passing=1;
  spk_old=-1.0;
  cycle_spk=0;



  y[1]=V_soma;
  y[2]=h_inf(V_soma);
  y[3]=n_inf(V_soma);
  y[4]=0.0;        /* initial Soma Ca concentration */
  y[5]=V_den;
  y[6]=0.0;        /* initial Dendrite Ca concentration */
  y[7]=Na_eq;

  dydx[1]=(-g_L*(y[1]-V_soma)       /* I_L */
	   -g_Na*power_3(m_inf(y[1]))*y[2]*(y[1]-V_Na) /* I_Na*/
	   -g_K*power_4(y[3])*(y[1]-V_K)    /* I_K */
	   -g_Ca_soma*power_2(m_Ca_inf(y[1]))*(y[1]-V_Ca)   /* I_Ca */
	   -g_AHP_soma*(y[4]/(y[4]+KD))*(y[1]-V_K)         /* I_AHP */
	   -(gc/p)*(y[1]-y[5]) 
	   -g_AHP_KNa*m_KNa_inf(y[7])*(y[1]-V_K)           /* I_KNa */
	   +I_soma
	   )/Cm;

  dydx[2]=0.0;
  dydx[3]=0.0;
  dydx[4]=-influx_Ca_soma*(g_Ca_soma*power_2(m_Ca_inf(y[1]))*(y[1]-V_Ca))
    -y[4]/tau_Ca_soma;
  
  dydx[5]=(-g_L*(y[5]-V_den)                               /* I_L */
	   -g_Ca_den*power_2(m_Ca_inf(y[5]))*(y[5]-V_Ca)   /* I_Ca */
	   -g_AHP_den*(y[6]/(y[6]+KD))*(y[5]-V_K)         /* I_AHP */
	   -(gc/(1-p))*(y[5]-y[1])
	   )/Cm;

  dydx[6]=-influx_Ca_den*(g_Ca_den*power_2(m_Ca_inf(y[5]))*(y[5]-V_Ca))
    -y[6]/tau_Ca_den;

  dydx[7]=(-influx_Na_soma*(g_Na*power_3(m_inf(y[1]))*y[2]*(y[1]-V_Na))      
	   -3.0*A*R_pump*(phi_Na(y[7])-phi_Na(Na_eq)))*phi_factor;

  /* start simulation */
  for(t=0.0; t<RUNTIME; t+=dt){
    if(t>=20000.0 && t<=40000.0) 
      I_soma=I_base+hi_fluc_scale*sin(-0.5*pi + 2.0*pi*2.0*(t/1000.0));
    else I_soma=I_base+low_fluc_scale*sin(-0.5*pi + 2.0*pi*2.0*(t/1000.0));
    
    derivatives(t,y,dydx); 

    rk4(y,dydx,NUM_var,t,dt,y_next,derivatives);
    
    if(within_proxy(t,0.0,dt)){
      printf("%f: %d\n",t,cycle_spk);
      fprintf(cycles,"%f %d\n",t,cycle_spk);
      cycle_spk=0;
    }
    
    /* compose the instantaneous firing rate */
    if(first_passing && y_next[1]>=-30.0){
      cycle_spk++; /*printf("spiking\n");*/
      if(spk_old<0.0)	spk_old=t;      
      else{
	spk_new=t;
	fprintf(f_inst,"%f %f\n",spk_old,1000.0/(spk_new-spk_old));
	fprintf(f_inst,"%f %f\n",spk_new,1000.0/(spk_new-spk_old));
	spk_old=spk_new;
      }
      first_passing=0;
    }else if(y_next[1]<-30.0){
      first_passing=1;
    }

    
    if((c%selectprint)==0){
      fprintf(outfp,"%f %f %f %f %f %f %f %f %f\n",t,y_next[1],
	      y_next[2],y_next[3],y_next[4],
	      y_next[5],y_next[6],y_next[7],I_soma);
      c=1;
    }else c++;
    
    copy_step(NUM_var,y,y_next);
  }
  fclose(outfp);
  fclose(f_inst);
  fclose(cycles);
  return;
}

void derivatives(double x, double y[], double dydx[])
{
  dydx[1]=(-g_L*(y[1]-V_soma)       /* I_L */
	   -g_Na*power_3(m_inf(y[1]))*y[2]*(y[1]-V_Na) /* I_Na*/
	   -g_K*power_4(y[3])*(y[1]-V_K)    /* I_K */
	   -g_Ca_soma*power_2(m_Ca_inf(y[1]))*(y[1]-V_Ca)   /* I_Ca */
	   -g_AHP_soma*(y[4]/(y[4]+KD))*(y[1]-V_K)         /* I_AHP */
	   -(gc/p)*(y[1]-y[5])                        /* coupling */
	   -g_AHP_KNa*m_KNa_inf(y[7])*(y[1]-V_K)
	   +I_soma
	   )/Cm;

  dydx[2]=phi*(alpha_h(y[1])*(1-y[2])-beta_h(y[1])*y[2]);
  dydx[3]=phi*(alpha_n(y[1])*(1-y[3])-beta_n(y[1])*y[3]);
  dydx[4]=-influx_Ca_soma*(g_Ca_soma*power_2(m_Ca_inf(y[1]))*(y[1]-V_Ca))
    -y[4]/tau_Ca_soma;

  dydx[5]=(-g_L*(y[5]-V_den)       /* I_L */
	   -g_Ca_den*power_2(m_Ca_inf(y[5]))*(y[5]-V_Ca)   /* I_Ca */
	   -g_AHP_den*(y[6]/(y[6]+KD))*(y[5]-V_K)         /* I_AHP */
	   -(gc/(1-p))*(y[5]-y[1])                     /* coupling */
	   )/Cm;

  dydx[6]=-influx_Ca_den*(g_Ca_den*power_2(m_Ca_inf(y[5]))*(y[5]-V_Ca))
    -y[6]/tau_Ca_den;

  dydx[7]=(-influx_Na_soma*(g_Na*power_3(m_inf(y[1]))*y[2]*(y[1]-V_Na))   
	   -3.0*A*R_pump*(phi_Na(y[7])-phi_Na(Na_eq)))*phi_factor;
  return;
}


void copy_step(int n, double y[], double y_next[])
{
  int i;
  for(i=0; i<n; i++)
    y[1+i]=y_next[1+i];
  return;
}
