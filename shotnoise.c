#include <math.h>
#include <stdio.h>

#include "shotnoise.h"

#define tauca 0.0227
//#define cpre 0.56
//#define cpost 1.24
#define thetad 1
#define thetap 1.3
#define gammad 332.
#define gammap 725.
#define sigma 3.35
#define tau 346.
#define rhostar 0.5
#define D 0.0046
#define fup 0.5
#define cmich 5.41
#define cmin ((cpre<cpost)?cpre:cpost)
#define cmax 10.
#define dc 0.01
#define Nintc (int)(cmax/dc)
#define icpre (int)(cpre/dc)
#define icpost (int)(cpost/dc)
#define icmin (int)(cmin/dc)

float Pfirst(float c, float r) {
  return(pow(c,2.*r*tauca-1.));
}

float updateWeight(float rho_old, float stepsize, float rate, float c_pre, float c_post){
	float rho_new, d_rho, alphas[2];
	d_rho = 0;
	getAlphas(rate, c_pre, c_post, alphas);
	
	d_rho = ((gammap * alphas[1]) - rho_old * (gammap * alphas[1] + gammad * alphas[0])) / tau;
	rho_new = rho_old + (stepsize * d_rho);
	printf("rho_old: %f, d_rho: %f, rho_new: %f\n", rho_old, d_rho, rho_new);
	return rho_new;
}

void getAlphas(float rate, float c_pre, float c_post, float *alphas){
  float P[Nintc],intP,alphad,alphap,pcmcpre,pcmcpost,r;
  float Gammap,Gammad,c;
  float rhobar,sigmap,taueff,UP,DOWN,synchange;
  int ic;	
	float cpre = c_pre;
	float cpost = c_post;
  //for(r=1.;r<50.5;r+=1.) {
	r = rate;
    P[(int)(cmin/dc)] = Pfirst(cmin,r);
    intP = pow(cmin,2.*r*tauca)/(2.*r*tauca);
    if(thetad<cmin) 
		alphad = (pow(cmin,2.*r*tauca)-pow(thetad,2.*r*tauca))/(2.*r*tauca);
    else 
		alphad = 0.;
	
    if(thetap<cmin) 
		alphap = (pow(cmin,2.*r*tauca)-pow(thetap,2.*r*tauca))/(2.*r*tauca);
    else 
		alphap = 0.;
	
    /*printf("%f\n",cmin);*/
    for(ic=1;ic<icmin+1;ic++) {
      c = ic*dc;
      P[ic] = Pfirst(c,r);
    }
	
    for(ic=icmin+1;ic<Nintc;ic++) {
      c = ic*dc;
      if(ic<icpre+1) 
		  pcmcpre=0.;
      else if(ic==icpre+1) 
		  pcmcpre=0.5*pow(dc/cpre,2*r*tauca)/(r*tauca);
      else 
		  pcmcpre=0.5*dc*(P[(int)((c-cpre)/dc)]/pow(c,2.*r*tauca)+P[(int)((c-cpre)/dc)-1]/pow(c-dc,2.*r*tauca));
      
	  if(ic<icpost+1) 
		  pcmcpost=0.;
      else if(ic==icpost+1) 
		  pcmcpost=0.5*pow(dc/cpost,2*r*tauca)/(r*tauca);
      else 
		  pcmcpost=0.5*dc*(P[(int)((c-cpost)/dc)]/pow(c,2.*r*tauca)+P[(int)((c-cpost)/dc)-1]/pow(c-dc,2.*r*tauca));
      
	  P[ic]=pow(c,2*r*tauca-1)*(P[ic-1]/pow(c-dc,2*r*tauca-1)-r*tauca*(pcmcpre+pcmcpost));
      intP+=0.5*dc*(P[ic-1]+P[ic]);
      /*printf("%f %f %f %f %f\n",c,pcmcpre,pcmcpost,P[ic],intP);*/
      if(c>thetad) 
		  alphad += 0.5*dc*(P[ic-1]+P[ic]);
      if(c>thetap) 
		  alphap += 0.5*dc*(P[ic-1]+P[ic]);
    }
	
    alphad /= intP;
    alphap /= intP;
    Gammad = gammad*alphad;
    Gammap = gammap*alphap;
    /*for(ic=1;ic<Nintc;ic++) {
      c=ic*dc;
      printf("%f %f\n",c,P[ic]/intP);
    }
    printf("\n");*/
    /*printf("%f %f %f %f\n",r,alphad,alphap,Gammap/(Gammad+Gammap));*/
    rhobar = Gammap/(Gammad+Gammap);
    sigmap = sigma*sqrt((alphap+alphad)/(Gammap+Gammad));
    taueff = tau/(Gammap+Gammad);
    UP = 0.5*(1.-erf((rhostar-rhobar+rhobar*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
    DOWN = 0.5*(1.+erf((rhostar-rhobar+(rhobar-1.)*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
    synchange = ((fup*(1.-DOWN)+(1.-fup)*UP)*cmich+fup*DOWN+(1.-fup)*(1.-UP))/(fup*cmich+1.-fup);
    printf("rate: %f, alphad: %f, alphap: %f, rhobar: %f, %f %f %f\n",r, alphad, alphap, rhobar, UP, DOWN, synchange);
  //}
	if (alphad < 0)
		alphad = 0;
	if (alphap < 0)
		alphap = 0;
	alphas[0] = alphad;
	alphas[1] = alphap;
	printf("alphas[0]: %f, alphas[1]: %f\n", alphas[0], alphas[1]);
}
    
//int main(void){
//	float rho = 0.5;
//	float stepsize = 0.01;
//	float rate = 12;
//	float c_pre = 1.24;
//	float c_post = 0.56;
//	float alphas[2];
//	
//	for(int i = 1; i < 5000; i++){
//		rho = updateWeight(rho, stepsize, rate, c_pre, c_post);
//		printf("i: %d, rho: %f\n", i, rho);
//		
//		//rate = (float)5.;
//		//getAlphas(rate, c_pre, c_post, alphas);
//	}
//	
//	return 0;
//}
