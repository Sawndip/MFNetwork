#include <math.h>
#include <stdio.h>

#define EPSILLON 0.000001
#define SQPI sqrt(4.*atan(1.))  /*square root of PI*/
#define TWOPI 8.*atan(1.)  /* two times PI */
#define a1  -1.26551223
#define a2  1.00002368
#define a3  .37409196
#define a4  .09678418
#define a5  -.18628806
#define a6  .27886087
#define a7  -1.13520398
#define a8 1.48851587
#define a9 -.82215223
#define a10 .17087277

/***** Time constants ****/
//#define taui 20.e-3       /*** membrane time constant in seconds ***/ /* now passed in as parameter */
#define taurp 0. //20.e-3 //2.e-3 //5.e-3       /*** absolute refractory period in seconds ***/

float nerf(float z)           /* function exp(z^2)(erf(z)+1)  */
{
        float t,ef,at;
        float w;
	//float bt,fex;
        w = fabs(z);
        t = 1.0e0/(1.0e0 + 0.5e0 * w);
        at=a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10))))))));
        ef=t*exp(at);
        if(z>0.0e0)
                ef = 2.0e0*exp(w*w)-ef;  
        return(ef);
}


double trans(double x,double y,double taui)     /* transduction function 
x=(threshold-mu)/sigma, y=(reset-mu)/sigma */ 
{       
  double w,z,cont,ylow;
  int i,N=10000;
	//int n;
  w=0.0e0;
  if(x<-100.&&y<-100.) {
    w=log(y/x)-0.25/pow(x,2.)+0.25/pow(y,2.);
    w=1./(taurp+taui*w);
  }
  else if(y<-100.) {
    ylow=-100.; 
    N=(int)(100.*(x-ylow) + EPSILLON); /* add epsillon to avoid numerical errors in cast */
    for(i=0;i<=N;i++) {
      z=ylow+(x-ylow)*(double)(i)/(double)(N);
      cont=nerf(z);
      if(i==0||i==N) w+=0.5*cont;
      else w+=cont;
    }
    w*=(x-ylow)*SQPI/(double)(N);
    w+=log(-y/100.)-0.000025+0.25/pow(y,2.);
    w=1.0e0/(taurp+taui*w);
  }
  else {
    ylow=y;
    N=(int)(100.*(x-ylow));
    for(i=0;i<=N;i++) {
      z=ylow+(x-ylow)*(double)(i)/(double)(N);
      cont=nerf(z);
      if(i==0||i==N) w+=0.5*cont;
      else w+=cont;
    }
    w*=(x-ylow)*SQPI/(double)(N);
    w=1.0e0/(taurp+taui*w);
  }
  return(w);
}

double function(double x) /* called from within cv() */
{
  float w;
  float y,ymin=-20.;
  int i,N=10000;
  w=0.0e0;
  for(i=0;i<=N;i++) {       
    y=ymin+(x-ymin)*(double)(i)/(double)(N);
    if(i==0||i==N) w+=0.5*pow(nerf(y)*exp(0.5*(x*x-y*y)),2.);
    else w+=pow(nerf(y)*exp(0.5*(x*x-y*y)),2.);
  }
  w*=(x-ymin)/(double)(N);
  return(w);
}

double cv(double x,double y,double taui)     /* coefficient of variation  
x=(threshold-mu)/sigma, y=(reset-mu)/sigma */ 
{       
  double w,z,cont,ylow;
  int i,N=10000;
	//double zmin,v;
	//int n,j;
  w=0.0e0;
  if(x<-100.&&y<-100.) {
    w=0.5/pow(x,2.)-0.5/pow(y,2.);
  }
  else if(y<-100.) {
    ylow=-100; 
    N=(int)(100.*(x-ylow));
    for(i=0;i<=N;i++) {
      z=ylow+(x-ylow)*(double)(i)/(double)(N);
      if(z>-10.) cont=function(z);
      else cont=1./(-TWOPI*pow(z,3.))+1.5/(TWOPI*pow(z,5.));
      if(i==0||i==N) w+=0.5*cont;
      else w+=cont;
    }
    w*=(x-ylow)*TWOPI/(double)(N);
    w+=0.00005-0.5/pow(y,2.);
  }
  else {
    ylow=y;
    N=(int)(100.*(x-ylow));
    for(i=0;i<=N;i++) {
      z=ylow+(x-ylow)*(double)(i)/(double)(N);
      if(z>-10.) cont=function(z);
      else cont=1./(-TWOPI*pow(z,3.))+1.5/(TWOPI*pow(z,5.));
      if(i==0||i==N) w+=0.5*cont;
      else w+=cont;
    }
    w*=(x-ylow)*TWOPI/(double)(N);
  }
  w*=pow(trans(x,y,taui)*taui,2.);
  return(sqrt(w));
}



//main() { 
//  double theta,hvr,desiredrate,rate,oldrate,coefvar;
//  double mu,si,muext;
//  float x,y;
//  int iter;
//
//  theta=20.;
//  hvr=17.;
//  si=5.;
//  for(mu=10.;mu<20.1;mu+=2.5) {
//    rate=trans((theta-mu)/si,(hvr-mu)/si);
//    coefvar=cv((theta-mu)/si,(hvr-mu)/si);
//    printf("mu: %f, rate: %f, cv: %f\n",mu,rate,coefvar);
//  }
//}
