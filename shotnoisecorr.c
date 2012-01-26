#include <math.h>
#include <stdio.h>

#define tauca 0.0227
#define cpre 0.56
#define cpost 1.24
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
#define cmax 10.
#define dc 0.001
#define Nintc (int)(cmax/dc)
#define icpre (int)(cpre/dc)
#define icpost (int)(cpost/dc)

FILE *fopen(),*cadistrib;

float Pfirst(float c, float r) {
  return(pow(c,r*tauca-1.));
}

main() {
  float P[Nintc],intP,alphad,alphap,pcmcpre,pcmcpost,r,T,pcmcfoot,pcmctop;
  float Gammap,Gammad,c,cfoot,ctop,cmin;
  float rhobar,sigmap,taueff,UP,DOWN,synchange;
  int ic,icmin,ictop,icfoot,nT;
  for(r=1.;r<1.5;r+=1.) {
    /**** post before pre ****/
    /*cadistrib=fopen("pca.dat","w");*/
    for(nT=10;nT<icpost+1;nT+=10) {
      T=tauca*log(cpost/(nT*dc));
      cfoot=nT*dc;
      icfoot=nT;
      cmin=cfoot;
      icmin=nT;
      ctop=cfoot+cpre;
      ictop=icmin+icpre;
      for(ic=1;ic<nT+1;ic++) {
	c=ic*dc;
	P[ic]=Pfirst(c,r);
	/*fprintf(cadistrib,"%f %f %f %f %f %f\n",c,pcmcpost,pcmcfoot,pcmctop,P[ic],intP); */
      }
      intP=pow(cmin,r*tauca)/(r*tauca);
      if(thetad<cmin) alphad=(pow(cmin,r*tauca)-pow(thetad,r*tauca))/(r*tauca);
      else alphad=0.;
      if(thetap<cmin) alphap=(pow(cmin,r*tauca)-pow(thetap,r*tauca))/(r*tauca);
      else alphap=0.;
      for(ic=icmin+1;ic<Nintc;ic++) {
	c=ic*dc;
	if(ic<icpost+1) pcmcpost=0.;
	else if(ic==icpost+1) pcmcpost=pow(dc/cpost,r*tauca)/(r*tauca);
	else pcmcpost=0.5*dc*(P[(int)((c-cpost)/dc)]/pow(c,r*tauca)+P[(int)((c-cpost)/dc)-1]/pow(c-dc,r*tauca));
	if(ic<icfoot+1) pcmcfoot=0.;
	else if(ic==icfoot+1) pcmcfoot=pow(dc/cfoot,r*tauca)/(r*tauca);
	else pcmcfoot=0.5*dc*(P[(int)((c-cfoot)/dc)]/pow(c,r*tauca)+P[(int)((c-cfoot)/dc)-1]/pow(c-dc,r*tauca));
	if(ic<ictop+1) pcmctop=0.;
	else if(ic==ictop+1) pcmctop=pow(dc/ctop,r*tauca)/(r*tauca);
	else pcmctop=0.5*dc*(P[(int)((c-ctop)/dc)]/pow(c,r*tauca)+P[(int)((c-ctop)/dc)-1]/pow(c-dc,r*tauca));
	P[ic]=pow(c,r*tauca-1)*(P[ic-1]/pow(c-dc,r*tauca-1)-r*tauca*(pcmcpost-pcmcfoot+pcmctop));
	intP+=0.5*dc*(P[ic-1]+P[ic]);
	/*if(ic<200) fprintf(cadistrib,"%f %f %f %f %f %f\n",c,pcmcpost,pcmcfoot,pcmctop,P[ic],intP);*/
	if(c>thetad) alphad+=0.5*dc*(P[ic-1]+P[ic]);
	if(c>thetap) alphap+=0.5*dc*(P[ic-1]+P[ic]);
      }
      alphad/=intP;
      alphap/=intP;
      Gammad=gammad*alphad;
      Gammap=gammap*alphap;
      /*for(ic=1;ic<Nintc;ic++) {
	c=ic*dc;
	printf("%f %f\n",c,P[ic]/intP);
	}
	printf("\n");*/
      /*printf("%d %f %f %f %f %f %f ",nT,cpost,cfoot,ctop,intP,alphad,alphap,Gammap/(Gammad+Gammap));*/
      rhobar=Gammap/(Gammad+Gammap);
      sigmap=sigma*sqrt((alphap+alphad)/(Gammap+Gammad));
      taueff=tau/(Gammap+Gammad);
      UP=0.5*(1.-erf((rhostar-rhobar+rhobar*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
      DOWN=0.5*(1.+erf((rhostar-rhobar+(rhobar-1.)*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
      synchange=((fup*(1.-DOWN)+(1.-fup)*UP)*cmich+fup*DOWN+(1.-fup)*(1.-UP))/(fup*cmich+1.-fup);
      printf("%f %f %f %f\n",-T,UP,DOWN,synchange);
    }
    /*** pre before post ***/
    for(nT=icpre-10;nT>10;nT-=10) {
      T=tauca*log((cpre/(nT*dc)));
      cfoot=nT*dc;
      icfoot=nT;
      cmin=cfoot;
      icmin=nT;
      ctop=cfoot+cpost;
      ictop=icmin+icpost;
      for(ic=1;ic<nT+1;ic++) {
	c=ic*dc;
	P[ic]=Pfirst(c,r);
      }
      intP=pow(cmin,r*tauca)/(r*tauca);
      if(thetad<cmin) alphad=(pow(cmin,r*tauca)-pow(thetad,r*tauca))/(r*tauca);
      else alphad=0.;
      if(thetap<cmin) alphap=(pow(cmin,r*tauca)-pow(thetap,r*tauca))/(r*tauca);
      else alphap=0.;
    /*printf("%f\n",cmin);*/
      for(ic=icmin+1;ic<Nintc;ic++) {
	c=ic*dc;
	if(ic<icpre+1) pcmcpre=0.;
	else if(ic==icpre+1) pcmcpre=pow(dc/cpre,r*tauca)/(r*tauca);
	else pcmcpre=0.5*dc*(P[(int)((c-cpre)/dc)]/pow(c,r*tauca)+P[(int)((c-cpre)/dc)-1]/pow(c-dc,r*tauca));
		  
	if(ic<icfoot+1) pcmcfoot=0.;
	else if(ic==icfoot+1) pcmcfoot=pow(dc/cfoot,r*tauca)/(r*tauca);
	else pcmcfoot=0.5*dc*(P[(int)((c-cfoot)/dc)]/pow(c,r*tauca)+P[(int)((c-cfoot)/dc)-1]/pow(c-dc,r*tauca));
		  
	if(ic<ictop+1) pcmctop=0.;
	else if(ic==ictop+1) pcmctop=pow(dc/ctop,r*tauca)/(r*tauca);
	else pcmctop=0.5*dc*(P[(int)((c-ctop)/dc)]/pow(c,r*tauca)+P[(int)((c-ctop)/dc)-1]/pow(c-dc,r*tauca));
		  
	P[ic]=pow(c,r*tauca-1)*(P[ic-1]/pow(c-dc,r*tauca-1)-r*tauca*(pcmcpre-pcmcfoot+pcmctop));
	intP+=0.5*dc*(P[ic-1]+P[ic]);
	/*printf("%f %f %f %f %f\n",c,pcmcpre,pcmcpost,P[ic],intP);*/
	if(c>thetad) alphad+=0.5*dc*(P[ic-1]+P[ic]);
	if(c>thetap) alphap+=0.5*dc*(P[ic-1]+P[ic]);
      }
      alphad/=intP;
      alphap/=intP;
      Gammad=gammad*alphad;
      Gammap=gammap*alphap;
      /*for(ic=1;ic<Nintc;ic++) {
	c=ic*dc;
	printf("%f %f\n",c,P[ic]/intP);
	}
	printf("\n");*/
      /*printf("%f %f %f %f\n",r,alphad,alphap,Gammap/(Gammad+Gammap));*/
      rhobar=Gammap/(Gammad+Gammap);
      sigmap=sigma*sqrt((alphap+alphad)/(Gammap+Gammad));
      taueff=tau/(Gammap+Gammad);
      UP=0.5*(1.-erf((rhostar-rhobar+rhobar*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
      DOWN=0.5*(1.+erf((rhostar-rhobar+(rhobar-1.)*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
      synchange=((fup*(1.-DOWN)+(1.-fup)*UP)*cmich+fup*DOWN+(1.-fup)*(1.-UP))/(fup*cmich+1.-fup);
      printf("%f %f %f %f\n",T,UP,DOWN,synchange);
    }
  }
}
    
