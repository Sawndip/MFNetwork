#include <math.h>
#include <stdio.h>

#include "shotnoise.h"

#define EPSILLON 0.00001 /* variable added to floats to be cast as ints in order to avoid numerical problems */
#define tauca 0.023 /*0.0227*/
//#define cpre 0.56 /* moved to local variables to allow changing via function argument */
//#define cpost 1.24
#define thetad 1
#define thetap 1.3
#define gammad 332.
#define gammap 725.
#define sigma 3.35
#define tau 346.361 /* 346 */
#define rhostar 0.5
//#define D 0.0046 /* unused */
#define fup 0.5
#define cmich 5.41
#define cmin ((cpre<cpost)?cpre:cpost)
#define cmax 30. /* cmax=30 and dc=0.001 make the computations somewhat safe up to r=100Hz */
#define dc 0.001
#define Nintc (int)(cmax/dc)
#define icpre (int)(cpre/dc)
#define icpost (int)(cpost/dc)
#define icmin (int)(cmin/dc)

char filename[100];

float Pfirst(float c, float r) {
	//printf("Pfirst: %f ", pow(c,2.*r*tauca-1.) );
  return(pow(c,2.*r*tauca-1.));
}

float updateWeight(float rho_old, float stepsize, float rate, float c_pre, float c_post){
	float rho_new, d_rho, alphas[2];
	d_rho = 0.;
	
	getAlphas(rate, c_pre, c_post, alphas);
	
	d_rho = ((gammap * alphas[1]) - rho_old * (gammap * alphas[1] + gammad * alphas[0])) / tau;
	rho_new = rho_old + (stepsize * d_rho);
	printf("rho_old: %g, d_rho: %g, rho_new: %g\n", rho_old, d_rho, rho_new);
	printf("alphad: %g, alphap: %g\n", alphas[0], alphas[1]);
	return rho_new;
}

void getAlphas(float rate, float c_pre, float c_post, float *alphas){
	/* This function has been modified to handle a zero firing rate by returning alphas=0,
		it does not necessarily handle other variables correctly in that case */
  float P[Nintc],intP,alphad,alphap,pcmcpre,pcmcpost,r;
  float Gammap,Gammad,c;
  float rhobar,sigmap,taueff,UP,DOWN,synchange;
  int ic;
	
  float cpre = c_pre;
  float cpost = c_post;
	float interP;
	FILE *fp;
	fp = fopen("shot_output.dat", "w");
	//fp = fopen(filename, "w");
  //for(r=1.;r<50.5;r+=1.) {
  if (rate > 0){ /* rate should always be > 0 */
	r = rate;
	P[0] = 0.; // avoids problems caused by numerical errors in ((int)((c-cmin)/dc)-1)
    P[(int)(cmin/dc + EPSILLON)] = Pfirst(cmin,r);
	intP = 0.;
    intP = pow(cmin,2.*r*tauca)/(2.*r*tauca);
    if(thetad < cmin){
		alphad = (pow(cmin,2.*r*tauca)-pow(thetad,2.*r*tauca))/(2.*r*tauca);
	}
    else 
		alphad = 0.;
	
    if(thetap < cmin){
		alphap = (pow(cmin,2.*r*tauca)-pow(thetap,2.*r*tauca))/(2.*r*tauca);
	}
    else 
		alphap = 0.;
	
    /*printf("%f\n",cmin);*/
    for(ic=1;ic<icmin+1;ic++) {
      c = ic*dc;
      P[ic] = Pfirst(c,r);
	  if(( ic % (1) ) == 0){
		  //printf("P1[%d]: %f, \n", ic, P[ic]);
	  }
	  fprintf(fp, "%f %f\n", c, P[ic]);
    }
	//printf("P[%d]: %f, intP(%d): %f\n", icmin, P[icmin], icmin, intP);
	
    for(ic=icmin+1;ic<Nintc;ic++) {
      c = ic*dc;
	  if(ic<icpre+1){ 
		  pcmcpre=0.;
	  }
      else if(ic==icpre+1){
		  //printf("ic==icpre+1\n");
		  pcmcpre=0.5*pow(dc/cpre,2*r*tauca)/(r*tauca);
		  //printf("pcmcpre: %f, ", pcmcpre);
	  }
      else{
		  pcmcpre=0.5*dc*( P[(int)((c-cpre)/dc + EPSILLON)]/pow(c,2.*r*tauca) + P[(int)((c-cpre)/dc + EPSILLON)-1]/pow(c-dc,2.*r*tauca) );
		  /*if(ic % 58 ==0){
			  printf("pcmcpre: %f, c: %f, cpre: %f, (c-cpre): %f, (c-cpre)/dc: %d, P[0]: %f\n", pcmcpre, c, cpre, (c-cpre), (int)((c-cpre)/dc), P[0]);
		  }*/
      }
		
	  if(ic<icpost+1){ 
		  pcmcpost=0.;
	  }
      else if(ic==icpost+1){
		  //printf("ic==icpost+1\n");
		  pcmcpost=0.5*pow(dc/cpost,2*r*tauca)/(r*tauca);
		  //printf("pcmcpost: %f, ", pcmcpost);
	  }
      else{
		  pcmcpost=0.5*dc*( P[(int)((c-cpost)/dc + EPSILLON)]/pow(c,2.*r*tauca) + P[(int)((c-cpost)/dc + EPSILLON)-1]/pow(c-dc,2.*r*tauca) );
		  //printf("pcmcpost: %f, ", pcmcpost);
	  }
		
	  P[ic] = pow(c,2*r*tauca-1)*(P[ic-1]/pow(c-dc,2*r*tauca-1)-r*tauca*(pcmcpre+pcmcpost));
	  if(P[ic] < 0){
		  /*if(( ic % (1) ) == 0){
			  printf("Yowzah, P[%d]: %f\n", ic, P[ic]);
		  }*/
		  P[ic] = 0.;
	  }
	  fprintf(fp, "%f %f\n", c, P[ic]);
	  //printf("P[%d]: %f ", ic, P[ic]);
	  //printf("comp: %f ", (r*tauca*(pcmcpre+pcmcpost)) );
      intP += 0.5*dc*(P[ic-1] + P[ic]);
	  if(ic == 172){
		  //printf("----> intP[17200]: %f, \n", intP);
		  interP = intP;
	  }
	  //printf("intP: %f\n", intP);
      /*printf("%f %f %f %f %f\n",c,pcmcpre,pcmcpost,P[ic],intP);*/
	  if(c > thetad){ 
		  alphad += 0.5*dc*(P[ic-1] + P[ic]);
	  }
	  if(c > thetap){
		  //printf("alphap += %f, P[%d-1]:%f, P[%d]:%f,", (0.5*dc*(P[ic-1] + P[ic])), ic, P[ic-1], ic, P[ic]);
		  alphap += 0.5*dc*(P[ic-1] + P[ic]);
		  //printf(" alphap: %f\n", alphap);
	  }
	  //printf("Alphap: %f\n", alphap);
	  if(( ic % (1) ) == 0){
		  //printf("P2[%d]: %f, \n", ic, P[ic]);
	  }
    }
	
	//printf("1)alphap: %f, ", alphap);
	if(intP == 0){
		printf("Error: intP==0, you need to increase cmax\n");
	}
	else{
		//printf("intP: %f, ", intP);
		alphad /= intP;
		alphap /= intP;
		interP /= intP;
		//printf("alphap: %f, alphad: %f, interP: %f\n", alphap, alphad, interP);
	}
	//printf("2)alphap: %f, ", alphap);
    Gammad = gammad * alphad;
    Gammap = gammap * alphap;
	/*for(ic=1;ic<Nintc;ic++) {
      c=ic*dc;
      printf("%f %f\n",c,P[ic]/intP);
     }
     printf("\n");*/
    /*printf("%f %f %f %f\n",r,alphad,alphap,Gammap/(Gammad+Gammap));*/
	if ((Gammad == 0) && (Gammap == 0)){
		printf("Gammap and Gammad == 0, alphap: %.10f, alphad: %f, intP(%d): %f\n", alphap, alphad, ic, intP);
	}
	else{
		rhobar = Gammap / (Gammad + Gammap);
		sigmap = sigma * sqrt( (alphap + alphad) / (Gammap + Gammad) );
		taueff = tau / (Gammap + Gammad);
		UP = 0.5*(1.-erf((rhostar-rhobar+rhobar*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
		DOWN = 0.5*(1.+erf((rhostar-rhobar+(rhobar-1.)*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
		synchange = ( (fup * (1.-DOWN) + (1.-fup) * UP) * cmich + fup * DOWN + (1.-fup) * (1.-UP) ) / (fup * cmich + 1. - fup);
		printf("rate: %f, alphad: %f, alphap: %f, rhobar: %f, UP:%f,  DOWN:%f, synchange:%f\n",r, alphad, alphap, rhobar, UP, DOWN, synchange);
    }
	//}
	if (alphad < 0){ //This case should no longer occur! (when modification for P(I)<0 is enabled above)
		//alphad = 0.;
		printf("changing alphad\n");
	}
	if (alphap < 0){
		//alphap = 0.;
		printf("changing alphap\n");
	}
	alphas[0] = alphad;
	alphas[1] = alphap;
	//printf("P[%d]: %f, alphas[0]: %f, alphas[1]: %f, interP: %f\n", (ic-1), P[ic-1], alphas[0], alphas[1], interP);
	fprintf(fp, "\n\n\n\n\n");
	fclose(fp);
  }
  else{ // Rate == 0, no activity, set alphas to 0
	printf("Rate == 0, manually setting alphas to 0.\n");
	alphas[0] = 0.;
	alphas[1] = 0.;
  }
}
    
/*int main(void){
	//float rho = 0.5;
	//float stepsize = 0.01;
	float rate = 0.01;
	float c_pre = 0.56;
	float c_post = 1.24;
	float alphas[2];
	
	float rhobar;
	
	FILE* fp;	
	strcpy(filename, "fine_rate_dep_alphas.dat");
	fp = fopen(filename, "a");
	fprintf(fp, "#rate alpha_d alpha_p (alpha_d - alpha_p) rhobar GammaD GammaP abs(GammaP-GammaD)\n");
	for(float i = 0.; i < 0.5; i+=.001){
		//rho = updateWeight(rho, stepsize, rate, c_pre, c_post);
		//printf("i: %d, rho: %f\n", i, rho);
		
		rate = (float) i;
		//sprintf(filename, "shot_out_rate_%f.dat", rate);
		printf("outfile: %s\n", filename);
		getAlphas(rate, c_pre, c_post, alphas);
		rhobar = (gammap * alphas[1]) / ((gammap * alphas[1]) + (gammad * alphas[0]));
		printf("i: %f, alphas[0]: %f, alphas[1]: %f\n\n", i, alphas[0], alphas[1]);
		fprintf(fp, "%f %f %f %f %f, %f, %f, %f\n", i, alphas[0], alphas[1], (alphas[0]-alphas[1]), rhobar, (alphas[0]*gammad), (alphas[1]*gammap), fabs((alphas[1]*gammap)-(alphas[0]*gammad)));
	}
	fprintf(fp, "\n\n\n\n\n");
	fclose(fp);
	
	printf("cmax: %f, dc: %f, Nintc: %f\n", cmax, dc, (float)Nintc);
	
	return 0;
}*/

