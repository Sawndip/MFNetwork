#include <math.h>
#include <stdio.h>
#include <stdlib.h> // malloc()

#include "shotnoise.h"

#define EPSILLON 0.0000001 /* variable added to floats to be cast as ints in order to avoid numerical problems */
#define tauca 0.0226936 /*0.0227*/
//#define cpre 0.56 /* moved to local variables to allow changing via function argument */
//#define cpost 1.24
#define thetad 1
#define thetap 1.3
#define gammad 331.909
#define gammap 725.085
#define sigma 0. /*3.3501*/
#define tau 346.3615 /* 346 */
#define rhostar 0.5
//#define D 0.0046 /* unused */
#define fup 0.5
#define cmich 5.41
/* cmax=30 makes the computations somewhat safe up to r=100Hz */
#define cmax 30.
#define dc 0.00001 /*0.001*/
#define Nintc (int)(cmax/dc) //EPSILLON

// Setup intervals for analytical solutions
#define cfirst ((cpre<cpost)?cpre:cpost) /*previously called cmin*/
#define icfirst (int)(cfirst/dc)//EPSILLON
#define icpre (int)(cpre/dc)//EPSILLON
#define csecondpre ((1.5*cpre<cpre+cpost)?1.5*cpre:cpre+cpost)
#define isecondpre (int)(csecondpre/dc)//EPSILLON
#define icpost (int)(cpost/dc)//EPSILLON
#define csecondpost ((1.5*cpost<cpre+cpost)?1.5*cpost:cpre+cpost)
#define isecondpost (int)(csecondpost/dc)//EPSILLON

/*#define cmin ((cpre<cpost)?cpre:cpost)
#define icpre (int)(cpre/dc)
#define icpost (int)(cpost/dc)
#define icmin (int)(cmin/dc)
//Additions for hypergeometric solution
#define cother ((cpre<cpost)?cpost:cpre)
//#define cimin (((2*cmin)<cother)?(2*cmin):cother)
//#define c2min (cmin + 0.5*(cimin - cmin))
#define c2min (((2*cmin)<cother)?(2*cmin):cother)
#define ic2min (int)(c2min/dc)*/

#define LARGE_NO 20 // factorial(150) is upperlimit for double type so keep this below 150


char filename[100];

// Hypergeometric function
//double HyperFunc21(float a, float b, float c, float z){
double myHyperFunc21(float a, float z){ //a=2*r*tauca, z = cpre-c/cpre
	// Only static arrays (al,bl,cl) would help with efficiency here
	double al = 1.;//, bl = 1., cl = 1.;
	double nl = 1.;
	double coef = 1.0;
	double s;
	int n;
	
	s = 1.; // n=0 case
	
	for(n = 1; n < LARGE_NO /*&& fabs(coef)>0.000001*/; n++){
		//TODO: the following loop is not really necessary, is it?
		//for(; i < (n+1); i++){ // Reminder: i always exits this loop with i=(n+1)
			//printf("n: %d, i: %d, ", n, i);
		
		coef = al * pow(a, 2.) * pow(z, n) / ((a + n) * nl);
		s += coef;
		al *= (a + n);
		nl *= (n + 1);
		
		/*al *= (a + n - 1);
		bl *= (b + n - 1);
		cl *= (c + n - 1);
		//}
		
		nl *= n; // running computation of factorial(n)
		coef = ((al * bl) / cl) * (pow(z, n) / nl);
		s += coef;*/
		//printf("nl: %ld, s: %f\n", nl, (float)s); 
	}
		
	//printf("hyperfn: n final: %d, s: %f\n", n, s);
	
	//printf("s: %f, ", (float)s);
	return s;
}

// Analytical solution for first interval [0,cfirst]
double Pfirst(double c, double r) {
	//printf("Pfirst: %f ", pow(c,2.*r*tauca-1.) );
	return(pow(c,2.*r*tauca-1.));
}

// Analytical solution for second interval (cmin,c2min]
/*float Psecond(float c, float r, float cpre, float cpost){
	float a, b, d;
	
	a = pow(c, 2.*r*tauca-1.);
	b = pow((c-cfirst)/cfirst, (2.*r*tauca)) / 2.; //pow((c-cmin), (2.*r*tauca-1));
	//d = HyperFunc21(2*r*tauca, 2*r*tauca, 2*r*tauca+1, 1-c);
	d = HyperFunc21(2*r*tauca, 2*r*tauca, 2*r*tauca + 1, (cfirst-c)/cfirst);
	
	return (a * (1 - b * d));
}*/

double updateWeight(double rho_old, double stepsize, double rate, double c_pre, double c_post){
	double rho_new, d_rho, alphas[2];
	d_rho = 0.;
	
	getAlphas(rate, c_pre, c_post, alphas);
	
	d_rho = ((gammap * alphas[1]) - rho_old * (gammap * alphas[1] + gammad * alphas[0])) / tau;
	rho_new = rho_old + (stepsize * d_rho);
	printf("rho_old: %g, d_rho: %g, rho_new: %g\n", rho_old, d_rho, rho_new);
	printf("alphad: %g, alphap: %g\n", alphas[0], alphas[1]);
	return rho_new;
}

void getAlphas(double rate, double c_pre, double c_post, double *alphas){
	/* This function has been modified to handle a zero firing rate by returning alphas=0,
		it does not necessarily handle other variables correctly in that case */
	double intP,alphad,alphap,r, Qpre, Qpost;//pcmcpre,pcmcpost
	double Gammap,Gammad,c;
	double rhobar,sigmap,taueff,UP,DOWN,synchange;
	int ic;

	//double P[Nintc]
	double *P = calloc(Nintc, sizeof(double));
	
	Qpre = 0.;
	Qpost = 0.;
	
	double cpre = c_pre;
	double cpost = c_post;
	//float interP; //debugging of a crossover point
	FILE *fp;
	fp = fopen("shot_output.dat", "w");
	//fp = fopen(filename, "w");
	//for(r=1.;r<50.5;r+=1.) {
	if (rate > 0){ /* rate should always be > 0 */
	  r = rate;
	  P[0] = 0.; // avoids problems caused by numerical errors in ((int)((c-cmin)/dc)-1) (probably no longer needed since introduction of EPSILLON)
	  //P[(int)(cfirst/dc + EPSILLON)] = Pfirst(cfirst,r);
	  intP = 0.;
	  intP = pow(cfirst,2.*r*tauca)/(2.*r*tauca);
	  if(thetad < cfirst){
		  alphad = (pow(cfirst,2.*r*tauca)-pow(thetad,2.*r*tauca))/(2.*r*tauca);
	  }
	  else 
		  alphad = 0.;
	
	  if(thetap < cfirst){
		  alphap = (pow(cfirst,2.*r*tauca)-pow(thetap,2.*r*tauca))/(2.*r*tauca);
	  }
	  else 
		  alphap = 0.;
	
	  /*printf("%f\n",cmin);*/
	  // Solve first interval analytically
	  for(ic=1;ic<icfirst+1;ic++) {
		  c = ic*dc;
		  P[ic] = Pfirst(c,r);
		  /*if(( ic % (1) ) == 0){
			printf("P1[%d]: %f, \n", ic, P[ic]);
		   }*/
		  fprintf(fp, "%f %f\n", c, P[ic]);
	  }
	  //printf("P[%d]: %f, intP(%d): %f\n", icmin, P[icmin], icmin, intP);
	  
	  // Solve second interval analytically partly using hypergeometric function
	  for(ic =icfirst+1; ic < Nintc; ic++){
		  c = ic*dc;
		  if (ic <= icpost)
			  Qpost = 0.;
		  else if (ic <= isecondpost)
			  Qpost = pow((c/cpost) - 1, 2*r*tauca) * myHyperFunc21(2*r*tauca, (cpost-c)/cpost) / (2*r*tauca);
		  else
			  Qpost += 0.5 * dc * (P[(int)((c-cpost)/dc)] / pow(c, 2*r*tauca) + P[(int)((c-cpost)/dc) - 1] / pow(c-dc, 2*r*tauca));//EPSILLON
		  
		  if (ic <= icpre)
			  Qpre = 0.;
		  else if (ic <= isecondpre)
			  Qpre = pow((c/cpre) - 1, 2*r*tauca) * myHyperFunc21(2*r*tauca, (cpre-c)/cpre) / (2*r*tauca);
		  else
			  Qpre += 0.5 * dc * (P[(int)((c-cpre)/dc)] / pow(c, 2*r*tauca) + P[(int)((c-cpre)/dc) - 1] / pow(c-dc, 2*r*tauca));//EPSILLON
		  
		  /*if(ic % 100 == 0)
			  printf("ic: %d, Qpre: %f, Qpost: %f\n", ic, Qpre, Qpost);*/
		  
		  P[ic] = Pfirst(c, r) * (1. - (0.5* 2*r * tauca * Qpre) - (0.5* 2*r * tauca * Qpost));

		  if(P[ic] < 0){
			  // P(c) should never be less than 0, reset to zero.
			  /*if(( ic % (1) ) == 0){
				printf("Yowzah, P[%d]: %f\n", ic, P[ic]);
				}*/
			  P[ic] = 0.;
		  }
		  
		  intP += 0.5*dc*(P[ic-1] + P[ic]);

		  if(c > thetad){
			  // c>thetad, so update value of alphad (amount of time above depression threshold)
			  alphad += 0.5*dc*(P[ic-1] + P[ic]);
			  /*if (P[ic] > 0)
				  printf("alphad: %f \n", alphad);*/
		  }
		  if(c > thetap){
			  // c>thetap, so update value of alphap (amount of tima above potentiation threshold)
			  //printf("alphap += %f, P[%d-1]:%f, P[%d]:%f,", (0.5*dc*(P[ic-1] + P[ic])), ic, P[ic-1], ic, P[ic]);
			  alphap += 0.5*dc*(P[ic-1] + P[ic]);
			  /*if (P[ic] > 0)
				  printf(" alphap: %f\n", alphap);*/
		  }
		  /*if(( ic % (1) ) == 0){
				printf("P2[%d]: %f, \n", ic, P[ic]);
		   }*/
		  fprintf(fp, "%f %f\n", c, P[ic]);
	  }
	
		
	if(intP == 0){
		printf("Error: intP==0, you need to increase cmax\n");
	}
	else{
		printf("intP: %lf, ", intP);
		alphad /= intP;
		alphap /= intP;
		//interP /= intP;
		//printf("alphap: %f, alphad: %f, interP: %f\n", alphap, alphad, interP);
	}
	//printf("2)alphap: %f, ", alphap);
	// Gamma values are alpha values weighted by learning rates
    Gammad = gammad * alphad;
    Gammap = gammap * alphap;
	/*for(ic=1;ic<Nintc;ic++) {
      c=ic*dc;
      printf("%f %f\n",c,P[ic]/intP);
     }
     printf("\n");*/
    /*printf("%f %f %f %f\n",r,alphad,alphap,Gammap/(Gammad+Gammap));*/
	if ((Gammad == 0) && (Gammap == 0)){
		// Calcium has not crossed either threshold (thetap and thetad) so no change can occur
		printf("Gammap and Gammad == 0, alphap: %.10lf, alphad: %lf, intP(%d): %lf\n", alphap, alphad, ic, intP);
	}
	else{
		rhobar = Gammap / (Gammad + Gammap);
		sigmap = sigma * sqrt( (alphap + alphad) / (Gammap + Gammad) );
		taueff = tau / (Gammap + Gammad);
		UP = 0.5*(1.-erf((rhostar-rhobar+rhobar*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
		DOWN = 0.5*(1.+erf((rhostar-rhobar+(rhobar-1.)*exp(-75./(r*taueff)))/(sigmap*sqrt(1.-exp(-150./(r*taueff))))));
		synchange = ( (fup * (1.-DOWN) + (1.-fup) * UP) * cmich + fup * DOWN + (1.-fup) * (1.-UP) ) / (fup * cmich + 1. - fup);
		printf("rate: %lf, alphad: %lf, alphap: %lf, rhobar: %lf, UP:%lf,  DOWN:%lf, synchange:%lf\n",r, alphad, alphap, rhobar, UP, DOWN, synchange);
    }
	//}
	if (alphad < 0){ //This case should no longer occur! (when modification for P(I)<0 is enabled above)
		alphad = 0.;
		printf("changing alphad\n");
	}
	if (alphap < 0){
		alphap = 0.;
		printf("changing alphap\n");
	}
	alphas[0] = alphad;
	alphas[1] = alphap;
	//printf("P[%d]: %f, alphas[0]: %f, alphas[1]: %f, interP: %f\n", (ic-1), P[ic-1], alphas[0], alphas[1], interP);
	fprintf(fp, "\n\n\n\n\n");
	fclose(fp);
		
	printf("DEBUG cfirst: %lf, csecondpre: %lf, csecondpost: %lf, icpre: %d, icpost: %d\n", cfirst, csecondpre, csecondpost, icpre, icpost);
  }
  else{ // Rate == 0, no activity, set alphas to 0
	printf("Rate == 0, manually setting alphas to 0.\n");
	alphas[0] = 0.;
	alphas[1] = 0.;
  }
	
	free(P);
}
    
int main(void){
	//float rho = 0.5;
	//float stepsize = 0.01;
	double rate = 0.01;
	rate = 1;
	double c_pre = 0.33705;//0.337;//0.562;//0;//5; //0.33705; //0.5617539;
	double c_post = 0.74378;//0.744; //1.24; //7;//8; //0.74378; //1.23964;
	double alphas[2];
	
	double rhobar;
	
	FILE* fp;	
	strcpy(filename, "fine_rate_dep_alphas.dat");
	fp = fopen(filename, "a");
	fprintf(fp, "#rate alpha_d alpha_p (alpha_d - alpha_p) rhobar GammaD GammaP abs(GammaP-GammaD)\n");
	//for(float i = 0.1; i < 100; i+=1){
	//for(float i = 1.0; i < 1.1; i+=1){
	//for(double i = -4; i < 2.001; i+=0.005){
	for(double i = 1; i < 2.001; i+=10.005){
		//rho = updateWeight(rho, stepsize, rate, c_pre, c_post);
		//printf("i: %d, rho: %f\n", i, rho);
		
		//rate = (float) i;
		//rate = pow(10, i);
		//sprintf(filename, "shot_out_rate_%f.dat", rate);
		printf("outfile: %s\n", filename);
		getAlphas(rate, c_pre, c_post, alphas);
		rhobar = (gammap * alphas[1]) / ((gammap * alphas[1]) + (gammad * alphas[0]));
		printf("i: %f, rate: %lf, alphas[0]: %lf, alphas[1]: %lf\n\n", i, rate, alphas[0], alphas[1]);
		fprintf(fp, "%lf %0.10lf %0.10lf %lf %lf, %0.8lf, %0.8lf, %f\n", rate, alphas[0], alphas[1], (alphas[0]-alphas[1]), rhobar, (alphas[0]*gammad), (alphas[1]*gammap), fabs((alphas[1]*gammap)-(alphas[0]*gammad)));
	}
	fprintf(fp, "\n\n\n\n\n");
	fclose(fp);
	
	printf("cmax: %lf, dc: %lf, Nintc: %lf\n", cmax, dc, (float)Nintc);
	
	return 0;
}

