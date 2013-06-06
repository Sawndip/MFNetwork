#include <math.h>
#include <stdio.h>
#include <stdlib.h> // malloc()

#include "newcv.h"
#include "shotnoise.h"

/* gives nu_e=5.83Hz, nu_i=11.53Hz
 #define JEE 0.1
 #define JEI 0.4
 #define JIE 0.2
 #define JII 0.8
 #define sigma 5. 
 */
#define c_ee (395) //(400) //7999 //400 //1000 //4000 /* connectivity: each neuron receives X incoming connections */
#define c_ie (395) //(400) //8000 //400 //1000 //4000 
#define c_ei (100) //(100) //2000 //100 //250 //1000 
#define c_ii (100) //(100) //1999 //100 //250 //1000

#define c_ss (5)
#define c_es (5)
#define c_se (395)
#define c_is (5)
#define c_si (100)

//#define w_ie 0.2e-3 //(0.2e-3 * 0.05) //0.2e-3 //0.1e-3 //0.2e-3 /* weights */
//#define w_ei 0.4e-3 //(0.4e-3 * 0.05) //0.4e-3 //0.6e-3  
//#define w_ii 0.8e-3 //(0.8e-3 * 0.05) //0.8e-3 //0.6e-3 
//#define w_min 0.0
//#define w_len 0.2e-3 //(0.2e-3 * 0.05) //0.2e-3 //1.0e-3

#define w_ie (0.1e-3) //(0.2e-3 * 0.05) //0.2e-3 //0.1e-3 //0.2e-3 /* weights */
#define w_ei (0.4e-3) //(0.4e-3 * 0.05) //0.4e-3 //0.6e-3  
#define w_ii (0.4e-3) //(0.8e-3 * 0.05) //0.8e-3 //0.6e-3 
#define w_min 0.0
#define w_len (0.2e-3) //(0.2e-3 * 0.05) //0.2e-3 //1.0e-3
#define RHO_INIT 0.164844 //0.019 //0.164855 //0.019 //0.203584 //0.26172
//0.232253
//(0.371888) 4hz in-vitro
//(0.231090) 2hz in-vitro
//(0.203584) 1hz in-vivo
//0.16492 1hz in-vitro
 //0.165171 //0.164871 //0.165171 //0.502946 //0.502964 /*0.185*/ /*0.176923 1hz net*/ //165 //0.35 //0.38 //0.5 //1.

#define w_es
#define w_se
#define w_ss

#define w_is w_ie
#define w_si w_ei

#define RHO_FIXED RHO_INIT //0.502964

#define J_EXT (0.011046) /*(0.0108731)*/ /*(0.0108731)*/ /*(0.006966)*/ /*(0.00707)*/ /*(0.00695)*/ /*(0.006944)*/ /*(0.0113628)*/ /*(0.00695)*/ /*(0.00686) in-vivo*/ /*(0.0072) 1hz network*/ /*(0.0066)*/
#define J_STIM (0.011046)

#define NU_E_INIT (10.2) /*(1.25)*/
#define NU_I_INIT (10.2) /*(1.25)*/
#define NU_S_INIT (10.2)

/*double c_pre = 0.33705;//0.337;//0.562;//0;//5; //0.33705; //0.5617539;
double c_post = 0.74378;//0.744; //1.24; //7;//8; //0.74378; //1.23964;*/
#define cpre 0.56175 //0.33705 //0.56175
#define cpost 1.23964 //0.74378 //1.23964

#define tau_e 0.02 //0.01 /* excitatory population time constant (seconds) */
#define tau_i 0.02 //0.01 /* inhibitory population time constant (seconds) */
#define tau_me 0.02 /* excitatory membrane time constant (seconds) */
#define tau_mi 0.02 //0.02 //0.01 /* inhibitory membrane time constant (seconds) */
#define theta 0.02 //0.016 //0.02 //0.016 //0.02 /* threshold potential */
#define v_r 0.01 //0.002 //0.01 //0.02 //0.01 //0. /* reset potential */
#define sigma_ext 5.e-3 //0.5e-3 //5.e-3 //0.5e-3 /* noise */

#define tau_s 0.02 /* s-population time constant (seconds)*/
#define tau_ms 0.02 /* s-population membrane time constant (seconds)*/

#define tmax 3.0 //30. //0.1 //3.1002 //10. /* seconds */
#define dt 0.0001
#define dwt 0.100
#define NintT ((int)(tmax/dt))
#define wNintT ((int)(tmax/dwt))
#define mfNintT ((int)(dwt/dt))
#define CONVERGENCE_CRITERION (0.00000001)
/* tau_rp defined in newcv.c */


FILE *fopen(),*output_file;


int main(void) {
	double mu_e, mu_i;

	double *nu_e = malloc((NintT+1) * sizeof(double));
	double *nu_i = malloc((NintT+1) * sizeof(double));
	double *w_ee = malloc((NintT+1) * sizeof(double));
	double *rho = malloc((wNintT+1) * sizeof(double));
	
	int it, jt;
	int mt, nt;
	
	rho[0] = RHO_INIT;
	//TODO: reenable weight update
	//w_ee[0] = w_min + (w_len * rho[0]); //0.1e-3; //0.1e-3; //0.2e-3; //0.7;
	w_ee[0] = w_min + (w_len * RHO_FIXED);
	nu_e[0] = NU_E_INIT; //5.25; /* excitatory population initial rate */
	nu_i[0] = NU_I_INIT; //5.25; /* inhibitory population initial rate */
	
	mu_e = J_EXT; //0.02; //0.025; /* external input to excitatory population */
	mu_i = J_EXT; //0.02; //0.025; /* external input to inhibitory population */
	
	if(! (output_file = fopen("output_MF.dat", "a")) ){
		perror("output_MF.dat, error opening output file.");
		return 1;
	}
	
	it=0; // Mean-field index variable
	jt=0; // Synaptic efficacy index variable
	
	printf("Inputs: mu_e: %gV, mu_i: %gV. Equivalent to: nu_ext_e: %gHz, nu_ext_i: %gHz\n", mu_e, mu_i, (mu_e/(c_ee * w_ee[0] * tau_me)), (mu_i/(c_ie * w_ie * tau_me)));
	printf("time: 0.00, nu_e[0]: %f, nu_i[0]: %f, w_ee[0]: %g\n", nu_e[it], nu_i[it], w_ee[it]);
	
	fprintf(output_file, "%f %f %f %f %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
	//for(it = 1; it < NintT; it++){
	for(mt = 0; mt < wNintT; mt++){
		double phi_mu_e, phi_mu_i;
		double x_local, y_local;
		double e_trans, i_trans;
		double d_nu_e, d_nu_i;
		double convergence = 1000.;
		double sigma_e, sigma_i;
		
		for(nt = 1; (nt <= mfNintT); nt++){
			// Update Mean-field equations, until next scheduled synapse update (or stop early upon convergence)
			it = (mt * mfNintT) + nt;
			
			//TODO: add a third population to account for stimulated subpopulation
			phi_mu_e = mu_e + (c_ee * w_ee[it-1] * tau_me * nu_e[it-1]) - (c_ei * w_ei * tau_me * nu_i[it-1]);
			//phi_mu_e = phi_mu_e + (5. * (0.8 * w_len) * tau_me * 2);
			phi_mu_i = mu_i + (c_ie * w_ie * tau_mi * nu_e[it-1]) - (c_ii * w_ii * tau_mi * nu_i[it-1]);
			printf("pme: %f, pmi: %f, ", phi_mu_e, phi_mu_i);
			
			sigma_e = sigma_ext + (c_ee * w_ee[it-1] * w_ee[it-1] * tau_me * nu_e[it-1]) + (c_ei * w_ei * w_ei * tau_me * nu_i[it-1]);
			sigma_i = sigma_ext + (c_ie * w_ie * w_ie * tau_mi * nu_e[it-1]) + (c_ii * w_ii * w_ii * tau_mi * nu_i[it-1]);
			//TODO: enable only external sigma values here
			//sigma_e = sigma_ext;
			//sigma_i = sigma_ext;
			printf("sig_e: %f, sig_i: %f, ", sigma_e, sigma_i);
			
			x_local = (theta - phi_mu_e)/sigma_e;
			y_local = (v_r - phi_mu_e)/sigma_e;
			e_trans = trans(x_local, y_local, tau_me);
			d_nu_e = (-nu_e[it-1] + e_trans) / tau_e;
			// Original update function (stable)
			nu_e[it] = nu_e[it-1] + (dt * d_nu_e);
			// TODO: Quick and dirty jump to transfer function value
			//nu_e[it] = e_trans;
			printf("x: %.2f, y: %.2f, trans: %.6f, de: %.2f, ", x_local, y_local, e_trans, d_nu_e);
			//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, e_trans);
		
			x_local = (theta - phi_mu_i)/sigma_e;
			y_local = (v_r - phi_mu_i)/sigma_e;
			i_trans = trans(x_local, y_local, tau_mi);
			d_nu_i = (-nu_i[it-1] + i_trans) / tau_i;
			// Original update function (stable)
			nu_i[it] = nu_i[it-1] + (dt * d_nu_i);
			// TODO: Quick and dirty jump to transfer function value
			//nu_i[it] = i_trans;
			printf("x: %.2f, y: %.2f, trans: %.6f, di: %.2f, ", x_local, y_local, i_trans, d_nu_i);
			//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);
			
			if(nt < mfNintT){ // On last pass through loop only update MF equations, don't update weight or output to file
				w_ee[it] = w_ee[it-1];
				
				printf("time: %g, nu_e[%d]: %.6f, nu_i[%d]: %.6f, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
				fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
				
				// Debug these guys
				//printf("d_nu_e: %0.8f, d_nu_i: %0.8f, ", d_nu_e, d_nu_i);
				
				convergence = fmax(fabs(d_nu_e), fabs(d_nu_i));
				if(convergence < CONVERGENCE_CRITERION){
					printf("convergence at %f, it: %d\n", convergence, it);
					nu_e[((mt+1)*mfNintT)] = nu_e[it];
					nu_i[((mt+1)*mfNintT)] = nu_i[it];
					w_ee[((mt+1)*mfNintT)-1] = w_ee[it-1];
					it = ((mt+1)*mfNintT);
					nt = mfNintT;
				}
			}
		}
		
		//if (it % mfNintT == 0){
		jt++;
		/* update weights here */
		// updateWeight(float rho_old, float stepsize, float rate, float c_pre, float c_post)
		//TODO: reenable weight updates
		//rho[jt] = updateWeight(rho[jt-1], dwt, nu_e[it], cpre, cpost);
		rho[jt] = rho[jt-1];
		//TODO: reenable weight updates
		//w_ee[it] = w_min + (w_len * rho[jt]); // -1 due to extra increment on loop above
		w_ee[it] = w_min + (w_len * RHO_FIXED);
		/*}
		else{
			w_ee[it] = w_ee[it-1];
		}*/
		if(fabs(rho[jt-1] - rho[jt]) < 0.00001){
			printf("BAZINGA!\n");
		}
		
		printf("time: %g, nu_e[%d]: %g, nu_i[%d]: %g, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
		fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
	}
	
	fprintf(output_file, "\n\n\n\n\n");
	fclose(output_file);
}
    
