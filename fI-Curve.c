/*
 *  fI-Curve.c
 *  MFNetwork
 *
 *  Created by David Higgins on 13/03/2013.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h> // malloc()
#include <math.h> // fmax()

#include "newcv.h"
#include "shotnoise.h"

/* gives nu_e=5.83Hz, nu_i=11.53Hz
 #define JEE 0.1
 #define JEI 0.4
 #define JIE 0.2
 #define JII 0.8
 #define sigma 5. 
 */
#define c_ee (0) //400 //400 //7999 //400 //1000 //4000 /* connectivity: each neuron receives X incoming connections */
#define c_ie (0) //400 //8000 //400 //1000 //4000
#define c_ei (0) //100 //2000 //100 //250 //1000
#define c_ii (0) //100 //1999 //100 //250 //1000

#define w_ie (0.1e-3) //(0.2e-3 * 0.05) //0.2e-3 //0.1e-3 //0.2e-3 /* weights */
#define w_ei (0.4e-3) //(0.4e-3 * 0.05) //0.4e-3 //0.6e-3  
#define w_ii (0.4e-3) //(0.8e-3 * 0.05) //0.8e-3 //0.6e-3 
#define w_min 0.0
#define w_len (0.2e-3) //(0.2e-3 * 0.05) //0.2e-3 //1.0e-3

#define RHO_INIT 0.164844 //0.019 //0.164855 //0.019 //0.203584 //0.26172

#define RHO_FIXED 0.5 //0.164871


#define J_EXT (0.006944) /*(0.00695)*/ /*(0.0113628)*/ /*(0.00695)*/ /*(0.00686) in-vivo*/ /*(0.0072) 1hz network*/ /*(0.0066)*/
#define NU_E_INIT (1) /*(1.25)*/
#define NU_I_INIT (1) /*(1.25)*/

/*double c_pre = 0.33705;//0.337;//0.562;//0;//5; //0.33705; //0.5617539;
 double c_post = 0.74378;//0.744; //1.24; //7;//8; //0.74378; //1.23964;*/
#define cpre 0.56175 //0.33705 //0.56175
#define cpost 1.23964 //0.74378 //1.23964

#define tau_e 0.02 //0.01 /* excitatory population time constant (seconds) */
#define tau_i 0.02 //0.01 /* inhibitory population time constant (seconds) */
#define tau_me 0.02 /* excitatory membrane time constant (seconds) */
#define tau_mi 0.02 //0.02 //0.01 /* inhibitory membrane time constant (seconds) */
#define theta 0.02 //0.016 //0.02 //0.016 //0.02 /* threshold potential */
#define v_r 0.01 //0.01 //0.002 //0.01 //0.02 //0.01 //0. /* reset potential */
#define sigma_ext 5.e-3 //0.5e-3 //5.e-3 //0.5e-3 /* noise */

#define tmax 3.0 //30. //0.1 //3.1002 //10. /* seconds */
#define dt 0.0001
#define dwt 0.100
#define NintT ((int)(tmax/dt))
#define wNintT ((int)(tmax/dwt))
#define mfNintT ((int)(dwt/dt))
#define CONVERGENCE_CRITERION (0.000001)
/* tau_rp defined in newcv.c */


FILE *fopen(),*output_file;


int main(void) {
	printf("hello\n");
	
	
	//TODO: begin loop over mu_tot with fixed synapses from here
	if(! (output_file = fopen("simple_fI_MF.dat", "a")) ){
		perror("simple_fI_MF.dat, error opening output file.");
		return 1;
	}
	
	// Simple frequency as a function of total drive (\mu)
	//for(float mu_tot = 0; mu_tot < 31e-3; mu_tot += 0.3e-3){
    for(float mu_tot = 0.005; mu_tot < 16e-3; mu_tot += 0.1e-3){
		float x_local, y_local, e_trans;
		float sigma_e, sigma_i;
	 
		//sigma_e = sigma_ext + (c_ee * w_ee * w_ee * tau_me * nu_e[it-1]) + (c_ei * w_ei * w_ei * tau_me * nu_i[it-1]);
		//sigma_i = sigma_ext + (c_ie * w_ie * w_ie * tau_mi * nu_e[it-1]) + (c_ii * w_ii * w_ii * tau_mi * nu_i[it-1]);
		//TODO: enable only external sigma values here
		sigma_e = sigma_ext;
		sigma_i = sigma_ext;
		
		x_local = (theta - mu_tot)/sigma_e;
		y_local = (v_r - mu_tot)/sigma_e;
		e_trans = trans(x_local, y_local, tau_me);
		printf("x: %.2f, y: %.2f, trans: %.6f\n", x_local, y_local, e_trans);
		fprintf(output_file, "%f %f\n", mu_tot, e_trans);
	}
	

	
	//TODO: loop over mu_ext and rho from here
	/*if(! (output_file = fopen("multi_fI_MF.dat", "a")) ){
		perror("multi_fI_MF.dat, error opening output file.");
		return 1;
	}
	

	float *nu_e = malloc((NintT+1) * sizeof(float));
	float *nu_i = malloc((NintT+1) * sizeof(float));
	
	
	nu_e[0] = NU_E_INIT; //5.25;
	nu_i[0] = NU_I_INIT; //5.25; 
	
	
	fprintf(output_file, "# mu_ext, e_trans, nu_e[NintT], nu_i[NintT], mu_tot_e, mu_tot_i\n");
	// Frequency as a function of external drive, requires convergence. Try alternate rho values.
	//for (float mu_ext = 0; mu_ext <= 14e-3; mu_ext += 0.5e-3){
    for (float mu_ext = 15.315e-3; mu_ext <= 16e-3; mu_ext += 05.5e-3){
		float mu_tot_e, mu_tot_i;
		float x_local, y_local;
		float e_trans, i_trans;
		float d_nu_e, d_nu_i;
		float convergence = 1000.;		
		float w_ee;
		float sigma_e, sigma_i;
		
		
		sigma_e = sigma_ext;
		sigma_i = sigma_ext;
		
		// Loop over rho values
		for (float rho = 0; rho <= 1; rho+=0.02){
			w_ee = (rho * 0.2e-3);
			
			// Find nu_final for this mu_ext
			for (int it = 1; it <= NintT; it++){
				// Update drives
				mu_tot_e = mu_ext + (c_ee * w_ee * tau_me * nu_e[it-1]) - (c_ei * w_ei * tau_me * nu_i[it-1]);
				mu_tot_i = mu_ext + (c_ie * w_ie * tau_mi * nu_e[it-1]) - (c_ii * w_ii * tau_mi * nu_i[it-1]);
				
				sigma_e = sigma_ext + (c_ee * w_ee * w_ee * tau_me * nu_e[it-1]) + (c_ei * w_ei * w_ei * tau_me * nu_i[it-1]);
				sigma_i = sigma_ext + (c_ie * w_ie * w_ie * tau_mi * nu_e[it-1]) + (c_ii * w_ii * w_ii * tau_mi * nu_i[it-1]);
				//TODO: enable only external sigma values here
				//sigma_e = sigma_ext;
				//sigma_i = sigma_ext;
				//printf("pme: %f, pmi: %f, ", mu_tot_e, mu_tot_i);
		
				// Excitatory population
				x_local = (theta - mu_tot_e)/sigma_e;
				y_local = (v_r - mu_tot_e)/sigma_e;
				e_trans = trans(x_local, y_local, tau_me);
				d_nu_e = (-nu_e[it-1] + e_trans) / tau_e;
				// Original update function (stable)
				nu_e[it] = nu_e[it-1] + (dt * d_nu_e);
				// TODO: Quick and dirty jump to transfer function value
				//nu_e[it] = e_trans;
				//printf("x: %.2f, y: %.2f, trans: %.6f, de: %.2f, ", x_local, y_local, e_trans, d_nu_e);
				//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, e_trans);
		
				// Inhibitory population
				x_local = (theta - mu_tot_i)/sigma_i;
				y_local = (v_r - mu_tot_i)/sigma_i;
				i_trans = trans(x_local, y_local, tau_mi);
				d_nu_i = (-nu_i[it-1] + i_trans) / tau_i;
				// Original update function (stable)
				nu_i[it] = nu_i[it-1] + (dt * d_nu_i);
				// TODO: Quick and dirty jump to transfer function value
				//nu_i[it] = i_trans;
				//printf("x: %.2f, y: %.2f, trans: %.6f, di: %.2f, ", x_local, y_local, i_trans, d_nu_i);
				//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);

			
			}
			convergence = fmax(fabs(d_nu_e), fabs(d_nu_i));
			if (convergence < CONVERGENCE_CRITERION){
				printf("Converged\n");
			}
			else{
				printf("Not converged: %f, %f\n", mu_ext, convergence);
			}
		
			fprintf(output_file, "%f %f %f %f %f %f %f %f %f %f\n", mu_ext, e_trans, nu_e[NintT], nu_i[NintT], w_ee, rho, mu_tot_e, mu_tot_i, sigma_e, sigma_i);
		}
	}*/
	
	
	
	//TODO: begin loop over mu_ext with plastic synapses from here
	/*if(! (output_file = fopen("plastic_fI_MF.dat", "a")) ){
		perror("plastic_fI_MF.dat, error opening output file.");
		return 1;
	}
	
	
	float *nu_e = malloc((NintT+1) * sizeof(float));
	float *nu_i = malloc((NintT+1) * sizeof(float));
	float *w_ee = malloc((NintT+1) * sizeof(float));
	float *rho = malloc((wNintT+1) * sizeof(float));

	fprintf(output_file, "# mu_ext, e_trans, nu_e[NintT], nu_i[NintT], mu_tot_e, mu_tot_i\n");
	// Frequency as a function of external drive, requires convergence
	for (float mu_ext = 0.1e-3; mu_ext <= 12e-3; mu_ext += 1e-3) //0.5e-3)
	{
		float mu_tot_e, mu_tot_i;
		float x_local, y_local;
		float e_trans, i_trans;
		float d_nu_e, d_nu_i;
		float convergence = 1000.;		
		float sigma_e, sigma_i;
		
		int it, jt;
		
		it=jt=0;
		
		//float mu_ext = 0.1e-3;
		
		rho[0] = RHO_INIT;
		w_ee[0] = w_min + (w_len * rho[0]);
		
		sigma_e = sigma_ext;
		sigma_i = sigma_ext;
		
		nu_e[0] = NU_E_INIT; //5.25;
		nu_i[0] = NU_I_INIT; //5.25; 
	
		
		for(int mt = 0; mt < wNintT; mt++){
			//float phi_mu_e, phi_mu_i;
			//float x_local, y_local;
			//float e_trans, i_trans;
			//float d_nu_e, d_nu_i;
			convergence = 1000.;
			
			for(int nt = 1; (nt <= mfNintT); nt++){
				// Update Mean-field equations, until next scheduled synapse update (or stop early upon convergence)
				it = (mt * mfNintT) + nt;
				
				// Update drives
				mu_tot_e = mu_ext + (c_ee * w_ee[it-1] * tau_me * nu_e[it-1]) - (c_ei * w_ei * tau_me * nu_i[it-1]);
				mu_tot_i = mu_ext + (c_ie * w_ie * tau_mi * nu_e[it-1]) - (c_ii * w_ii * tau_mi * nu_i[it-1]);
				
				//sigma_e = sigma_ext + (c_ee * w_ee[it-1] * w_ee[it-1] * tau_me * nu_e[it-1]) + (c_ei * w_ei * w_ei * tau_me * nu_i[it-1]);
				//sigma_i = sigma_ext + (c_ie * w_ie * w_ie * tau_mi * nu_e[it-1]) + (c_ii * w_ii * w_ii * tau_mi * nu_i[it-1]);
				//TODO: enable only external sigma values here
				sigma_e = sigma_ext;
				sigma_i = sigma_ext;
				
				//printf("pme: %f, pmi: %f, ", mu_tot_e, mu_tot_i);
				
				x_local = (theta - mu_tot_e)/sigma_e;
				y_local = (v_r - mu_tot_e)/sigma_e;
				e_trans = trans(x_local, y_local, tau_me);
				d_nu_e = (-nu_e[it-1] + e_trans) / tau_e;
				// Original update function (stable)
				nu_e[it] = nu_e[it-1] + (dt * d_nu_e);
				// TODO: Quick and dirty jump to transfer function value
				//nu_e[it] = e_trans;
				//printf("x: %.2f, y: %.2f, trans: %.6f, de: %.2f, ", x_local, y_local, e_trans, d_nu_e);
				//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, e_trans);
				
				x_local = (theta - mu_tot_i)/sigma_i;
				y_local = (v_r - mu_tot_i)/sigma_i;
				i_trans = trans(x_local, y_local, tau_mi);
				d_nu_i = (-nu_i[it-1] + i_trans) / tau_i;
				// Original update function (stable)
				nu_i[it] = nu_i[it-1] + (dt * d_nu_i);
				// TODO: Quick and dirty jump to transfer function value
				//nu_i[it] = i_trans;
				//printf("x: %.2f, y: %.2f, trans: %.6f, di: %.2f, ", x_local, y_local, i_trans, d_nu_i);
				//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);
				
				if(nt < mfNintT){ // On last pass through loop only update MF equations, don't update weight or output to file
					w_ee[it] = w_ee[it-1];
					
					//printf("time: %g, nu_e[%d]: %.6f, nu_i[%d]: %.6f, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
					//fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
					
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
			// update weights here /
			// updateWeight(float rho_old, float stepsize, float rate, float c_pre, float c_post)
			//TODO: reenable weight updates
			rho[jt] = updateWeight(rho[jt-1], dwt, nu_e[it], cpre, cpost);
			//rho[jt] = rho[jt-1];
			//TODO: reenable weight updates
			w_ee[it] = w_min + (w_len * rho[jt]); // -1 due to extra increment on loop above
			//w_ee[it] = w_min + (w_len * RHO_FIXED);
			//}
			 //else{
			 //w_ee[it] = w_ee[it-1];
			 //}/
			if(fabs(rho[jt-1] - rho[jt]) < 0.00001){
				printf("BAZINGA!\n");
			}
			
			//printf("time: %g, nu_e[%d]: %g, nu_i[%d]: %g, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
			//fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
		}
		//fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
		printf("%f %f %f %f %f %f %f %f %f %f\n", mu_ext, e_trans, nu_e[NintT], nu_i[NintT], w_ee[NintT], rho[wNintT], mu_tot_e, mu_tot_i, sigma_e, sigma_i);
		fprintf(output_file, "%f %f %f %f %f %f %f %f %f %f\n", mu_ext, e_trans, nu_e[NintT], nu_i[NintT], w_ee[NintT], rho[wNintT], mu_tot_e, mu_tot_i, sigma_e, sigma_i);
		
	}*/
    
    
	
	
	//TODO: begin loop over nu with fixed synapses from here
	/*if(! (output_file = fopen("firing_fI_MF.dat", "a")) ){
	 perror("firing_fI_MF.dat, error opening output file.");
	 return 1;
	 }
	*/
    
	/*double mu_ext = J_EXT;
	double w_ee = RHO_FIXED * w_len;
	*/
	
	/*float *nu_i = malloc((NintT+1) * sizeof(float));
	
	// Loop over nu_e, what is the effect on nu_i
	for(float nu_e = 0.5; nu_e < 1.5; nu_e += 0.05){
		float convergence = 1000.;
		float e_trans;
		float mu_tot_e, mu_tot_i;
		float sigma_e, sigma_i;
		
		nu_i[0] = NU_I_INIT; //5.25; 
		
		//for(float nu_i = 0; nu_i < 25; nu_i += 0.5)
		// Find nu_i final for this nu_e, mu_ext, w_ee combo 
		for (int it = 1; it <= NintT; it++)
		{
			//float nu_i = 1;
			//float mu_tot_e, mu_tot_i;
			float x_local, y_local, i_trans;
			float d_nu_e, d_nu_i;
	 
			// Update drives
			mu_tot_e = mu_ext + (c_ee * w_ee * tau_me * nu_e) - (c_ei * w_ei * tau_me * nu_i[it-1]);
			mu_tot_i = mu_ext + (c_ie * w_ie * tau_mi * nu_e) - (c_ii * w_ii * tau_mi * nu_i[it-1]);
			
			//sigma_e = sigma_ext + (c_ee * w_ee * w_ee * tau_me * nu_e[it-1]) + (c_ei * w_ei * w_ei * tau_me * nu_i[it-1]);
			//sigma_i = sigma_ext + (c_ie * w_ie * w_ie * tau_mi * nu_e[it-1]) + (c_ii * w_ii * w_ii * tau_mi * nu_i[it-1]);
			//TODO: enable only external sigma values here
			sigma_e = sigma_ext;
			sigma_i = sigma_ext;
	 
			x_local = (theta - mu_tot_e)/sigma_e;
			y_local = (v_r - mu_tot_e)/sigma_e;
			e_trans = trans(x_local, y_local, tau_me);
			// Holding nu_e constant here!
			//d_nu_e = (-nu_e[it-1] + e_trans) / tau_e;
			// Original update function (forward Euler)
			//nu_e[it] = nu_e[it-1] + (dt * d_nu_e);
			// TODO: Quick and dirty jump to transfer function value
			//nu_e[it] = e_trans;
			//printf("x: %.2f, y: %.2f, trans: %.6f\n", x_local, y_local, e_trans);
			//fprintf(output_file, "%f %f %f ", nu_e, nu_i[it], e_trans);
			
			x_local = (theta - mu_tot_i)/sigma_i;
			y_local = (v_r - mu_tot_i)/sigma_i;
			i_trans = trans(x_local, y_local, tau_mi);
			d_nu_i = (-nu_i[it-1] + i_trans) / tau_i;
			// Original update function (stable)
			nu_i[it] = nu_i[it-1] + (dt * d_nu_i);
			// TODO: Quick and dirty jump to transfer function value
			//nu_i[it] = i_trans;
			//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);
			//fprintf(output_file, "%f\n", i_trans);

			//printf("x: %.2f, y: %.2f, trans: %.6f, di: %.2f, ", x_local, y_local, i_trans, d_nu_i);
			//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);
			
			convergence = fmax(fabs(d_nu_e), fabs(d_nu_i));
			if (convergence < CONVERGENCE_CRITERION){
				//printf("Converged\n");
			}
			else{
				//printf("Not converged: %f, %f\n", mu_ext, convergence);
			}
		}
		printf("%f %f %f %f\n", nu_e, e_trans, nu_i[NintT], convergence);
		fprintf(output_file, "%f %f %f %f %f %f %f %f\n", nu_e, e_trans, nu_i[NintT], convergence, mu_tot_e, mu_tot_i, sigma_e, sigma_i);
	}
     */
	
    /*
	printf("goodbye 1\n");
	fprintf(output_file,"\n\n\n\n\n");
	 
	free(nu_i);
     */
	
	/*
	double *nu_e = malloc((NintT+1) * sizeof(double));
	// Loop over nu_i, what is the effect on nu_e
	for(double nu_i = 1.0; nu_i < 1.5; nu_i += 0.005)
	{
		//double nu_i = 2.313;
		float convergence = 1000.;
		double i_trans;
		double mu_tot_e, mu_tot_i;
		double sigma_e, sigma_i;
		
		nu_e[0] = 1.1;nu_i; //NU_E_INIT; //5.25; 
		
		
		//fprintf(output_file, "%f %f %f %f %f %f %f %f\n", nu_i, i_trans, nu_e[NintT], convergence, mu_tot_e, mu_tot_i, sigma_e, sigma_i);
		
		//for(float nu_i = 0; nu_i < 25; nu_i += 0.5)
		// Find nu_i final for this nu_e, mu_ext, w_ee combo 
		for (int it = 1; it <= NintT; it++)
		{
			//float nu_i = 1;
			//float mu_tot_e, mu_tot_i;
			double x_local, y_local, e_trans;
			double d_nu_e, d_nu_i;
			
			// Update drives
			mu_tot_e = mu_ext + (c_ee * w_ee * tau_me * nu_e[it-1]) - (c_ei * w_ei * tau_me * nu_i);
			mu_tot_i = mu_ext + (c_ie * w_ie * tau_mi * nu_e[it-1]) - (c_ii * w_ii * tau_mi * nu_i);
			
			//sigma_e = sigma_ext + (c_ee * w_ee * w_ee * tau_me * nu_e[it-1]) + (c_ei * w_ei * w_ei * tau_me * nu_i[it-1]);
			//sigma_i = sigma_ext + (c_ie * w_ie * w_ie * tau_mi * nu_e[it-1]) + (c_ii * w_ii * w_ii * tau_mi * nu_i[it-1]);
			//TODO: enable only external sigma values here
			sigma_e = sigma_ext;
			sigma_i = sigma_ext;
			
			x_local = (theta - mu_tot_e)/sigma_e;
			y_local = (v_r - mu_tot_e)/sigma_e;
			e_trans = trans(x_local, y_local, tau_me);
			d_nu_e = (-nu_e[it-1] + e_trans) / tau_e;
			// Original update function (forward Euler)
			nu_e[it] = nu_e[it-1] + (dt * d_nu_e);
			if(nu_e[it] < 0)
				nu_e[it] = 0;
			// TODO: Quick and dirty jump to transfer function value
			//nu_e[it] = e_trans;
			//printf("x: %.2f, y: %.2f, trans: %.6f\n", x_local, y_local, e_trans);
			//fprintf(output_file, "%f %f %f ", nu_e, nu_i[it], e_trans);
			
			x_local = (theta - mu_tot_i)/sigma_i;
			y_local = (v_r - mu_tot_i)/sigma_i;
			i_trans = trans(x_local, y_local, tau_mi);
			// Holding nu_i constant here!
			//d_nu_i = (-nu_i[it-1] + i_trans) / tau_i;
			// Original update function (stable)
			//nu_i[it] = nu_i[it-1] + (dt * d_nu_i);
			// TODO: Quick and dirty jump to transfer function value
			//nu_i[it] = i_trans;
			//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);
			//fprintf(output_file, "%f\n", i_trans);
			
			//printf("x: %.2f, y: %.2f, trans: %.6f, di: %.2f, ", x_local, y_local, i_trans, d_nu_i);
			//printf("x: %.2f, y: %.2f, trans: %.2f, ", x_local, y_local, i_trans);
			
			convergence = fmax(fabs(d_nu_e), fabs(d_nu_i));
			if (convergence < CONVERGENCE_CRITERION){
				//printf("Converged\n");
			}
			else{
				//printf("Not converged: %f, %f\n", mu_ext, convergence);
			}
			//printf("%f %f %f %f\n", nu_i, i_trans, nu_e[NintT], convergence);
			//fprintf(output_file, "%f %f %f %f %f %f %f %f %f\n", nu_i, i_trans, nu_e[it], convergence, mu_tot_e, mu_tot_i, sigma_e, sigma_i, e_trans);
		}
		printf("%f %f %f %f\n", nu_i, i_trans, nu_e[NintT], convergence);
		fprintf(output_file, "%f %f %f %f %f %f %f %f\n", nu_i, i_trans, nu_e[NintT], convergence, mu_tot_e, mu_tot_i, sigma_e, sigma_i);
	 }*/
	
	
	
	fprintf(output_file, "\n\n\n\n\n");
	fclose(output_file);
	printf("goodbye\n");
}


