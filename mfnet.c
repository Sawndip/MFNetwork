#include <math.h>
#include <stdio.h>

#include "newcv.h"
#include "shotnoise.h"

/* gives nu_e=5.83Hz, nu_i=11.53Hz
 #define JEE 0.1
 #define JEI 0.4
 #define JIE 0.2
 #define JII 0.8
 #define sigma 5. 
 */
#define c_ee 400 //1000 //4000 /* connectivity: each neuron receives X incoming connections */
#define c_ie 400 //1000 //4000 
#define c_ei 100 //250 //1000 
#define c_ii 100 //250 //1000

#define w_ie 0.2e-3 //0.2e-3 //0.1e-3 //0.2e-3 /* weights */
#define w_ei 0.4e-3 //0.4e-3 //0.6e-3  
#define w_ii 0.8e-3 //0.8e-3 //0.6e-3 
#define w_min 0.0
#define w_len 0.2e-3 //0.1e-3 //1.0e-3
#define RHO_INIT 0.5 //1.

#define J_EXT (0.01)
#define NU_E_INIT (1.25)
#define NU_I_INIT (1.25)

#define cpre 0.56
#define cpost 1.24

#define tau_e 0.01 /* excitatory population time constant (seconds) */
#define tau_i 0.01 //0.01 /* inhibitory population time constant (seconds) */
#define tau_me 0.02 /* excitatory membrane time constant (seconds) */
#define tau_mi 0.01 //0.01 /* inhibitory membrane time constant (seconds) */
#define theta 0.016 //0.02 //0.016 //0.02 /* threshold potential */
#define v_r 0.002 //0.01 //0.02 //0.01 //0. /* reset potential */
#define sigma 5.e-3 //0.5e-3 //5.e-3 //0.5e-3 /* noise */

#define tmax 300.0 //0.1 //3.1002 //10. /* seconds */
#define dt 0.0001
#define dwt 0.100
#define NintT ((int)(tmax/dt))
#define wNintT ((int)(tmax/dwt))
#define mfNintT ((int)(dwt/dt))
/* tau_rp defined in newcv.c */


FILE *fopen(),*output_file;


int main(void) {
	float mu_e, mu_i;
	
	float nu_e[NintT+1], nu_i[NintT+1], w_ee[NintT+1], rho[wNintT+1];
	
	int it, jt;
	int mt, nt;
	
	rho[0] = RHO_INIT;
	w_ee[0] = w_min + (w_len * rho[0]); //0.1e-3; //0.1e-3; //0.2e-3; //0.7;
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
		float phi_mu_e, phi_mu_i;
		float x_local, y_local;
		float e_trans, i_trans;
		float d_nu_e, d_nu_i;
		float convergence = 1000.;
		
		for(nt = 1; (nt <= mfNintT); nt++){
			// Update Mean-field equations, until next scheduled synapse update (or stop early upon convergence)
			it = (mt * mfNintT) + nt;
			
			phi_mu_e = mu_e + (c_ee * w_ee[it-1] * tau_me * nu_e[it-1]) - (c_ei * w_ei * tau_me * nu_i[it-1]);
			phi_mu_i = mu_i + (c_ie * w_ie * tau_mi * nu_e[it-1]) - (c_ii * w_ii * tau_mi * nu_i[it-1]);
			printf("pme: %f, pmi: %f, ", phi_mu_e, phi_mu_i);
			
			x_local = (theta - phi_mu_e)/sigma;
			y_local = (v_r - phi_mu_e)/sigma;
			e_trans = trans(x_local, y_local, tau_me);
			d_nu_e = (-nu_e[it-1] + e_trans) / tau_e;
			nu_e[it] = nu_e[it-1] + (dt * d_nu_e);
			printf("x: %.2f, y: %.2f, trans: %.2f, de: %.2f, ", x_local, y_local, e_trans, d_nu_e);
		
			x_local = (theta - phi_mu_i)/sigma;
			y_local = (v_r - phi_mu_i)/sigma;
			i_trans = trans(x_local, y_local, tau_mi);
			d_nu_i = (-nu_i[it-1] + i_trans) / tau_i;
			nu_i[it] = nu_i[it-1] + (dt * d_nu_i);
			printf("x: %.2f, y: %.2f, trans: %.2f, di: %.2f, ", x_local, y_local, i_trans, d_nu_i);
			
			if(nt < mfNintT){ // On last pass through loop only update MF equations, don't update weight or output to file
				w_ee[it] = w_ee[it-1];
				
				printf("time: %g, nu_e[%d]: %g, nu_i[%d]: %g, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
				fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
				
				convergence = fmax(fabs(d_nu_e), fabs(d_nu_i));
				if(convergence < 0.5){
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
		rho[jt] = updateWeight(rho[jt-1], dwt, nu_e[it], cpre, cpost);
		w_ee[it] = w_min + (w_len * rho[jt]); // -1 due to extra increment on loop above
		/*}
		else{
			w_ee[it] = w_ee[it-1];
		}*/
		
		printf("time: %g, nu_e[%d]: %g, nu_i[%d]: %g, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
		fprintf(output_file, "%f %f %f %g %f\n", (it*dt), nu_e[it], nu_i[it], w_ee[it], rho[jt]);
	}
	
	fprintf(output_file, "\n\n\n\n\n");
	fclose(output_file);
}
    
