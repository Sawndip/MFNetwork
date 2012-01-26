#include <math.h>
#include <stdio.h>

#include "newcv.h"

#define c_ee 4000 /* connectivity: each neuron receives X incoming connections */
#define c_ie 4000 
#define c_ei 1000 
#define c_ii 1000
#define w_ie 0.2e-3 /* weights */
#define w_ei 0.6e-3  
#define w_ii 0.6e-3 
#define tau_e 0.01 /* excitatory population time constant (seconds) */
#define tau_i 0.01 //0.01 /* inhibitory population time constant (seconds) */
#define tau_me 0.02 /* excitatory membrane time constant (seconds) */
#define tau_mi 0.02 //0.01 /* inhibitory membrane time constant (seconds) */
#define theta 0.02 /* threshold potential */
#define v_r 0.01 //0. /* reset potential */
#define sigma 0.5e-3 /* noise */
#define tmax 1 //10. /* seconds */
#define dt 0.01
#define NintT ((int)(tmax/dt))
/* tau_rp defined in newcv.c */



FILE *fopen(),*cadistrib;


int main(void) {
	float mu_e, mu_i;
	
	float nu_e[NintT], nu_i[NintT], w_ee[NintT];
	
	int it;
	
	w_ee[0] = 0.2e-3; //0.7;
	nu_e[0] = 5.25; //5.25; /* excitatory population initial rate */
	nu_i[0] = 5.25; //5.25; /* inhibitory population initial rate */
	
	mu_e = 0.01; //0.025; /* external input to excitatory population */
	mu_i = 0.01; //0.025; /* external input to inhibitory population */
	
	it=0;
	printf("Inputs: mu_e: %g, mu_i: %g. Equivalent to: nu_ext_e: %gHz, nu_ext_i: %gHz\n", mu_e, mu_i, (mu_e/(c_ee * w_ee[0] * tau_me)), (mu_i/(c_ie * w_ie * tau_me)));
	printf("time: 0.00, nu_e[0]: %f, nu_i[0]: %f, w_ee[0]: %g\n", nu_e[it], nu_i[it], w_ee[it]);
	
	for(it = 1; it < NintT; it++){
		float phi_mu_e, phi_mu_i;
		float x_local, y_local;
		float e_trans, i_trans;
		float d_nu_e, d_nu_i;
		
		phi_mu_e = mu_e + (c_ee * w_ee[it-1] * tau_me * nu_e[it-1]) - (c_ei * w_ei * tau_mi * nu_i[it-1]);
		phi_mu_i = mu_i + (c_ie * w_ie * tau_me * nu_e[it-1]) - (c_ii * w_ii * tau_mi * nu_i[it-1]);
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
		
		w_ee[it] = w_ee[it-1];
		
		printf("time: %.2f, nu_e[%d]: %.3f, nu_i[%d]: %.3f, w_ee[%d]: %g\n", (it*dt), it, nu_e[it], it, nu_i[it], it, w_ee[it]);
	}
}
    
