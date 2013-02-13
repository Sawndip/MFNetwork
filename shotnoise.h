/*
 *  shotnoise.h
 *  MFNetwork
 *
 *  Created by David Higgins on 01/02/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SHOTNOISE_H
#define SHOTNOISE_H

#include <string.h>

void getAlphas(double rate, double c_pre, double c_post, double *alphas); /* calculates alphas[1]=alphap and alphas[0]=alphad */

double updateWeight(double rho_old, double stepsize, double rate, double c_pre, double c_post); /* returns a new value for rho */

#endif /* SHOTNOISE_H */