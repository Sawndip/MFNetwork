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

void getAlphas(float rate, float c_pre, float c_post, float *alphas); /* calculates alphas[1]=alphap and alphas[0]=alphad */

float updateWeight(float rho_old, float stepsize, float rate, float c_pre, float c_post); /* returns a new value for rho */

#endif /* SHOTNOISE_H */