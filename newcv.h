/*
 *  newcv.h
 *  MFNetwork
 *
 *  Created by David Higgins on 24/01/2012.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef NEWCV_H
#define NEWCV_H

float nerf(float z); /* function exp(z^2)(erf(z)+1)  */

double trans(double x,double y,double taui, double taurp);     /* transduction function 
									 x=(threshold-mu)/sigma, y=(reset-mu)/sigma */ 

//double function(double x);  /* called from within cv() */
double cv(double x,double y,double taui, double taurp);     /* coefficient of variation  
								  x=(threshold-mu)/sigma, y=(reset-mu)/sigma */ 

#endif /* NEWCV_H */