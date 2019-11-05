//
//  Utils.hpp
//  October2019
//
//  Created by António Alberto Santos on 15/10/2019.
//  Copyright © 2019 Antonio Santos. All rights reserved.
//

#ifndef Utils_hpp
#define Utils_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <gsl/gsl_sf.h>


float* readArrayFloat(const char* file,int *n);
double* readArrayDouble(const char* file,int *n);
void writeArrayDouble(double *x,const char* file,int n);
//
double digamma(double x);
double trigamma(double x);

double norm(double x[2]);

double logpdfbivariatenormal(double x[2],double mm[2],double vcv[4]);

#endif /* Utils_hpp */
