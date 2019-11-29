//
//  func.hpp
//  Nov20
//
//  Created by António Alberto Santos on 20/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef func_hpp
#define func_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <gsl/gsl_sf.h>

double digamma(double x);
double trigamma(double x);
double norm(double x[2]);


//
float* readArrayFloat(const char* file,int *n);
double* readArrayDouble(const char* file,int *n);
void writeArrayDouble(double *x,const char* file,int n);

#endif /* func_hpp */
