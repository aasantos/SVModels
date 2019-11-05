//
//  SVDUR.cpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "SVDUR.hpp"

void svdursimul(double sigmav,double phi1,double sigma1,double lambda,double gam,double *y,double *dur,double *alpha,int n,Random *random)
{
    double a = 0.0;
    for(int i=0;i<100;i++) a = phi1*a + sigma1*random->normal();
    for(int i=0;i<n;i++){
        a = phi1*a + sigma1*random->normal();
        alpha[i] = a;
        y[i] = sigmav*exp(0.5*a)*random->normal();
        dur[i] = random->gamma(lambda*exp(-0.5*a), gam);
    }
}
