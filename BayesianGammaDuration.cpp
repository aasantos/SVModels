//
//  BayesianGammaDuration.cpp
//  October2019
//
//  Created by António Alberto Santos on 25/10/2019.
//  Copyright © 2019 Antonio Santos. All rights reserved.
//

#include <stdio.h>
#include "BayesianGammaDuration.hpp"

double BayesianGammaDuration::lognormalpdf(double x,double m,double s)
{
return -0.5*(x - m)*(x - m)/(s*s);
}

void BayesianGammaDuration::simulatebeta()
{
    double sumd = 0.0;
    double sumexpa = 0.0;
    for(int i=0;i<this->n;i++){
        sumd += dur[i];
        sumexpa += 1.0/sqrt(expf(alpha[i]));
    }
    double mbeta = this->lambda*sumexpa/sumd;
    double sbeta = sqrt(mbeta*mbeta/(this->lambda*sumexpa));
    this->beta = mbeta + sbeta*random->normal();
}


double BayesianGammaDuration::dfa(double x)
{
    double res = 0.0;
    for(int i=0;i<this->n;i++) res+= (log(this->beta) + log(this->dur[i]) - digamma(x/sqrt(expf(alpha[i]))))/sqrt(expf(alpha[i]));
    return res;
}

double BayesianGammaDuration::ddfa(double x)
{
    double res = 0.0;
    for(int i=0;i<this->n;i++) res += -1.0*trigamma(x/sqrt(expf(alpha[i])))/expf(alpha[i]);
    return res;
}

double BayesianGammaDuration::newtownlambda(double l0)
{
    double x0 = l0;
    int niter = 0;
    double g = dfa(x0);
    while(fabs(g) > 0.00001 && niter < 20){
        x0 = x0 - g/ddfa(x0);
        g = dfa(x0);
        niter++;
    }
    return x0;
}

void BayesianGammaDuration::lambdameanvar(double l0,double out[2])
{
    out[0] = newtownlambda(l0);
    out[1] = -1.0/ddfa(out[0]);
}

void BayesianGammaDuration::simulatelambda()
{
    double mv[2];
    lambdameanvar(this->lambda,mv);
    this->lambda = mv[0] + sqrt(mv[1])*random->normal();
}

void BayesianGammaDuration::simulateparameters()
{
    simulatebeta();
    simulatelambda();
}

double BayesianGammaDuration::getlambda()
{
    return this->lambda;
}

double BayesianGammaDuration::getbeta()
{
    return this->beta;
}

