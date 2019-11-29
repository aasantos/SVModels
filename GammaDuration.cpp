//
//  GammaDuration.cpp
//  Nov20
//
//  Created by António Alberto Santos on 20/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "GammaDuration.hpp"

double GammaDuration::lognormalpdf(double x,double m,double s)
{
return -0.5*(x - m)*(x - m)/(s*s);
}

double GammaDuration::dfa(double x)
{
    double res = 0.0;
    for(int i=0;i<this->n;i++){
        double expa = exp(-0.5*this->alpha[i]);
        res += expa*log(this->dur[i]) - expa*digamma(x*expa);
    }
    return res;
}

double GammaDuration::ddfa(double x)
{
    double res = 0.0;
    for(int i=0;i<this->n;i++) {
        double expa = exp(-0.5*this->alpha[i]);
        res += -1.0*expa*expa*trigamma(x*expa);
    }
    return res;
}

double GammaDuration::newtownlambda(double l0)
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

void GammaDuration::lambdameanvar(double l0,double out[2])
{
    out[0] = newtownlambda(l0);
    out[1] = -1.0/ddfa(out[0]);
}

double GammaDuration::loglikdur(double ll)
{
    double lk = 0.0;
    for(int i=0;i<this->n;i++){
        double expa = exp(-0.5*this->alpha[i]);
        lk += ll*expa*log(this->dur[i]) - lgamma(ll*expa);
    }
    return lk;
}


double GammaDuration::loglikdurprior(double ll,double a,double b)
{
    double lkp = this->loglikdur(ll) + (a - 1.0)*log(ll) - b*ll;
    return lkp;
}


double GammaDuration::metroproblambda(double llnew,double ll,double a,double b,double mm,double ss)
{
    double ll1 = loglikdurprior(llnew, a, b);
    double ll2 = loglikdurprior(ll,a, b);
    double ll3 = lognormalpdf(ll,mm, ss);
    double ll4 = lognormalpdf(llnew,mm, ss);
    double temp = ll1 - ll2 + ll3 - ll4;
    if(temp < -15.0) return 0.0;
    if(temp > 0.0) return 1.0;
    return exp(temp);
}
 

void GammaDuration::simulatelambda()
{
    double mv[2];
    lambdameanvar(this->lambda,mv);
    double lambdanew = mv[0] + sqrt(mv[1])*random->normal();
    double mprob = this->metroproblambda(lambdanew, this->lambda,
                                         this->lambdaprior[0], this->lambdaprior[1] ,
                                         mv[0], sqrt(mv[1]));
    if(random->uniform() < mprob) this->lambda = lambdanew;
    if(this->lambda < this->lambdainf) this->lambda = this->lambdainf;
    if(this->lambda > this->lambdasup) this->lambda = this->lambdasup;
}


double GammaDuration::getlambda()
{
    return this->lambda;
}


