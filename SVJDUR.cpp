//
//  SVJDUR.cpp
//  November2019
//
//  Created by António Alberto Santos on 08/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include <stdio.h>
#include "SVJDUR.hpp"

//
double  SVJDUR::df(double a,double yy,double dd,double a0,double a1,double sz0,double sz1)
{
    double result = 0.0;
    double x1 = 0.5*this->lambda*exp(-0.5*a);
    double x2 = -1.0*x1*log(dd) + x1*digamma(2.0*x1);
    double t1 = a1 - this->phi1*a - sz1;
    double t2 = a - this->phi1*a0 - sz0;
    result = x2 - 0.5 + 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a)) 
    + t1*this->phi1/(this->sigma1*this->sigma1) - t2/(this->sigma1*this->sigma1);
    return result;
}
//
//
double  SVJDUR::ddf(double a,double yy,double dd)
{
    double result = 0.0;
    double x1 =  this->lambda*exp(-0.5*a);
    double x2 = 0.25*log(dd) - 0.25*x1*digamma(x1);
    double x3 = x2 - 0.25*x1*x1*trigamma(x1);
    result = x3 - 0.5*yy*yy/(sigmav*sigmav*exp(a)) - (1.0 + this->phi1*this->phi1)/(this->sigma1*this->sigma1);
    return result;
}
//
//
double  SVJDUR::meanAlpha(double yy,double dd,double a0,double a1,double sz0,double sz1)
{
    double a = 0.5*(a0 + a1);
    int niter = 0;
    double g = df(a,yy,dd,a0,a1,sz0,sz1);
    while (fabs(g) > 0.00001 && niter < 20) {
        double h = ddf(a,yy,dd);
        a = a - g/h;
        g = df(a,yy,dd,a0,a1,sz0,sz1);
        niter++;
    }
    return a;
}
//
//
double  SVJDUR::varAlpha(double a,double yy,double dd)
{
    return -1.0/ddf(a,yy,dd);
}
//
//
double  SVJDUR::logpdfsv(double a,double yy,double dd,double a0,double a1,double sz0,double sz1)
{
    double t1 = a1 - this->phi1*a - sz1;
    double t2 = a - this->phi1*a0 - sz0;
    double x1 = this->lambda*exp(-0.5*a);
    double x2 = x1*log(dd) - lgamma(x1);
    return  x2 - 0.5*a - 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
    - 0.5*t1*t1/(this->sigma1*this->sigma1) - 0.5*t2*t2/(this->sigma1*this->sigma1);
}
//
double  SVJDUR::logpdfnormal(double a,double m,double s)
{
    return -0.5*(a - m)*(a - m)/(s*s);
}
//
double  SVJDUR::metroprob(double anew,double a,double yy,double dd,double a0,double a1,double sz0,double sz1,double m,double s)
{
    double t1 = logpdfsv(anew,yy,dd,a0,a1,sz0,sz1);
    double t2 = logpdfsv(a,yy,dd,a0,a1,sz0,sz1);
    double t3 = logpdfnormal(a,m,s);
    double t4 = logpdfnormal(anew,m,s);
    double ee = t1 - t2 + t3 - t4;
    if(ee > 0.0) ee = 0.01;
    if(ee < -15.0) ee = -15.0;
    return fminf(1.0,expf(ee));
}
//
double  SVJDUR::singlestepstate(double alphaold,double yy,double dd,double a0,double a1,double sz0,double sz1)
{
    double meana = meanAlpha(yy,dd,a0,a1,sz0,sz1);
    double vara = varAlpha(meana,yy,dd);
    double alphatemp = meana + sqrt(vara)*random->normal();
    if(random->uniform() < metroprob(alphatemp, alphaold, yy,dd, a0, a1, sz0, sz1, meana, sqrt(vara))){
        return alphatemp;
    }else{
        return alphaold;
    }
}
//
void  SVJDUR::simulatestates()
{
    double a0 = 0.0;
    double a1 = 0.0;
    double sz0 = 0.0;
    double sz1 = 0.0;
    alphafirst = (this->sigma1/sqrt(1.0 - this->phi1*this->phi1))*random->normal();
    for(int i=0;i<n;i++){
        double yy = this->y[i];
        double dd = this->dur[i];
        if(i==0){
            a0 = alphafirst;
            a1 = alpha[i+1];
            sz0 = 0.0;
            sz1 = (double)this->S[i]*Z[i];
        }else if(i==(n-1)){
            a0 = alpha[i-1];
            a1 = alphalast;
            sz0 = (double)this->S[i-1]*Z[i-1];
            sz1 = (double)this->S[i]*Z[i];
        }else{
            a0 = alpha[i-1];
            a1 = alpha[i+1];
            sz0 = (double)this->S[i-1]*Z[i-1];
            sz1 = (double)this->S[i]*Z[i];
        }
        this->alpha[i] = singlestepstate(this->alpha[i],yy,dd,a0,a1,sz0,sz1);
    }
    alphalast = this->phi1*alpha[n-1] + this->sigma1*random->normal();
}
//
void SVJDUR::simulatestatesadaptation(int k)
{
    for(int i=0;i<k;i++) this->simulatestates();
}

//
void  SVJDUR::simulateparameters()
{
    SVJ::simulateparameters();
    //
    /*
    GammaDuration *gammadur = new GammaDuration(this->dur,this->alpha,this->n,this->lambda,this->random);
    gammadur->simulatelambda();
    this->lambda = gammadur->getlambda();
    delete gammadur;
    */
}
//
double  SVJDUR::getlambda()
{
    return this->lambda;
}
//

void SVJDUR::setprior(struct PriorStructSVJDUR prior)
{
    this->sigmavpriortype = prior.sigmavpriortype;
    this->phi1priortype = prior.phi1priortype;
    this->sigma1priortype = prior.sigma1priortype;
    this->mujpriortype = prior.mujpriortype;
    this->sigmajpriortype = prior.sigmajpriortype;
    //
    if(this->sigmavpriortype == 0){
        this->sigmavprior[0] = prior.sigmavprior[0];
        this->sigmavprior[1] = prior.sigmavprior[1];
    }else if(this->sigmavpriortype == 1){
        this->sigmavlognormalprior[0] = prior.sigmavprior[0];
        this->sigmavlognormalprior[1] = prior.sigmavprior[1];
    }else{
        this->sigmavinvgaussianprior[0] = prior.sigmavprior[0];
        this->sigmavinvgaussianprior[1] = prior.sigmavprior[1];
    }
    //
    if(this->phi1priortype == 0){
        this->phi1prior[0] = prior.phi1prior[0];
        this->phi1prior[1] = prior.phi1prior[1];
    }else{
        this->phi1betaprior[0] = prior.phi1prior[0];
        this->phi1betaprior[1] = prior.phi1prior[1];
    }
    //
    if(this->sigma1priortype == 0){
        this->sigma1prior[0] = prior.sigma1prior[0];
        this->sigma1prior[1] = prior.sigma1prior[1];
    }else if(this->sigma1priortype == 1){
        this->sigma1lognormalprior[0] = prior.sigma1prior[0];
        this->sigma1lognormalprior[1] = prior.sigma1prior[1];
    }else{
        this->sigma1invgaussianprior[0] = prior.sigma1prior[0];
        this->sigma1invgaussianprior[1] = prior.sigma1prior[1];
    }
    //
    if(this->mujpriortype == 0){
        this->mujprior[0] = prior.mujprior[0];
        this->mujprior[1] = prior.mujprior[1];
    }else{
        this->mujbetaprior[0] = prior.mujprior[0];
        this->mujbetaprior[1] = prior.mujprior[1];
    }
    //
    if(this->sigmajpriortype == 0){
        this->sigmajprior[0] = prior.sigmajprior[0];
        this->sigmajprior[1] = prior.sigmajprior[1];
    }else if(this->sigmajpriortype == 1){
        this->sigmajlognormalprior[0] = prior.sigmajprior[0];
        this->sigmajlognormalprior[1] = prior.sigmajprior[1];
    }else{
        this->sigmajinvgaussianprior[0] = prior.sigmajprior[0];
        this->sigmajinvgaussianprior[1] = prior.sigmajprior[1];
    }
    //
    this->thetaprior[0] = prior.thetaprior[0];
    this->thetaprior[1] = prior.thetaprior[1];
    //
    this->lambdaprior[0] = prior.lambdaprior[0];
    this->lambdaprior[1] = prior.lambdaprior[1];
    //
}
