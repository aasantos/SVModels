//
//  SV.cpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "SV.hpp"


double SV::df(double a,double yy,double a0,double a1)
{
    double result = 0.0;
    double t1 = a1 - this->phi1*a;
    double t2 = a - this->phi1*a0;
    result = -0.5 + 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
    + t1*this->phi1/(this->sigma1*this->sigma1) - t2/(this->sigma1*this->sigma1);
    return result;
}

double SV::ddf(double a,double yy)
{
    double result = 0.0;
    result = -0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
    - (1.0 + this->phi1*this->phi1)/(this->sigma1*this->sigma1);
    return result;
}

double SV::meanAlpha(double yy,double a0,double a1)
{
    double a = 0.5*(a0 + a1);
    int niter = 0;
    double g = df(a,yy,a0,a1);
    while (fabs(g) > 0.00001 && niter < 20) {
        double h = ddf(a,yy);
        a = a - g/h;
        g = df(a,yy,a0,a1);
        niter++;
    }
    return a;
}


double SV::varAlpha(double a,double yy)
{
    return -1.0/ddf(a, yy);
}


double SV::logpdfsv(double a,double yy,double a0,double a1)
{
    double t1 = a1 - this->phi1*a;
    double t2 = a - this->phi1*a0;
    return -0.5*a - 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
    - 0.5*t1*t1/(this->sigma1*this->sigma1) - 0.5*t2*t2/(this->sigma1*this->sigma1);
}

double SV::logpdfnormal(double a,double m,double s)
{
    return -0.5*(a - m)*(a - m)/(s*s);
}


double SV::metroprob(double anew,double a,double yy,double a0,double a1,double m,double s)
{
    double result = 0.0;
    double t1 = logpdfsv(anew, yy, a0, a1);
    double t2 = logpdfsv(a, yy, a0, a1);
    double t3 = logpdfnormal(a, m, s);
    double t4 = logpdfnormal(anew, m, s);
    double tt = t1 - t2 + t3 - t4;
    if(tt > 0.0) tt = 0.01;
    if(tt < -15.0) tt = -15.0;
    result = fminf(1.0, expf(tt));
    return result;
}

double SV::singlestepstate(double alphaold,double yy,double a0,double a1)
{
    double meana = meanAlpha(yy, a0, a1);
    double vara = varAlpha(meana, yy);
    double alphatemp = meana + sqrt(vara)*random->normal();
    if(random->uniform() < metroprob(alphatemp, alphaold, yy, a0, a1, meana, sqrt(vara))){
        return alphatemp;
    }else{
        return alphaold;
    }
}

void SV::simulatestates()
{
    double a0 = 0.0;
    double a1 = 0.0;
    this->alphafirst = (this->sigma1/sqrt(1.0 - this->phi1*this->phi1))*random->normal();
    for(int i=0;i<n;i++){
        double yy = this->y[i];
        if(i==0){
            a0 = this->alphafirst;
            a1 = this->alpha[i+1];
        }else if(i==(n-1)){
            a0 = this->alpha[i-1];
            a1 = this->alphalast;
        }else{
            a0 = this->alpha[i-1];
            a1 = this->alpha[i+1];
        }
        this->alpha[i] = singlestepstate(this->alpha[i], yy, a0, a1);
    }
    alphalast = this->phi1*alpha[n-1] + this->sigma1*random->normal();
}


void SV::simulatestatesadaptation(int k)
{
    for(int i=0;i<k;i++) this->simulatestates();
}



void SV::simulateparameters()
{
    AR1Model *ar1 = new AR1Model(this->alpha,this->n,this->phi1,this->sigma1,this->random);
    //
    ar1->setupdatemu(this->updatephi1);
    ar1->setmudiffuse(this->phi1diffuse);
    ar1->setmupriortype(this->phi1priortype);
    ar1->setmuprior(this->phi1prior);
    ar1->setmubetaprior(this->phi1betaprior);
    //
    ar1->setupdatesigma(this->updatesigma1);
    ar1->setsigmadiffuse(this->sigma1diffuse);
    ar1->setsigmapriortype(this->sigma1priortype);
    ar1->setsigmaprior(this->sigma1prior);
    ar1->setsigmalognormalprior(this->sigma1lognormalprior);
    ar1->setsigmainvgaussianprior(this->sigma1invgaussianprior);
    //
    ar1->simulateparameters();
    this->phi1 = ar1->getmu();
    this->sigma1 = ar1->getsigma();
    //
    delete ar1;
    //
    if(this->updatesigmav){
        double *err = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++) err[i] = this->y[i]/exp(0.5*this->alpha[i]);
        NormalBayesian *nb = new NormalBayesian(err,this->n,0.0,this->sigmav,this->random);
        nb->setupdatemu(false);
        nb->setupdatesigma(true);
        nb->setsigmadiffuse(this->sigmavdiffuse);
        nb->setsigmapriortype(this->sigmavpriortype);
        nb->setsigmaprior(this->sigmavprior);
        nb->setsigmalognormalprior(this->sigmavlognormalprior);
        nb->setsigmainvgaussianprior(this->sigmavinvgaussianprior);
        nb->simulateparameters();
        this->sigmav = nb->getsigma();
        delete nb;
        free(err);
    }
    
}

void SV::setprior(struct PriorStructSV prior)
{
    this->sigmavpriortype = prior.sigmavpriortype;
    this->phi1priortype = prior.phi1priortype;
    this->sigma1priortype = prior.sigma1priortype;
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
}


