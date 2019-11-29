//
//  SVJ.cpp
//  November2019
//
//  Created by António Alberto Santos on 08/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include <stdio.h>
#include "SVJ.hpp"


double SVJ::df(double a,double yy,double a0,double a1,double sz0,double sz1)
{
    double result = 0.0;
    double t1 = a1 - this->phi1*a - sz1;
    double t2 = a - this->phi1*a0 - sz0;
    result =  - 0.5 + 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
    + t1*this->phi1/(this->sigma1*this->sigma1) - t2/(this->sigma1*this->sigma1);
    return result;
}

//
//
double SVJ::ddf(double a,double yy)
{
    double result = 0.0;
    result =  - 0.5*yy*yy/(sigmav*sigmav*exp(a))
    - (1.0 + this->phi1*this->phi1)/(this->sigma1*this->sigma1);
    return result;
}

//
//
double SVJ::meanAlpha(double yy,double a0,double a1,double sz0,double sz1)
{
    double a = 0.5*(a0 + a1);
    int niter = 0;
    double g = df(a,yy,a0,a1,sz0,sz1);
    while (fabs(g) > 0.00001 && niter < 20) {
        double h = ddf(a,yy);
        a = a - g/h;
        g = df(a,yy,a0,a1,sz0,sz1);
        niter++;
    }
    return a;
}
//
//
double SVJ::varAlpha(double a,double yy)
{
    return -1.0/ddf(a,yy);
}
//
//
double SVJ::logpdfsv(double a,double yy,double a0,double a1,double sz0,double sz1)
{
    double t1 = a1 - this->phi1*a - sz1;
    double t2 = a - this->phi1*a0 - sz0;
    return  - 0.5*a - 0.5*yy*yy/(sigmav*sigmav*exp(a))
    - 0.5*t1*t1/(this->sigma1*this->sigma1) - 0.5*t2*t2/(this->sigma1*this->sigma1);
}
//
//
double SVJ::logpdfnormal(double a,double m,double s)
{
    return -0.5*(a - m)*(a - m)/(s*s);
}
//
//
double SVJ::metroprob(double anew,double a,double yy,double a0,double a1,double sz0,double sz1,double m,double s)
{
    double t1 = logpdfsv(anew,yy,a0,a1,sz0,sz1);
    double t2 = logpdfsv(a,yy,a0,a1,sz0,sz1);
    double t3 = logpdfnormal(a,m,s);
    double t4 = logpdfnormal(anew,m,s);
    double ee = t1 - t2 + t3 - t4;
    if(ee > 0.0) ee = 0.01;
    if(ee < -15.0) ee = -15.0;
    return fminf(1.0,expf(ee));
}
//
double SVJ::singlestepstate(double alphaold,double yy,double a0,double a1,double sz0,double sz1)
{
    double meana = meanAlpha(yy,a0,a1,sz0,sz1);
    double vara = varAlpha(meana,yy);
    double alphatemp = meana + sqrt(vara)*random->normal();
    if(random->uniform() < metroprob(alphatemp, alphaold, yy, a0, a1, sz0, sz1, meana, sqrt(vara))){
        return alphatemp;
    }else{
        return alphaold;
    }
}
//
void SVJ::simulatestates()
{
    double a0 = 0.0;
    double a1 = 0.0;
    double sz0 = 0.0;
    double sz1 = 0.0;
    this->alphafirst = (this->sigma1/sqrt(1.0 - this->phi1*this->phi1))*random->normal();
    for(int i=0;i<n;i++){
        double yy = this->y[i];
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
        this->alpha[i] = singlestepstate(this->alpha[i],yy,a0, a1,sz0,sz1);
    }
    this->alphalast = this->phi1*alpha[n-1] + this->sigma1*random->normal();
}

void SVJ::simulatestatesadaptation(int k)
{
    //SV::statesadaptation(k);
    for(int i=0;i<k;i++) this->simulatestates();
}

//
void SVJ::simulateparameters()
{
    AR1JModel *ar1j = new AR1JModel(this->alpha,this->n,this->phi1,this->sigma1,this->muj,this->sigmaj,this->theta,this->random);
    //
    ar1j->setupdatemu(this->updatephi1);
    ar1j->setmudiffuse(this->phi1diffuse);
    ar1j->setmupriortype(this->phi1priortype);
    ar1j->setmuprior(this->phi1prior);
    ar1j->setmubetaprior(this->phi1betaprior);
    //
    ar1j->setupdatesigma(this->updatesigma1);
    ar1j->setsigmadiffuse(this->sigma1diffuse);
    ar1j->setsigmapriortype(this->sigma1priortype);
    ar1j->setsigmaprior(this->sigma1prior);
    ar1j->setsigmalognormalprior(this->sigma1lognormalprior);
    ar1j->setsigmainvgaussianprior(this->sigma1invgaussianprior);
    //
    ar1j->setupdatemuj(this->updatemuj);
    ar1j->setmujdiffuse(this->mujdiffuse);
    ar1j->setmujpriortype(this->mujpriortype);
    ar1j->setmujprior(this->mujprior);
    ar1j->setmujbetaprior(this->mujbetaprior);
    //
    ar1j->setupdatesigmaj(this->updatesigmaj);
    ar1j->setsigmajdiffuse(this->sigmajdiffuse);
    ar1j->setsigmajpriortype(this->sigmajpriortype);
    ar1j->setsigmajprior(this->sigmajprior);
    ar1j->setsigmajlognormalprior(this->sigmajlognormalprior);
    ar1j->setsigmajinvgaussianprior(this->sigmajinvgaussianprior);
    //
    ar1j->setupdatetheta(this->updatetheta);
    ar1j->setthetaprior(this->thetaprior);
    //
    ar1j->simulateparameters();
    //
    this->phi1   =  ar1j->getmu();
    this->sigma1 =  ar1j->getsigma();
    this->muj    =  ar1j->getmuj();
    this->sigmaj =  ar1j->getsigmaj();
    this->theta  =  ar1j->gettheta();
    //
    for(int i=0;i<(this->n-1);i++){
        this->S[i] = ar1j->getS(i);
        this->Z[i] = ar1j->getZ(i);
    }
    this->S[this->n - 1] = 0;
    this->Z[this->n - 1] = 0.0;
    //
    delete ar1j;
    //
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
    //
}
//
//
double SVJ::getmuj()
{
    return this->muj;
}
double SVJ::getsigmaj()
{
    return this->sigmaj;
}
double SVJ::gettheta()
{
    return this->theta;
}
//
void SVJ::setupdatemuj(bool flag)
{
    this->updatemuj = flag;
}
void SVJ::setupdatesigmaj(bool flag)
{
    this->updatesigmaj = flag;
}
void SVJ::setupdatetheta(bool flag)
{
    this->updatetheta = flag;
}
//
void SVJ::setmujdiffuse(bool flag)
{
    this->mujdiffuse = flag;
}
void SVJ::setsigmajdiffuse(bool flag)
{
    this->sigmajdiffuse = flag;
}
void SVJ::setthetadiffuse(bool flag)
{
    this->thetadiffuse = flag;
}

void SVJ::setpriror(struct PriorStructSVJ prior)
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
    
    this->thetaprior[0] = prior.thetaprior[0];
    this->thetaprior[1] = prior.thetaprior[1];
}
