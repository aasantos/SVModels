//
//  SV.hpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef SV_hpp
#define SV_hpp

#include <stdio.h>
#include "DLM.hpp"

void svsimul(double sigmav,double phi1,double sigma1,double *y,double *alpha,int n,Random *random);
//
//
class SV: public DLM{
    
public:
    SV(double *y,int n,double sigmav,double phi1,double sigma1,Random *random): DLM(y,n,sigmav,phi1,sigma1,random)
    {
        
    }
    
    
    double df(double a,double yy,double a0,double a1)
    {
        double result = 0.0;
        double t1 = a1 - this->phi1*a;
        double t2 = a - this->phi1*a0;
        result = -0.5 + 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
        + t1*this->phi1/(this->sigma1*this->sigma1) - t2/(this->sigma1*this->sigma1);
        return result;
    }
    
    double ddf(double a,double yy)
    {
        double result = 0.0;
        result = -0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
        - (1.0 + this->phi1*this->phi1)/(this->sigma1*this->sigma1);
        return result;
    }
    
    double meanAlpha(double yy,double a0,double a1)
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

    
    double varAlpha(double a,double yy)
    {
        return -1.0/ddf(a, yy);
    }

    
    double logpdfsv(double a,double yy,double a0,double a1)
    {
        double t1 = a1 - this->phi1*a;
        double t2 = a - this->phi1*a0;
        return -0.5*a - 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
        - 0.5*t1*t1/(this->sigma1*this->sigma1) - 0.5*t2*t2/(this->sigma1*this->sigma1);
    }
    
    double logpdfnormal(double a,double m,double s)
    {
        return -0.5*(a - m)*(a - m)/(s*s);
    }
    
    
    double metroprob(double anew,double a,double yy,double a0,double a1,double m,double s)
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

    double singlestepstate(double alphaold,double yy,double a0,double a1)
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

    void simulatestates()
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

    
    
    void simulateparameters()
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

    
    
    
    
    

    
};

#endif /* SV_hpp */
