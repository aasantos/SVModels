//
//  SVDUR.hpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef SVDUR_hpp
#define SVDUR_hpp

#include <stdio.h>
#include "SV.hpp"
#include "BayesianGammaDuration.hpp"


void svdursimul(double sigmav,double phi1,double sigma1,
                double lambda,double gam,double *y,double *dur,double *alpha,int n,Random *random);

class SVDUR: public SV{

protected:
    double *dur;
    double lambda;
    double gam;

public:
    SVDUR(double *y,int n,double sigmav,double phi1,double sigma1,Random *random,double *dur,double lambda,double gam): SV(y,n,sigmav,phi1,sigma1,random)
    {
        this->dur = dur;
        this->lambda = lambda;
        this->gam = gam;
    }
    
    double dfdur(double dd,double a)
    {
           double tt = exp(-0.5*a);
           return -0.5*this->lambda*log(this->gam*dd)*tt + 0.5*this->lambda*digamma(this->lambda*tt)*tt;
    }
    
    
    double ddfdur(double dd,double a)
    {
        double tt = exp(-0.5*a);
        return 0.25*this->lambda*log(this->gam*dd)*tt - 0.25*this->lambda*digamma(this->lambda*tt)*tt - 0.25*this->lambda*this->lambda*trigamma(this->lambda*tt)*tt*tt;
    }
    
    double loglikdur(double dd,double a)
    {
        double tt = exp(-0.5*a);
        return this->lambda*log(this->gam*dd)*tt - lgamma(this->lambda*tt);
    }
    
    double df(double a,double yy,double dd,double a0,double a1)
    {
        double t1 = a1 - this->phi1*a;
        double t2 = a - this->phi1*a0;
        return  -0.5 + 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a)) + dfdur(dd, a)
        + t1*this->phi1/(this->sigma1*this->sigma1) - t2/(this->sigma1*this->sigma1);
    }
    // calculate the second derivative
    double ddf(double a,double yy,double dd)
    {
        return  -0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
        + ddfdur(dd, a) - (1.0 + this->phi1*this->phi1)/(this->sigma1*this->sigma1);
        
    }
    //
    // calculate the mean
    double meanAlpha(double yy,double dd,double a0,double a1)
    {
        double a = 0.5*(a0 + a1);
        int niter = 0;
        double g = df(a,yy,dd,a0,a1);
        while (fabs(g) > 0.00001 && niter < 20) {
            double h = ddf(a,yy,dd);
            a = a - g/h;
            g = df(a,yy,dd,a0,a1);
            niter++;
        }
        return a;
    }
    //
    //
    double varAlpha(double a,double yy,double dd)
    {
        return -1.0/ddf(a, yy,dd);
    }
    //
    //
    double logpdfsv(double a,double yy,double dd,double a0,double a1)
    {
        double t1 = a1 - this->phi1*a;
        double t2 = a - this->phi1*a0;
        return -0.5*a - 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a))
        + loglikdur(dd,a) - 0.5*t1*t1/(this->sigma1*this->sigma1) - 0.5*t2*t2/(this->sigma1*this->sigma1);
    }
    //
    double logpdfnormal(double a,double m,double s)
    {
        return -0.5*(a - m)*(a - m)/(s*s);
    }
    //
    //
    double metroprob(double anew,double a,double yy,double dd,double a0,double a1,double m,double s)
    {
        double result = 0.0;
        double t1 = logpdfsv(anew, yy,dd, a0, a1);
        double t2 = logpdfsv(a, yy, dd,a0, a1);
        double t3 = logpdfnormal(a, m, s);
        double t4 = logpdfnormal(anew, m, s);
        double tt = t1 - t2 + t3 - t4;
        if(tt > 0.0) tt = 0.01;
        if(tt < -15.0) tt = -15.0;
        result = fminf(1.0, expf(tt));
        return result;
    }
    //
    double singlestepstate(double alphaold,double yy,double dd,double a0,double a1)
     {
         double meana = meanAlpha(yy, dd,a0, a1);
         double vara = varAlpha(meana, yy,dd);
         double alphatemp = meana + sqrt(vara)*random->normal();
         if(random->uniform() < metroprob(alphatemp, alphaold, yy,dd, a0, a1, meana, sqrt(vara))){
             return alphatemp;
         }else{
             return alphaold;
         }
     }
     //
     void simulatestates()
     {
         double a0 = 0.0;
         double a1 = 0.0;
         this->alphafirst = (this->sigma1/sqrt(1.0 - this->phi1*this->phi1))*random->normal();
         for(int i=0;i<n;i++){
             double yy = this->y[i];
             double dd = this->dur[i];
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
             this->alpha[i] = singlestepstate(this->alpha[i], yy,dd, a0, a1);
         }
         alphalast = this->phi1*alpha[n-1] + this->sigma1*random->normal();
     }


     void simulateparameters()
     {
         SV::simulateparameters();
         //
         BayesianGammaDuration *bdur = new BayesianGammaDuration(this->dur,this->alpha,this->n,this->lambda,this->gam,this->random);
         bdur->simulateparameters();
         this->lambda = bdur->getlambda();
         this->gam = bdur->getbeta();
         delete bdur;
     }
     //
     double getlambda()
     {
         return this->lambda;
     }
    //
     double getgam()
     {
         return this->gam;
     }
    
};

#endif /* SVDUR_hpp */
