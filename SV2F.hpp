//
//  SV2F.hpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef SV2F_hpp
#define SV2F_hpp

#include <stdio.h>
#include "SV.hpp"


void sv2fsimul(double sigmav,double phi1,double sigma1,double phi2,double sigma2,double *y,double *alpha,double *beta,int n,Random *random);

struct PriorStruct2f : PriorStruct
{
    int phi2priortype;
    int sigma2priortype;
    double phi2prior[2];
    double sigma2prior[2];
};


struct InitsStruct2f : InitsStruct
{
    double phi2;
    double sigma2;
    double *beta;
};



void sv2fsimul(double sigmav,double phi1,double sigma1,
               double phi2,double sigma2,double *y,
               double *alpha,double *beta,int n,Random *random);


void svfirstsecondderiv(double sigmav,double phi1,double sigma1,
                        double phi2,double sigma2,
                        double y,double a,double b,
                        double a0,double a1,double b0,double b1,
                        double fd[2],double sdinv[4]);

void svnewton(double sigmav,double phi1,double sigma1,
              double phi2,double sigma2,double y,double a,
              double b,double a0,double a1,double b0,
              double b1,double sol[2],double sdinv[4]);


class SV2F: public SV{
protected:
    double phi2;
    double sigma2;
    double *beta;
    //
    //
    bool updatephi2;
    bool updatesigma2;
    //
    bool phi2diffuse;
    bool sigma2diffuse;
    //
    int phi2priortype; // 0-normal 1-beta
    double phi2prior[2];
    double phi2betaprior[2];
    //
    int sigma2priortype; //0-invgamma 1-lognormal 2-invgaussian
    double sigma2prior[2];
    double sigma2lognormalprior[2];
    double sigma2invgaussianprior[2];
    //
    ////
    double betafirst;
    double betalast;
    
public:
    SV2F(double *y,int n,double sigmav,double phi1,double sigma1,Random *random,double phi2,double sigma2):SV(y,n,sigmav,phi1,sigma1,random)
    {
        this->phi2 = phi2;
        this->sigma2 = sigma2;
        this->beta = (double*)calloc(n, sizeof(double));
        //
        this->updatephi2 = true;
        this->updatesigma2 = true;
        //
        this->phi2diffuse = false;
        this->sigma2diffuse = false;
        //
        this->phi2priortype = 1; // 0-normal 1-beta
        this->phi2prior[0] = 0.9;
        this->phi2prior[1] = 0.5;
        this->phi2betaprior[0] = 90.0;
        this->phi2betaprior[1] = 10.0;
        //
        this->sigma2priortype = 1; //0-invgamma 1-lognormal 2-invgaussian
        this->sigma2prior[0] = 2.5;
        this->sigma2prior[1] = 0.025;
        this->sigma2lognormalprior[0] = -2.0;
        this->sigma2lognormalprior[1] = 0.5;
        this->sigma2invgaussianprior[0]  = 0.5;
        this->sigma2invgaussianprior[1] = 0.5;
        
        //
        this->betafirst = 0.0;
        this->betalast = 0.0;
    }
    
    ~SV2F()
    {
        if(beta){
            free(beta);
        }
    }
    
    
     double logpdfbivariatenormal(double x[2],double mm[2],double vcv[4]);
     //
     double logpdfsv2f(double yy,double a0,double a,double a1,
                       double b0,double b,double b1);
     //
     double metroprob(double anew,double bnew, double yy,
                      double a0,double a,double a1,
                      double b0,double b,double b1,
                      double mm[2],double vcv[4]);
     //
     void singlestatesimulate(double y,double a0,double a,double a1,
                              double b0,double b,double b1,
                              double abnew[2]);
     //
     void simulatestates();
      //
     void singlestatesimulateAdaptation(double y,double a0,double a,double a1,
                              double b0,double b,double b1,
                              double abnew[2]);
     //
     void simulatestatesAdaptation(int k);
     void simulateparameters();
     //
     double getalpha(int k);
     double getbeta(int k);
    //
    double getphi2();
    double getsigma2();
    
     void setalpha(double *a);
     void setbeta(double *b);
    
     //
     
    
    void setprior(struct PriorStruct2f prior);
    void setinits(struct InitsStruct2f inits);
};

#endif /* SV2F_hpp */
