//
//  SV2FDUR.hpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//
#ifndef SV2FDUR_hpp
#define SV2FDUR_hpp

#include <stdio.h>
#include "SV2F.hpp"
#include "GammaDuration.hpp"

struct PriorStructSV2FDUR{
    int sigmavpriortype;
    int phi1priortype;
    int sigma1priortype;
    double sigmavprior[2];
    double phi1prior[2];
    double sigma1prior[2];
    int phi2priortype;
    int sigma2priortype;
    double phi2prior[2];
    double sigma2prior[2];
    double lambdaprior[2];
};

class SV2FDUR: public SV2F{
    
protected:
    double *dur;
    double lambda;
    double lambdaprior[2];


public:
    SV2FDUR(double *y,int n,double sigmav,double phi1,double sigma1,Random *random,double phi2,double sigma2,double *dur,double lambda):SV2F(y,n,sigmav,phi1,sigma1,random,phi2,sigma2)
    {
        this->dur = dur;
        this->lambda = lambda;
        this->lambdaprior[0] = Stats(this->dur, this->n).mean();
        this->lambdaprior[1] = 1.0;
    }
    //
    void firstsecondderiv(double yy,double dd,double a,double b,double a0,
                          double a1,double b0,double b1,double fd[2],double sdinv[4]);
    
    void newton(double yy,double dd,double a,double b, double a0,
                double a1,double b0,double b1, double sol[2],double sdinv[4]);
    
    double logpdfbivariatenormal(double x[2],double mm[2],double vcv[4]);
    //
    double logpdfsv2fdur(double yy,double dd,double a0,double a,double a1,
                         double b0,double b,double b1);
    //
    double metroprob(double anew,double bnew, double yy,double dd,
                     double a0,double a,double a1,
                     double b0,double b,double b1,
                     double mm[2],double vcv[4]);
    //
    void singlestatesimulate(double yy,double dd,double a0,double a,double a1,
                             double b0,double b,double b1,
                             double abnew[2]);
    //
    void singlestatesimulateadaptation(double yy,double dd,double a0,double a,double a1,
                                  double b0,double b,double b1, double abnew[2]);
    //
    void simulatestatesadaptation(int k);
    //
    void simulatestates();
    void simulateparameters();
    //
    void setprior(struct PriorStructSV2FDUR prior);
    //
    double getlambda();
    
};

#endif /* SV2FDUR_hpp */
