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
#include "GammaDuration.hpp"


struct PriorStructSVDUR
{
    int sigmavpriortype;
    int phi1priortype;
    int sigma1priortype;
    double sigmavprior[2];
    double phi1prior[2];
    double sigma1prior[2];
    double lambdaprior[2];
};



class SVDUR: public SV{

protected:
    double *dur;
    double lambda;
    double lambdaprior[2];

public:
    SVDUR(double *y,int n,double sigmav,double phi1,double sigma1,Random *random,double *dur,double lambda): SV(y,n,sigmav,phi1,sigma1,random)
    {
        this->dur = dur;
        this->lambda = lambda;
        this->lambdaprior[0] = Stats(this->dur,this->n).mean();
        this->lambdaprior[1] = 1.0;
    }
    
    double df(double a,double yy,double dd,double a0,double a1);
    double ddf(double a,double yy,double dd);
    double meanAlpha(double yy,double dd,double a0,double a1);
    double varAlpha(double a,double yy,double dd);
    double logpdfsvdur(double a,double yy,double dd,double a0,double a1);
    double logpdfnormal(double a,double m,double s);
    double metroprob(double anew,double a,double yy,double dd,double a0,double a1,double m,double s);
    double singlestepstate(double alphaold,double yy,double dd,double a0,double a1);
    void simulatestates();
    void simulatestatesadaptation(int k);
    void simulateparameters();
    void setprior(struct PriorStructSVDUR prior);
    double getlambda();
    
};

#endif /* SVDUR_hpp */
