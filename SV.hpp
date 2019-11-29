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


struct PriorStructSV
{
    int sigmavpriortype;
    int phi1priortype;
    int sigma1priortype;
    double sigmavprior[2];
    double phi1prior[2];
    double sigma1prior[2];
};


//
//
class SV: public DLM{
    
public:
    SV(double *y,int n,double sigmav,double phi1,double sigma1,Random *random): DLM(y,n,sigmav,phi1,sigma1,random)
    {
    }
        
    double df(double a,double yy,double a0,double a1);
    double ddf(double a,double yy);
    double meanAlpha(double yy,double a0,double a1);
    double varAlpha(double a,double yy);
    double logpdfsv(double a,double yy,double a0,double a1);
    double logpdfnormal(double a,double m,double s);
    double metroprob(double anew,double a,double yy,double a0,double a1,double m,double s);
    double singlestepstate(double alphaold,double yy,double a0,double a1);
    void simulatestates();
    void simulatestatesadaptation(int k);
    void simulateparameters();
    void setprior(struct PriorStructSV prior);
    
};

#endif /* SV_hpp */
