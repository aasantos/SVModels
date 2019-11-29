//
//  SVJ2FDUR.hpp
//  November2019
//
//  Created by António Alberto Santos on 07/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef SVJ2FDUR_hpp
#define SVJ2FDUR_hpp

#include <stdio.h>
#include "GammaDuration.hpp"
#include "SVJ2F.hpp"


struct PriorStructSVJ2FDUR
{
    int sigmavpriortype;
    int phi1priortype;
    int sigma1priortype;
    int phi2priortype;
    int sigma2priortype;
    int mujpriortype;
    int sigmajpriortype;
    double sigmavprior[2];
    double phi1prior[2];
    double sigma1prior[2];
    double phi2prior[2];
    double sigma2prior[2];
    double mujprior[2];
    double sigmajprior[2];
    double thetaprior[2];
    double lambdaprior[2];
};



class SVJ2FDUR: public SVJ2F{
//
protected:
    double *dur;
    double lambda;
    double lambdaprior[2];

//
public:
    SVJ2FDUR(double *y,int n,double sigmav,
             double phi1,double sigma1,
             Random *random,double phi2,double sigma2,
             double muj,double sigmaj,double theta,
             double *dur,double lambda):
    SVJ2F(y,n,sigmav,phi1,sigma1,random,phi2,sigma2,muj,sigmaj,theta)
    {
        this->dur = dur;
        this->lambda = lambda;
        this->lambdaprior[0] = Stats(this->dur, this->n).mean();
        this->lambdaprior[1] = 1.0;
    }
    
    void firstsecondderiv(double yy,double dd,double a,double b, double a0,double a1,
                          double b0,double b1,double sz0,double sz1,double fd[2],double sdinv[4]);
    
    void newton(double yy,double dd, double a,double b,double a0,double a1,double b0,
                double b1,double sz0,double sz1,double sol[2],double sdinv[4]);
    
    double logpdfbivariatenormal(double x[2],double mm[2],double vcv[4]);
    //
    double logpdfsv2fdur(double yy,double dd,double a0,double a,double a1,
                         double b0,double b,double b1,double sz0,double sz1);
    //
    double metroprob(double anew,double bnew, double yy,double dd,
                     double a0,double a,double a1,
                     double b0,double b,double b1,double sz0,double sz1,
                     double mm[2],double vcv[4]);
    //
    void singlestatesimulate(double yy,double dd,double a0,double a,double a1,
                             double b0,double b,double b1,double sz0,double sz1,
                             double abnew[2]);

    void simulatestates();
    //
    void singlestatesimulateadaptation(double yy,double dd,double a0,double a,double a1,
                             double b0,double b,double b1,double sz0,double sz1,
                             double abnew[2]);

    void simulatestatesadaptation(int k);
    
    //
    void simulateparameters();
    
    void setprior(struct PriorStructSVJ2FDUR prior);
    
    double getlambda();
    
};

#endif /* SVJ2FDUR_hpp */
