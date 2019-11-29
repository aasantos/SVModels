#ifndef SVJDUR_hpp
#define SVJDUR_hpp

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "GammaDuration.hpp"
#include "SVJ.hpp"


struct PriorStructSVJDUR
{
    int sigmavpriortype;
    int phi1priortype;
    int sigma1priortype;
    int mujpriortype;
    int sigmajpriortype;
    double sigmavprior[2];
    double phi1prior[2];
    double sigma1prior[2];
    double mujprior[2];
    double sigmajprior[2];
    double thetaprior[2];
    double lambdaprior[2];
};

class SVJDUR: public SVJ{
    
protected:
    double *dur;
    double lambda;
    double lambdaprior[2];

    
public:
    SVJDUR(double *x,int n,double sigmav,double phi1,double sigma1,Random *random,
        double muj,double sigmaj,double theta,
        double *dur,double lambda):
    SVJ(x,n,sigmav,phi1,sigma1,random,muj,sigmaj,theta)
    {
        this->dur = dur;
        this->lambda = lambda;
        this->lambdaprior[0] = Stats(this->dur, this->n).mean();
        this->lambdaprior[1] = 1.0;
    }
    //
    double dfdur(double dd,double a);
    double ddfdur(double dd,double a);
    double loglikdur(double dd,double a);
    double df(double a,double yy,double dd,double a0,double a1,double sz0,double sz1);
    double ddf(double a,double yy,double dd);
    double meanAlpha(double yy,double dd,double a0,double a1,double sz0,double sz1);
    double varAlpha(double a,double yy,double dd);
    double logpdfsv(double a,double yy,double dd,double a0,double a1,double sz0,double sz1);
    double logpdfnormal(double a,double m,double s);
    double metroprob(double anew,double a,double yy,double dd,double a0,double a1,double sz0,double sz1,double m,double s);
    double singlestepstate(double alphaold,double yy,double dd,double a0,double a1,double sz0,double sz1);
    void simulatestates();
    void simulatestatesadaptation(int k);
    void simulateparameters();
    double getlambda();
    
    void setprior(struct PriorStructSVJDUR prior);
    
};

#endif
