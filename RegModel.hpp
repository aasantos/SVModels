#ifndef RegModel_hpp
#define RegModel_hpp

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "Random.hpp"
#include "Stats.hpp"



class RegModel{
protected:
    int n;
    double *y;
    double *x;
    //
    double mu;
    double sigma;
    //
    bool updatemu;
    bool updatesigma;
    //
    bool mudiffuse;
    bool sigmadiffuse;
    //
    int mupriortype; // 0-normal 1-beta
    double muprior[2];
    double mubetaprior[2];
    //
    int sigmapriortype; //0-invgamma 1-lognormal 2-invgaussian
    double sigmaprior[2];
    double sigmalognormalprior[2];
    double sigmainvgaussianprior[2];
    //
    Random *random;
    
public:
    RegModel(){}
    
    ~RegModel(){
    }
    
    RegModel(double *y,double *x,int n,double mu,double sigma,Random *random)
    {
        this->n = n;
        this->mu = mu;
        this->sigma = sigma;
        //
        this->x = x;
        this->y = y;
        this->random = random;
        
        this->updatemu = true;
        this->updatesigma = true;
        //
        this->mudiffuse = false;
        this->sigmadiffuse = false;
        //
        this->mupriortype = 0; // 0-normal 1-beta
        this->muprior[0] = 0.0;
        this->muprior[1] = 1.0;
        this->mubetaprior[0] = 20.0;
        this->mubetaprior[1] = 2.0;
        //
        this->sigmapriortype = 0; //0-invgamma 1-lognormal 2-invgaussian
        this->sigmaprior[0] = 2.5;
        this->sigmaprior[1] = 0.0025;
        this->sigmalognormalprior[0] = -1.0;
        this->sigmalognormalprior[1] = 0.5;
        this->sigmainvgaussianprior[0] = 0.5;
        this->sigmainvgaussianprior[1] = 0.5;
    }
    
    void simulatemu()
    {
        if(this->updatemu){
            double num = 0.0;
            double den = 0.0;
            for(int i=0;i<n;i++){
                num += this->y[i]*this->x[i];
                den += this->x[i]*this->x[i];
            }
            double betaols = num/den;
            if(this->mudiffuse){
                //printf("Simulating from diffuse..\n");
                this->mu = betaols + (this->sigma/sqrt(den))*random->normal();
            }else{
                if(this->mupriortype == 0){
                    //printf("Simulating from the normal distribution..\n");
                    double meansample = betaols;
                    double varsample = this->sigma*this->sigma/den;
                    double meanpost = (meansample*this->muprior[1]*this->muprior[1] +
                                       this->muprior[0]*varsample)/(this->muprior[1]*this->muprior[1] + varsample);
                    double varpost = (varsample*this->muprior[1]*this->muprior[1])/(this->muprior[1]*this->muprior[1] + varsample);
                    this->mu = meanpost + sqrt(varpost)*random->normal();
                }else{
                    //printf("Simulating from the beta distribution in RegModel %2.3f %2.3f..\n",this->mubetaprior[0],this->mubetaprior[1]);
                    double meansample = betaols;
                    if(meansample > 0.9999) meansample = 0.9999;
                    if(meansample < 0.0001) meansample = 0.01;
                    double varsample = this->sigma*this->sigma/den;
                    double mutemp = normalbetasimul(this->mu,meansample,sqrt(varsample), this->mubetaprior[0], this->mubetaprior[1]);
                    if(mutemp > 0.9999) mutemp = 0.9999;
                    if(mutemp < 0.0001) mutemp = 0.01;
                    this->mu = mutemp;
                }
            }
        }
    }
    
    
    void simulatesigma()
    {
        if(this->updatesigma){
            double rss = 0.0;
            for(int i=0;i<this->n;i++){
                rss += (this->y[i] - this->mu*this->x[i])*(this->y[i] - this->mu*this->x[i]);
            }
            if(this->sigmadiffuse){
                this->sigma = 1.0/sqrt(random->gamma(0.5*(double)this->n, 0.5*rss));
            }else{
                if(this->sigmapriortype == 0){
                    double p1 = 0.5*(double)this->n + this->sigmaprior[0];
                    double p2 = 0.5*rss + this->sigmaprior[1];
                    this->sigma = 1.0/sqrt(random->gamma(p1, p2));
                }else if(this->sigmapriortype == 1){
                    double a = 0.5*(double)this->n;
                    double b = 0.5*rss;
                    this->sigma = invgammalognormalsimul(this->sigma, a, b, this->sigmalognormalprior[0], this->sigmalognormalprior[1]);
                }else{
                    double a = 0.5*(double)this->n;
                    double b = 0.5*rss;
                    this->sigma = invgammainvgaussiansimul(this->sigma, a, b, this->sigmainvgaussianprior[0], this->sigmainvgaussianprior[1]);
                }
            }
        }
    }
    
    void simulateparameters()
    {
        this->simulatemu();
        this->simulatesigma();
    }
    
    double getmu(){
        return this->mu;
    }
    
    double getsigma(){
        return this->sigma;
    }
    
    void setupdatemu(bool flag)
    {
        this->updatemu = flag;
    }
    
    void setupdatesigma(bool flag)
    {
        this->updatesigma = flag;
    }
    
    void setmudiffuse(bool flag)
    {
        this->mudiffuse = flag;
    }
    
    void setsigmadiffuse(bool flag)
    {
        this->sigmadiffuse = flag;
    }
    
    void setmupriortype(int type)
    {
        this->mupriortype = type;
    }
    
    void setsigmapriortype(int type)
    {
        this->sigmapriortype = type;
    }
    
    void setmuprior(double mup[2])
    {
        this->muprior[0] = mup[0];
        this->muprior[1] = mup[1];
    }
    
    void setmubetaprior(double mup[2])
    {
        this->mubetaprior[0] = mup[0];
        this->mubetaprior[1] = mup[1];
    }
    
    void setsigmaprior(double sigmap[2])
    {
        this->sigmaprior[0] = sigmap[0];
        this->sigmaprior[1] = sigmap[1];
    }
    
    void setsigmalognormalprior(double sigmap[2])
    {
        this->sigmalognormalprior[0] = sigmap[0];
        this->sigmalognormalprior[1] = sigmap[1];
    }
    
    void setsigmainvgaussianprior(double sigmap[2])
    {
        this->sigmainvgaussianprior[0] = sigmap[0];
        this->sigmainvgaussianprior[1] = sigmap[1];
    }
    
    //auxiliary functions
    double lognormalpdf(double x,double m,double s)
    {
        return -0.5*(x - m)*(x - m)/(s*s);
    }
    
    double logpdfbetaunnorm(double x,double a,double b)
    {
        double result = 0.0;
        if(x > 0.0 && x < 1.0){
            result =  (a - 1.0)*log(x) + (b - 1.0)*log(1.0 - x);
        }else{
            result =  -99999.0;
        }
        return result;
    }
    
    double normalbetasimul(double phihold,double mnorm,double stdnorm,double a,double b)
    {
        double phinew = mnorm + stdnorm*random->normal();
        if(phinew > 0.99999) return phinew = 0.9999;
        double t1 = logpdfbetaunnorm(phinew,a,b);
        double t2 = logpdfbetaunnorm(phihold,a,b);
        double t3 = lognormalpdf(phinew, mnorm, stdnorm);
        double t4 = lognormalpdf(phihold, mnorm, stdnorm);
        double tt = t1 - t2 + t3 - t4;
        if(tt > 0.0) tt = 0.1;
        if(tt < -15.0) tt = -15.0;
        double metroprob = fminf(1.0,expf(tt));
        if(random->uniform() < metroprob){
            return phinew;
        }else{
            return phihold;
        }
    }
    
    double logpdflognormalunnorm(double x,double m,double s)
    {
        double err = (log(x) - m);
        return -1.0*log(x) - 0.5*err*err/(s*s);
    }
    
    double invgammalognormalsimul(double sigmaold,double ag,double bg,double mulog,double stdlog)
    {
        double sigmanew = 1.0/sqrt(random->gamma(ag, bg));
        double t1 = logpdflognormalunnorm(sigmanew*sigmanew, mulog, stdlog);
        double t2 = logpdflognormalunnorm(sigmaold*sigmaold, mulog, stdlog);
        double tt = t1 - t2;
        if(tt > 0.0){
            tt = 0.1;
        }
        if(tt < -15.0){
            tt = -15.0;
        }
        double metroprob = fminf(1.0,expf(tt));
        if(random->uniform() < metroprob)
        {
            return sigmanew;
        }else{
            return sigmaold;
        }
    }
    
    double loginvgaussian(double x,double m,double lambda)
    {
        double err = (x - m);
        return -1.5*log(x) - 0.5*lambda*err*err/(m*m*x);
    }
    
    double invgammainvgaussiansimul(double sigmaold,double ag,double bg,double m,double lambda)
    {
        double sigmanew = 1.0/sqrt(random->gamma(ag, bg));
        double t1 = loginvgaussian(sigmanew*sigmanew, m, lambda);
        double t2 = loginvgaussian(sigmaold*sigmaold, m, lambda);
        double tt = t1 - t2;
        if(tt > 0.0){
            tt = 0.1;
        }
        if(tt < -15.0){
            tt = -15.0;
        }
        double metroprob = fminf(1.0,expf(tt));
        if(random->uniform() < metroprob)
        {
            return sigmanew;
        }else{
            return sigmaold;
        }
    }
        
};


#endif
