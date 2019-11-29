#ifndef SVJ_hpp
#define SVJ_hpp

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "NormalMixture.hpp"
#include "AR1JModel.hpp"
#include "SV.hpp"



struct PriorStructSVJ
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
};


class SVJ: public SV{
    
protected:
    //
    int    *S;
    double  *Z;
    //
    double muj;
    double sigmaj;
    double theta;
    //
    bool updatemuj;
    bool updatesigmaj;
    bool updatetheta;
    //
    bool mujdiffuse;
    bool sigmajdiffuse;
    bool thetadiffuse;
    //
    int mujpriortype; // 0-normal 1-beta
    double mujprior[2] ;
    double mujbetaprior[2];
    //
    int sigmajpriortype; //0-invgamma 1-lognormal 2-invgaussian
    double sigmajprior[2];
    double sigmajlognormalprior[2];
    double sigmajinvgaussianprior[2];
    //
    double thetaprior[2];
        
public:
    ~SVJ(){
        if(S){
            free(S);
        }
        if(Z){
            free(Z);
        }
    }
    
    SVJ(double *x,int n,double sigmav,double mu,double sigma,Random *random,double muj,double sigmaj,double theta):
    SV(x,n,sigmav,mu,sigma,random)
    {
        this->muj = muj;
        this->sigmaj = sigmaj;
        this->theta = theta;
        //
        this->S = (int*)malloc(this->n*sizeof(int));
        this->Z = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++){
            this->S[i] = 0;
            this->Z[i] = 0.0;
        }
        
        this->updatemuj = true;
        this->updatesigmaj = true;
        this->updatetheta = true;
        //
        this->mujdiffuse = false;
        this->sigmajdiffuse = false;
        this->thetadiffuse = false;
        //
        this->mujpriortype = 1; // 0-normal 1-beta
        this->mujprior[0] = 0.8;
        this->mujprior[1] = 0.25;
        this->mujbetaprior[0] = 4.0;
        this->mujbetaprior[1] = 4.0;
        //
        this->sigmajpriortype = 1; //0-invgamma 1-lognormal 2-invgaussian
        this->sigmajprior[0] = 2.5;
        this->sigmajprior[1] = 0.025;
        this->sigmajlognormalprior[0] = -0.5;
        this->sigmajlognormalprior[1] =  0.5;
        this->sigmajinvgaussianprior[0] = 0.5;
        this->sigmajinvgaussianprior[1] = 0.5;
        //
        this->thetaprior[0] = 2.0;
        this->thetaprior[1] = 198.0;
        
        // initial values of the states
    }
    
    double df(double a,double yy,double a0,double a1,double sz0,double sz1);
    double ddf(double a,double yy);
    double meanAlpha(double yy,double a0,double a1,double sz0,double sz1);
    double varAlpha(double a,double yy);
    double logpdfsv(double a,double yy,double a0,double a1,double sz0,double sz1);
    double logpdfnormal(double a,double m,double s);
    double metroprob(double anew,double a,double yy,double a0,double a1,double sz0,double sz1,double m,double s);
    double singlestepstate(double alphaold,double yy,double a0,double a1,double sz0,double sz1);
    void simulatestates();
    void simulatestatesadaptation(int k);
    void simulateparameters();
    double getmuj();
    double getsigmaj();
    double gettheta();
    void setupdatemuj(bool flag);
    void setupdatesigmaj(bool flag);
    void setupdatetheta(bool flag);
    void setmujdiffuse(bool flag);
    void setsigmajdiffuse(bool flag);
    void setthetadiffuse(bool flag);
    void setpriror(struct PriorStructSVJ prior);
       
};

#endif 

