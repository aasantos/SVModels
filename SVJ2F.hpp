//
//  SVJ2F.hpp
//  November2019
//
//  Created by António Alberto Santos on 07/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef SVJ2F_hpp
#define SVJ2F_hpp

#include <stdio.h>
#include "AR1JModel.hpp"
#include "SV2F.hpp"

struct PriorStructSVJ2F
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
};


class SVJ2F: public SV2F{
//
//
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
//
//
public:
    ~SVJ2F(){
        if(S){
            free(S);
        }
        if(Z){
            free(Z);
        }
    }
    
    SVJ2F(double *y,int n,double sigmav,double phi1,double sigma1,
          Random *random,double phi2,double sigma2,
          double muj,double sigmaj,double theta):
    SV2F(y,n,sigmav,phi1,sigma1,random,phi2,sigma2)
    {
               this->muj = muj;
               this->sigmaj = sigmaj;
               this->theta = theta;
               //
               this->S = (int*)calloc(this->n,sizeof(int));
               this->Z = (double*)calloc(this->n,sizeof(double));
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
               this->sigmajlognormalprior[0] =  -0.5;
               this->sigmajlognormalprior[1] =  0.5;
               this->sigmajinvgaussianprior[0] = 0.5;
               this->sigmajinvgaussianprior[1] = 0.5;
               //
               this->thetaprior[0] = 2.0;
               this->thetaprior[1] = 198.0;
    }
    
    void firstsecondderiv(double yy,double a,double b, double a0,double a1,
                          double b0,double b1,double sz0,double sz1,double fd[2],double sdinv[4]);
    
    void newton(double yy,double a,double b,double a0,double a1,double b0,
                double b1,double sz0,double sz1,double sol[2],double sdinv[4]);
    
    double logpdfbivariatenormal(double x[2],double mm[2],double vcv[4]);
    //
    double logpdfsv2f(double yy,double a0,double a,double a1,
                      double b0,double b,double b1,double sz0,double sz1);
    //
    double metroprob(double anew,double bnew, double yy,
                     double a0,double a,double a1,
                     double b0,double b,double b1,double sz0,double sz1,
                     double mm[2],double vcv[4]);
    
    void singlestatesimulate(double y,double a0,double a,double a1,
                             double b0,double b,double b1,double sz0,double sz1,
                             double abnew[2]);
    //
    void simulatestates();
    
    void singlestatesimulateadaptation(double y,double a0,double a,double a1,
                             double b0,double b,double b1,double sz0,double sz1,
                             double abnew[2]);
    //
    void simulatestatesadaption(int k);
    
    
    void simulateparameters();
    //
    double getmuj();
    double getsigmaj();
    double gettheta();
    
    void setprior(struct PriorStructSVJ2F prior);
};

#endif /* SVJ2F_hpp */
