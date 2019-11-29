//
//  DLM.hpp
//  October2019
//
//  Created by António Alberto Santos on 15/10/2019.
//  Copyright © 2019 Antonio Santos. All rights reserved.
//

#ifndef DLM_hpp
#define DLM_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Random.hpp"
#include "func.hpp"
#include "Stats.hpp"
#include "NormalBayesian.hpp"
#include "RegModel.hpp"
#include "AR1Model.hpp"

void dlmsimul(double sigmav,double phi1,double sigma1,
              double *y,double *alpha,int n,Random *random);



struct PriorStruct
{
    int sigmavpriortype;
    int phi1priortype;
    int sigma1priortype;
    double sigmavprior[2];
    double phi1prior[2];
    double sigma1prior[2];
};

class DLM{
//
protected:
    int n;
    double *y;
    double *alpha;
    double sigmav;
    double phi1;
    double sigma1;
    //
    bool updatesigmav;
    bool updatephi1;
    bool updatesigma1;
    //
    bool sigmavdiffuse;
    bool phi1diffuse;
    bool sigma1diffuse;
    //
    int phi1priortype; // 0-normal 1-beta
    double phi1prior[2];
    double phi1betaprior[2];
    //
    int sigma1priortype; //0-invgamma 1-lognormal 2-invgaussian
    double sigma1prior[2];
    double sigma1lognormalprior[2];
    double sigma1invgaussianprior[2];
    //
    int sigmavpriortype; //0-invgamma 1-lognormal 2-invgaussian
    double sigmavprior[2];
    double sigmavlognormalprior[2];
    double sigmavinvgaussianprior[2];
    //
    double alphafirst;
    double alphalast;
    //
    Random *random;
//
public:
    DLM(double *y,int n,double sigmav,double phi1,double sigma1,Random *random)
    {
        this->y = y;
        this->n = n;
        this->sigmav = sigmav;
        this->phi1 = phi1;
        this->sigma1 = sigma1;
        this->alpha = (double*)calloc(n,sizeof(double));
        this->random = random;
        //
        this->updatesigmav = true;
        this->updatephi1 = true;
        this->updatesigma1 = true;
        //
        this->sigmavdiffuse = false;
        this->phi1diffuse = false;
        this->sigma1diffuse = false;
        //
        this->phi1priortype = 1; // 0-normal 1-beta
        this->phi1prior[0] = 0.9;
        this->phi1prior[1] = 0.5;
        this->phi1betaprior[0] = 95.0;
        this->phi1betaprior[1] =  5.0;
        //
        this->sigma1priortype = 1; //0-invgamma 1-lognormal 2-invgaussian
        this->sigma1prior[0] = 2.5;
        this->sigma1prior[1] = 0.025;
        this->sigma1lognormalprior[0] = -2.0;
        this->sigma1lognormalprior[1] = 0.5;
        this->sigma1invgaussianprior[0]  = 0.5;
        this->sigma1invgaussianprior[1] = 0.5;
        //
        this->sigmavpriortype = 1; //0-invgamma 1-lognormal 2-invgaussian
        this->sigmavprior[0] = 2.5;
        this->sigmavprior[1] = 0.025;
        this->sigmavlognormalprior[0]  = -2.0;
        this->sigmavlognormalprior[1]  = 0.5;
        this->sigmavinvgaussianprior[0] = 0.5;
        this->sigmavinvgaussianprior[1] = 0.5;
        //
        this->alphafirst = 0.0;
        this->alphalast = 0.0;
    }
    
    ~DLM()
    {
        if(alpha){
            free(alpha);
        }
    }
    //
    void simulatestates();
    void simulateparamters();
    //
    void setprior(struct PriorStruct prior);
    //
    double getalpha(int k);
    //
    double getsigmav();
    double getphi1();
    double getsigma1();
    //
    void estimate(const char *file,int niter,int nwarmup,int nstep);
};


#endif /* DLM_hpp */
