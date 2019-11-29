#ifndef AR1JModel_hpp
#define AR1JModel_hpp

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "Random.hpp"
#include "Stats.hpp"
#include "NormalBayesian.hpp"
#include "NormalMixture.hpp"
#include "RegModel.hpp"


class AR1JModel{
protected:
    int n;
    int m;
    double *x2;
    double *x1;
    int    *S;
    double *Z;
    double mu;
    double sigma;
    double muj;
    double sigmaj;
    double theta;
    //
    bool updatemu;
    bool updatesigma;
    bool updatemuj;
    bool updatesigmaj;
    bool updatetheta;
    //
    bool mudiffuse;
    bool sigmadiffuse;
    bool mujdiffuse;
    bool sigmajdiffuse;
    bool thetadiffuse;
    
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
    int mujpriortype; // 0-normal 1-beta
    double mujprior[2];
    double mujbetaprior[2];
    //
    int sigmajpriortype; //0-invgamma 1-lognormal 2-invgaussian
    double sigmajprior[2];
    double sigmajlognormalprior[2];
    double sigmajinvgaussianprior[2];
    //
    double thetaprior[2];
    //
    Random *random;
    
public:
    AR1JModel(){}
    
    ~AR1JModel(){
        if(x2){
            free(x2);
        }
        if(x1){
            free(x1);
        }
        if(S){
            free(S);
        }
        if(Z){
            free(Z);
        }
    }
    
    AR1JModel(double *x,int m,double mu,double sigma,double muj,double sigmaj,double theta,Random *random)
    {
        this->m = m;
        this->n = m - 1;
        //
        this->mu = mu;
        this->sigma = sigma;
        this->muj = muj;
        this->sigmaj = sigmaj;
        this->theta = theta;
        this->random = random;
        //
        this->x2 = (double*)malloc(this->n*sizeof(double));
        this->x1 = (double*)malloc(this->n*sizeof(double));
        this->S = (int*)malloc(this->n*sizeof(int));
        this->Z = (double*)malloc(this->n*sizeof(double));
        //
        for(int i=0;i<this->n;i++){
            this->x2[i] = x[i+1];
            this->x1[i] = x[i];
            this->S[i] = 0;
            this->Z[i] = 0.0;
        }
        //
        this->updatemu = true;
        this->updatemuj = true;
        this->updatesigmaj = true;
        this->updatetheta = true;
        this->updatesigma = true;
        //
        this->mujdiffuse = false;
        this->sigmajdiffuse = false;
        this->thetadiffuse = false;
        this->sigmadiffuse = false;
        //
        this->mupriortype = 0;
        this->mujpriortype = 0;
        this->sigmajpriortype = 0;
        this->sigmapriortype = 0;
        //
        this->muprior[0] = 0.0;
        this->muprior[1] = 1.0;
        //
        this->mubetaprior[0] = 20.0;
        this->mubetaprior[1] = 2.0;
        //
        this->mujprior[0] = 0.0;
        this->mujprior[1] = 1.0;
        this->mujbetaprior[0] = 20.0;
        this->mujbetaprior[1] = 2.0;
        //
        this->sigmajprior[0] = 2.5;
        this->sigmajprior[1] = 0.025;
        this->sigmajlognormalprior[0]  = -1.0;
        this->sigmajlognormalprior[1]  = 0.5;
        this->sigmajinvgaussianprior[0] = 0.5;
        this->sigmajinvgaussianprior[1] = 0.5;
        //
        this->thetaprior[0] =  2.0;
        this->thetaprior[1] = 198.0;
        //
        this->sigmaprior[0] = 2.5;
        this->sigmaprior[1] = 0.025;
        this->sigmalognormalprior[0] = -1.0;
        this->sigmalognormalprior[1] = 0.5;
        this->sigmainvgaussianprior[0] = 0.5;
        this->sigmainvgaussianprior[1] = 0.5;
    }
    
    
    void simulateparameters()
    {
        double *xx = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++) xx[i] = this->x2[i] - this->mu*this->x1[i];
        NormalMixture *nm = new NormalMixture(xx,this->n,this->muj,this->sigmaj,this->theta,this->sigma,this->random);
        //
        nm->setupdatemuj(this->updatemuj);
        nm->setupdatesigmaj(this->updatesigmaj);
        nm->setupdatetheta(this->updatetheta);
        nm->setupdatesigma(this->updatesigma);
        //
        nm->setmujdiffuse(this->mujdiffuse);
        nm->setsigmajdiffuse(this->sigmajdiffuse);
        nm->setsigmadiffuse(this->sigmajdiffuse);
        //
        nm->setmujpriortype(this->mujpriortype);
        nm->setsigmajpriortype(this->sigmajpriortype);
        nm->setsigmapriortype(this->sigmapriortype);
        //
        nm->setmujprior(this->mujprior);
        nm->setmujbetaprior(this->mujbetaprior);
        //
        nm->setsigmajprior(this->sigmajprior);
        nm->setsigmajlognormalprior(this->sigmajlognormalprior);
        nm->setsigmajinvgaussianprior(this->sigmajinvgaussianprior);
        //
        nm->setsigmaprior(this->sigmaprior);
        nm->setsigmalognormalprior(this->sigmalognormalprior);
        nm->setsigmainvgaussianprior(this->sigmainvgaussianprior);
        //
        nm->setthetaprior(this->thetaprior);
        //
        nm->simulateparameters();
        this->muj = nm->getmuj();
        this->sigmaj = nm->getsigmaj();
        this->theta = nm->gettheta();
        this->sigma = nm->getsigma();
        //
        for(int i=0;i<this->n;i++){
            this->S[i] = nm->getS(i);
            this->Z[i] = nm->getZ(i);
        }
        //
        double *xx1 = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++) xx1[i] = this->x2[i] - (double)this->S[i]*this->Z[i];
        RegModel *rm = new RegModel(xx1,this->x1,this->n,this->mu,this->sigma,this->random);
        rm->setupdatemu(true);
        rm->setmudiffuse(this->mudiffuse);
        rm->setupdatesigma(false);
        rm->setmupriortype(this->mupriortype);
        rm->setmuprior(this->muprior);
        rm->setmubetaprior(this->mubetaprior);
        rm->simulateparameters();
        this->mu = rm->getmu();
        
        //
        delete rm;
        delete nm;
        free(xx);
        free(xx1);
    }
    
    
    //getters
    double getmuj()
    {
        return this->muj;
    }
    double getsigmaj()
    {
        return this->sigmaj;
    }
    double gettheta()
    {
        return this->theta;
    }
    
    double getmu(){
        return this->mu;
    }
    
    double getsigma(){
        return this->sigma;
    }
    
    
    int getS(int k){
        return this->S[k];
    }
    
    double getZ(int k)
    {
        return this->Z[k];
    }
    
    void setupdatemu(bool flag)
    {
        this->updatemu = flag;
    }
    
    void setupdatesigma(bool flag)
    {
        this->updatesigma = flag;
    }
    
    void setupdatemuj(bool flag)
    {
        this->updatemuj = flag;
    }
    void setupdatesigmaj(bool flag)
    {
        this->updatesigmaj = flag;
    }
    void setupdatetheta(bool flag)
    {
        this->updatetheta = flag;
    }
    
    void setmujdiffuse(bool flag)
    {
        this->mujdiffuse = flag;
    }
    
    void setsigmajdiffuse(bool flag)
    {
        this->sigmajdiffuse = flag;
    }
    
    void setthetadiffuse(bool flag)
    {
        this->thetadiffuse = flag;
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
    
    void setmujpriortype(int t)
    {
        this->mujpriortype = t;
    }
    void setsigmajpriortype(int t)
    {
        this->sigmajpriortype = t;
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
    
    void setmujprior(double prior[2])
    {
        this->mujprior[0] = prior[0];
        this->mujprior[1] = prior[1];
    }
    void setmujbetaprior(double prior[2])
    {
        this->mujbetaprior[0] = prior[0];
        this->mujbetaprior[1] = prior[1];
    }
    void setsigmajprior(double prior[2])
    {
        this->sigmajprior[0] = prior[0];
        this->sigmajprior[1] = prior[1];
    }
    void setsigmajlognormalprior(double prior[2])
    {
        this->sigmajlognormalprior[0] = prior[0];
        this->sigmajlognormalprior[1] = prior[1];
    }
    void setsigmajinvgaussianprior(double prior[2])
    {
        this->sigmajinvgaussianprior[0] = prior[0];
        this->sigmajinvgaussianprior[1] = prior[1];
    }
    void setthetaprior(double prior[2])
    {
        this->thetaprior[0] = prior[0];
        this->thetaprior[1] = prior[1];
    }
};


#endif
