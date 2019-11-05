#ifndef AR1Model_hpp
#define AR1Model_hpp

//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "Random.hpp"
#include "Stats.hpp"
#include "NormalBayesian.hpp"
#include "RegModel.hpp"


class AR1Model{    
protected:
    int n;
    int m;
    double *x2;
    double *x1;
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
    AR1Model(){}
    
    ~AR1Model(){
        if(x2){
            free(x2);
        }
        if(x1){
            free(x1);
        }
    }
    
    AR1Model(double *x,int m,double mu,double sigma,Random *random)
    {
        this->m = m;
        this->n = m - 1;
        //
        this->mu = mu;
        this->sigma = sigma;
        this->random = random;
        //
        this->x2 = (double*)malloc(this->n*sizeof(double));
        this->x1 = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++){
            this->x2[i] = x[i+1];
            this->x1[i] = x[i];
        }
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
    
    
    void simulateparameters()
    {
        RegModel *rm = new RegModel(this->x2,this->x1,this->n,this->mu,this->sigma,this->random);
        if(this->updatemu){
            rm->setupdatemu(true);
            if(this->mudiffuse){
                rm->setmudiffuse(true);
            }else{
                rm->setmupriortype(this->mupriortype);
                rm->setmuprior(this->muprior);
                rm->setmubetaprior(this->mubetaprior);
            }
        }else{
            rm->setupdatemu(false);
        }
        //
        if(this->updatesigma)
        {
            rm->setupdatesigma(true);
            if(this->sigmadiffuse)
            {
                rm->setsigmadiffuse(true);
            }else{
                rm->setsigmapriortype(this->sigmapriortype);
                rm->setsigmaprior(this->sigmaprior);
                rm->setsigmalognormalprior(this->sigmalognormalprior);
                rm->setsigmainvgaussianprior(this->sigmainvgaussianprior);
            }
        }else{
            rm->setupdatesigma(false);
        }
        rm->simulateparameters();
        this->mu = rm->getmu();
        this->sigma = rm->getsigma();
        delete rm;
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
    
};



#endif
