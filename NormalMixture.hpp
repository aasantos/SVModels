#ifndef NormalMixture_hpp
#define NormalMixture_hpp


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "Random.hpp"
#include "Stats.hpp"
#include "NormalBayesian.hpp"
#include "RegModel.hpp"

class NormalMixture{
protected:
    //
    int n;
    double *x;
    int    *S;
    double *Z;
    int njump;
    //
    double muj;
    double sigmaj;
    double theta;
    double sigma;
    //
    bool updatemuj;
    bool updatesigmaj;
    bool updatetheta;
    bool updatesigma;
    //
    bool mujdiffuse;
    bool sigmajdiffuse;
    bool sigmadiffuse;
    //
    int mujpriortype;
    int sigmajpriortype;
    int sigmapriortype;
    //
    double mujprior[2];
    double mujbetaprior[2];
    //
    double sigmajprior[2];
    double sigmajlognormalprior[2];
    double sigmajinvgaussianprior[2];
    //
    double thetaprior[2];
    //
    double sigmaprior[2] ;
    double sigmalognormalprior[2];
    double sigmainvgaussianprior[2];
    //
    Random *random;
    //
    //
public:
    NormalMixture(){}
    
    ~NormalMixture(){
         if(S){
            free(S);
        }
        if(Z){
            free(Z);
        }
    }
    
    
    NormalMixture(double *x,int n,double muj,double sigmaj,double theta,double sigma,Random *random)
    {
        this->n = n;
        //
        this->muj = muj;
        this->sigmaj = sigmaj;
        this->theta = theta;
        this->sigma = sigma;
        this->random = random;
        //
        this->x = x;
        this->S = (int*)malloc(this->n*sizeof(int));
        this->Z = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<n;i++){
            this->S[i] = 0;
            this->Z[i] = 0.0;
        }
        this->njump = 0;
        //
        this->updatemuj = true;
        this->updatesigmaj = true;
        this->updatetheta = true;
        this->updatesigma = true;
        //
        this->mujdiffuse = false;
        this->sigmajdiffuse = false;
        this->sigmadiffuse = false;
        //
        this->mujpriortype = 0;
        this->sigmajpriortype = 0;
        this->sigmapriortype = 0;
        //
        this->mujprior[0] = 0.0;
        this->mujprior[1] = 1.0;
        //
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
        //
    }
    
    double normalpdf(double x,double m,double s)
    {
        return (1.0/sqrt(2.0*3.141593*s*s))*exp(-0.5*(x - m)*(x-m)/(s*s));
    }
    
    
    void simultatetheta()
    {
        if(this->updatetheta){
            this->njump = 0;
            for(int i=0;i<this->n;i++){
                double w1 = this->theta*normalpdf(x[i],this->muj,sqrt(this->sigmaj*this->sigmaj + this->sigma*this->sigma));
                double w2 = (1.0 - this->theta)*normalpdf(x[i],0.0,this->sigma);
                double w = w1/(w1 + w2);
                this->S[i] = random->bernoulli(w);
                if(this->S[i] == 1){
                    this->njump++;
                    double mm = (this->x[i]*this->sigmaj*this->sigmaj + this->muj*this->sigma*this->sigma)/
                    (this->sigma*this->sigma + this->sigmaj*this->sigmaj);
                    double std = sqrt((this->sigma*this->sigma*this->sigmaj*this->sigmaj)/(this->sigma*this->sigma + this->sigmaj*this->sigmaj));
                    this->Z[i] = mm + std*random->normal();
                }
            }
            int n2 = this->n - this->njump;
            if(this->njump > 4){
                double a = (double)this->njump + thetaprior[0];
                double b = (double)n2 + thetaprior[1];
                this->theta = random->beta(a, b);
            }
        }
    }
    
    
    void simulatesigma()
    {
        if(this->updatesigma){
            double *err = (double*)malloc(this->n*sizeof(double));
            for(int i=0;i<this->n;i++) err[i] = this->x[i] - (double)this->S[i]*this->Z[i];
            NormalBayesian *nb = new NormalBayesian(err,this->n,0.0,this->sigma,this->random);
            nb->setupdatemu(false);
            nb->setsigmapriortype(this->sigmapriortype);
            nb->setsigmaprior(this->sigmaprior);
            nb->setsigmalognormalprior(this->sigmalognormalprior);
            nb->setsigmainvgaussianprior(this->sigmainvgaussianprior);
            nb->simulateparameters();
            this->sigma = nb->getsigma();
            delete nb;
            free(err);
        }
    }
    
    
    void simulatemujsigmaj()
    {
        if(this->njump > 4){
            double *xtemp = (double*)malloc(this->njump*sizeof(double));
            int j = 0;
            for(int i=0;i<this->n;i++){
                if(this->S[i] == 1){
                    xtemp[j] = this->Z[i];
                    j++;
                }
            }
            NormalBayesian *nb = new NormalBayesian(xtemp,this->njump,this->muj,this->sigmaj,this->random);
            nb->setupdatemu(this->updatemuj);
            nb->setupdatesigma(this->updatesigmaj);
            nb->setmudiffuse(this->mujdiffuse);
            nb->setsigmadiffuse(this->sigmajdiffuse);
            nb->setmupriortype(this->mujpriortype);
            nb->setsigmapriortype(this->sigmajpriortype);
            nb->setmuprior(this->mujprior);
            nb->setmubetaprior(this->mujbetaprior);
            nb->setsigmaprior(this->sigmajprior);
            nb->setsigmalognormalprior(this->sigmajlognormalprior);
            nb->setsigmainvgaussianprior(this->sigmajinvgaussianprior);
            nb->simulateparameters();
            this->muj = nb->getmu();
            this->sigmaj = nb->getsigma();
            free(xtemp);
            delete nb;
        }
    }
    
    
    void simulateparameters()
    {
        this->simultatetheta();
        this->simulatesigma();
        this->simulatemujsigmaj();
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
    double getsigma()
    {
        return this->sigma;
    }
    
    int getS(int k){
        return this->S[k];
    }
    
    double getZ(int k)
    {
        return this->Z[k];
    }
    
    int getnjump()
    {
        return this->njump;
    }
    //
    void setmujdiffuse(bool flag)
    {
        this->mujdiffuse = flag;
    }
    
    void setsigmajdiffuse(bool flag)
    {
        this->sigmajdiffuse = flag;
    }
    
    void setsigmadiffuse(bool flag)
    {
        this->sigmadiffuse = flag;
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
    void setupdatesigma(bool flag)
    {
        this->updatesigma = flag;
    }
    void setmujpriortype(int t)
    {
        this->mujpriortype = t;
    }
    void setsigmajpriortype(int t)
    {
        this->sigmajpriortype = t;
    }
    
    void setsigmapriortype(int t)
    {
        this->sigmapriortype = t;
    }
    //
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
    void setsigmaprior(double prior[2])
    {
        this->sigmaprior[0] = prior[0];
        this->sigmaprior[1] = prior[1];
    }
    void setsigmalognormalprior(double prior[2])
    {
        this->sigmalognormalprior[0] = prior[0];
        this->sigmalognormalprior[1] = prior[1];
    }
    void setsigmainvgaussianprior(double prior[2])
    {
        this->sigmainvgaussianprior[0] = prior[0];
        this->sigmainvgaussianprior[1] = prior[1];
    }

    
};


#endif
