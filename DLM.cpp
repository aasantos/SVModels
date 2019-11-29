//
//  DLM.cpp
//  October2019
//
//  Created by António Alberto Santos on 15/10/2019.
//  Copyright © 2019 Antonio Santos. All rights reserved.
//

#include "DLM.hpp"

void dlmsimul(double sigmav,double phi1,double sigma1,
              double *y,double *alpha,int n,Random *random)
{
    double a = 0.0;
    for(int i=0;i<100;i++) a = phi1*a + sigma1*random->normal();
    for(int i=0;i<n;i++){
        a = phi1*a + sigma1*random->normal();
        alpha[i] = a;
        y[i] = a + sigmav*random->normal();
    }
}

void DLM::simulatestates()
{
    double sigmavsq = this->sigmav*this->sigmav;
    double sigma1sq = this->sigma1*this->sigma1;
    //
    double a0 = 0.0;
    double a1 = 0.0;
    double yy = 0.0;
    //
    this->alphafirst = (this->sigma1/sqrt(1.0 - this->phi1*this->phi1))*random->normal();
    for(int i=0;i<this->n;i++){
        if(i !=0 && i != (this->n - 1)){
            a0 = this->alpha[i-1];
            a1 = this->alpha[i+1];
            yy = this->y[i];
        }else if(i==0){
            a0 = this->alphafirst;
            a1 = this->alpha[i+1];
            yy = this->y[i];
        }else{
            a0 = this->alpha[i-1];
            a1 = this->alphalast;
            yy = this->y[i];
        }
        double mm = (yy*sigma1sq + sigmavsq*this->phi1*(a0 + a1))/
        (sigmavsq*(1.0 + this->phi1*this->phi1) + sigma1sq);
        double ss = sqrt((sigmavsq*sigma1sq)/
                         (sigmavsq*(1.0 + this->phi1*this->phi1) + sigma1sq));
        this->alpha[i] = mm + ss*this->random->normal();
    }
    this->alphalast = this->phi1*this->alpha[this->n - 2] +
    this->sigma1*this->random->normal();
}

void DLM::simulateparamters(){
    AR1Model *ar1 = new AR1Model(this->alpha,this->n,this->phi1,this->sigma1,this->random);
    //
    ar1->setupdatemu(this->updatephi1);
    ar1->setmudiffuse(this->phi1diffuse);
    ar1->setmupriortype(this->phi1priortype);
    ar1->setmuprior(this->phi1prior);
    ar1->setmubetaprior(this->phi1betaprior);
    //
    ar1->setupdatesigma(this->updatesigma1);
    ar1->setsigmadiffuse(this->sigma1diffuse);
    ar1->setsigmapriortype(this->sigma1priortype);
    ar1->setsigmaprior(this->sigma1prior);
    ar1->setsigmalognormalprior(this->sigma1lognormalprior);
    ar1->setsigmainvgaussianprior(this->sigma1invgaussianprior);
    //
    ar1->simulateparameters();
    this->phi1 = ar1->getmu();
    this->sigma1 = ar1->getsigma();
    //
    delete ar1;
    //
    //
    if(this->updatesigmav){
        double *err = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++) err[i] = this->y[i] - this->alpha[i];
        NormalBayesian *nb = new NormalBayesian(err,this->n,0.0,this->sigmav,this->random);
        nb->setupdatemu(false);
        nb->setupdatesigma(true);
        nb->setsigmadiffuse(this->sigmavdiffuse);
        nb->setsigmapriortype(this->sigmavpriortype);
        nb->setsigmaprior(this->sigmavprior);
        nb->setsigmalognormalprior(this->sigmavlognormalprior);
        nb->setsigmainvgaussianprior(this->sigmavinvgaussianprior);
        nb->simulateparameters();
        this->sigmav = nb->getsigma();
        delete nb;
        free(err);
    }
}

void DLM::setprior(struct PriorStruct prior)
{
    this->sigmavpriortype = prior.sigmavpriortype;
    this->phi1priortype = prior.phi1priortype;
    this->sigma1priortype = prior.sigma1priortype;
    //
    if(this->sigmavpriortype == 0){
        this->sigmavprior[0] = prior.sigmavprior[0];
        this->sigmavprior[1] = prior.sigmavprior[1];
    }else if(this->sigmavpriortype == 1){
        this->sigmavlognormalprior[0] = prior.sigmavprior[0];
        this->sigmavlognormalprior[1] = prior.sigmavprior[1];
    }else{
        this->sigmavinvgaussianprior[0] = prior.sigmavprior[0];
        this->sigmavinvgaussianprior[1] = prior.sigmavprior[1];
    }
    //
    if(this->phi1priortype == 0){
        this->phi1prior[0] = prior.phi1prior[0];
        this->phi1prior[1] = prior.phi1prior[1];
    }else{
        this->phi1betaprior[0] = prior.phi1prior[0];
        this->phi1betaprior[1] = prior.phi1prior[1];
    }
    //
    if(this->sigma1priortype == 0){
        this->sigma1prior[0] = prior.sigma1prior[0];
        this->sigma1prior[1] = prior.sigma1prior[1];
    }else if(this->sigma1priortype == 1){
        this->sigma1lognormalprior[0] = prior.sigma1prior[0];
        this->sigma1lognormalprior[1] = prior.sigma1prior[1];
    }else{
        this->sigma1invgaussianprior[0] = prior.sigma1prior[0];
        this->sigma1invgaussianprior[1] = prior.sigma1prior[1];
    }
}

double DLM::getalpha(int k)
{
    return this->alpha[k];
}

double DLM::getsigmav()
{
    return this->sigmav;
}

double DLM::getphi1()
{
    return this->phi1;
}

double DLM::getsigma1()
{
    return this->sigma1;
}

void DLM::estimate(const char *file,int niter,int nwarmup,int nstep)
{
    for(int i=0;i<100;i++) this->simulatestates();
    
    //warmup
    for(int i=0;i<nwarmup;i++){
        if(i%100 == 0) printf("Warmup iteration: %d / %d\n",i,nwarmup);
        this->simulatestates();
        this->simulateparamters();
    }
    
    //main iterations
    double *out = (double*)calloc(niter*3,sizeof(double));
    int iter = 0;
    for(int i=0;i<niter;i++){
        if(i%100 == 0) printf("Main iteration: %d / %d\n",i,niter);
        for(int k=0;k<nstep;k++){
            this->simulatestates();
            this->simulateparamters();
        }
        out[iter] = this->getsigmav();
        iter++;
        out[iter] = this->getphi1();
        iter++;
        out[iter] = this->getsigma1();
        iter++;
    }
    
    FILE *fp;
    fp = fopen(file,"wa");
    fprintf(fp,"sigmav phi1 sigma1\n");
    for(int i=0;i<niter;i++){
        fprintf(fp,"%.4f %.4f %.4f\n",out[i*3],out[i*3 + 1],out[i*3 + 2]);
    }
    fclose(fp);
    free(out);
}
