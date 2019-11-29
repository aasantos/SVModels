//
//  SVJ2F.cpp
//  November2019
//
//  Created by António Alberto Santos on 07/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "SVJ2F.hpp"


void SVJ2F::firstsecondderiv(double yy,double a,double b, double a0,
                             double a1,double b0,double b1,double sz0,
                             double sz1,double fd[2],double sdinv[4])
{
    double sigmavsq = this->sigmav*this->sigmav;
    double sigma1sq = this->sigma1*this->sigma1;
    double sigma2sq = this->sigma2*this->sigma2;
    double dfa = -0.5 + 0.5*yy*yy*exp(-a - b)/sigmavsq +
    (a1 - this->phi1*a)*this->phi1/sigma1sq - (a - this->phi1*a0)/sigma1sq;
    double dfb = -0.5 + 0.5*yy*yy*exp(-a - b)/sigmavsq +
    (b1 - this->phi2*b - sz1)*this->phi2/sigma2sq - (b - this->phi2*b0 - sz0)/sigma2sq;
    fd[0] = dfa;
    fd[1] = dfb;
    //
    //
    double tt = yy*yy*exp(-a - b);
    double dfaa = 0.5*(tt*sigma1sq + 2.0*sigmavsq*(1.0 + this->phi1*this->phi1))/(sigmavsq*sigma1sq);
    double dfab = 0.5*tt/sigmavsq;
    double dfbb = 0.5*(tt*sigma2sq + 2.0*sigmavsq*(1.0 + this->phi2*this->phi2))/(sigmavsq*sigma2sq);
    double ddd = dfaa*dfbb - dfab*dfab;
    sdinv[0] = dfbb/ddd;
    sdinv[1] = -1.0*dfab/ddd;
    sdinv[2] = sdinv[1];
    sdinv[3] = dfaa/ddd;
}


void SVJ2F::newton(double yy,double a,double b,double a0,double a1,
                   double b0,double b1,double sz0,double sz1,double sol[2],double sdinv[4])
{
    double fd[2];
    sol[0] = a;
    sol[1] = b;
    this->firstsecondderiv(yy,sol[0], sol[1], a0, a1, b0, b1, sz0, sz1, fd, sdinv);
    int niter = 0;
    while(norm(fd) > 0.000001 && niter < 30){
        this->firstsecondderiv(yy, sol[0], sol[1], a0, a1, b0, b1, sz0, sz1, fd, sdinv);
        sol[0] = sol[0] + sdinv[0]*fd[0] + sdinv[1]*fd[1];
        sol[1] = sol[1] + sdinv[2]*fd[0] + sdinv[3]*fd[1];
        niter++;
    }
}

double SVJ2F::logpdfbivariatenormal(double x[2],double mm[2],double vcv[4])
{
    double rho = vcv[1]/sqrt(vcv[0]*vcv[3]);
    double mu1 = mm[0];
    double mu2 = mm[1];
    double s1sq = vcv[0];
    double s2sq = vcv[3];
    double err1 = x[0] - mu1;
    double err2 = x[1] - mu2;
    double zz = err1*err1/s1sq + err2*err2/s2sq -
    2.0*rho*err1*err2/sqrt(s1sq*s2sq);
    return -1.0*log(2.0*3.1415926535897932385*sqrt(s1sq)*sqrt(s1sq)*sqrt(1.0 - rho*rho))
    - 0.5*zz/(1.0 - rho*rho);
}

double SVJ2F::logpdfsv2f(double yy,double a0,double a,double a1,double b0,double b,double b1,double sz0,double sz1)
{
    double t1 = -0.5*a - 0.5*b - 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a + b));
    double t2 = -0.5*(a1 - this->phi1*a)*(a1 - this->phi1*a)/(this->sigma1*this->sigma1);
    double t3 = -0.5*(a - this->phi1*a0)*(a - this->phi1*a0)/(this->sigma1*this->sigma1);
    double t4 = -0.5*(b1 - this->phi2*b - sz1)*(b1 - this->phi2*b - sz1)/(this->sigma2*this->sigma2);
    double t5 = -0.5*(b - this->phi2*b0 - sz0)*(b - this->phi2*b0 - sz0)/(this->sigma2*this->sigma2);
    return t1 + t2 + t3 + t4 + t5;
}

double SVJ2F::metroprob(double anew,double bnew, double yy,
                 double a0,double a,double a1,
                 double b0,double b,double b1,double sz0,double sz1,double mm[2],double vcv[4])
{
    double t1 = logpdfsv2f(yy,a0,anew,a1,b0,bnew,b1,sz0,sz1);
    double t2 = logpdfsv2f(yy,a0,a,a1,b0,b,b1,sz0,sz1);
    double prop0[2] = {a,b};
    double prop1[2] = {anew,bnew};
    double t3 = logpdfbivariatenormal(prop0, mm, vcv);
    double t4 = logpdfbivariatenormal(prop1, mm, vcv);
    double tt = t1 - t2 + t3 - t4;
    if(tt < -15.0) tt = -15.0;
    if(tt > 0.0) tt = 0.0;
    return exp(tt);
}


void SVJ2F::singlestatesimulate(double yy,double a0,double a,double a1,
                         double b0,double b,double b1,double sz0,double sz1,
                         double abnew[2])
{
    double mm[2];
    double SS[4];
    this->newton(yy, a, b, a0, a1, b0, b1, sz0, sz1, mm, SS);
    double mmf[2];
    mmf[0] = (double)mm[0];
    mmf[1] = (double)mm[1];
    double SSf[4];
    SSf[0] = (double)SS[0];
    SSf[1] = (double)SS[1];
    SSf[2] = (double)SS[2];
    SSf[3] = (double)SS[3];
    double prop[2];
    random->bivariatenormal(mmf, SSf, prop);
    double mprob = metroprob((double)prop[0], (double)prop[1], yy, a0, a, a1,
                             b0, b, b1, sz0,sz1,
                             mm, SS);
    if(random->uniform() < mprob){
        abnew[0] = (double)prop[0];
        abnew[1] = (double)prop[1];
    }else{
        abnew[0] = a;
        abnew[1] = b;
    }
}

void SVJ2F::simulatestates()
{
    double a0 = 0.0;
    double a1 = 0.0;
    double b0 = 0.0;
    double b1 = 0.0;
    double yy = 0.0;
    double sz0 = 0.0;
    double sz1 = 0.0;
    this->alphafirst = this->sigma1*random->normal();
    this->betafirst = this->sigma2*random->normal();
    for(int i=0;i<this->n;i++){
        if(i !=0 && i != (this->n - 1)){
            a0 = this->alpha[i-1];
            a1 = this->alpha[i+1];
            sz0 = (double)this->S[i-1]*Z[i-1];
            sz1 = (double)this->S[i]*Z[i];
            b0 = this->beta[i-1];
            b1 = this->beta[i+1];
            yy = this->y[i];
        }else if(i==0){
            a0 = this->alphafirst;
            a1 = this->alpha[i+1];
            sz0 = 0.0;
            sz1 = (double)this->S[i]*Z[i];
            b0 = this->betafirst;
            b1 = this->beta[i+1];
            yy = this->y[i];
        }else{
            a0 = this->alpha[i-1];
            a1 = this->alphalast;
            sz0 = (double)this->S[i-1]*Z[i-1];
            sz1 = (double)this->S[i]*Z[i];
            b0 = this->beta[i-1];
            b1 = this->betalast;
            yy = this->y[i];
        }
        double abnew[2];
        this->singlestatesimulate(yy, a0, this->alpha[i], a1, b0, this->beta[i], b1, sz0, sz1, abnew);
        this->alpha[i] = abnew[0];
        this->beta[i] = abnew[1];
    }
    this->alphalast = this->phi1*this->alpha[this->n - 2] +
    this->sigma1*this->random->normal();
    this->betalast = this->phi2*this->beta[this->n - 2] +
    this->sigma2*this->random->normal();
}

void SVJ2F::singlestatesimulateadaptation(double yy,double a0,double a,double a1,
                         double b0,double b,double b1,double sz0,double sz1,
                         double abnew[2])
{
    double mm[2];
    double SS[4];
    this->newton(yy, a, b, a0, a1, b0, b1, sz0, sz1, mm, SS);
    double mmf[2];
    mmf[0] = (double)mm[0];
    mmf[1] = (double)mm[1];
    double SSf[4];
    SSf[0] = (double)SS[0];
    SSf[1] = (double)SS[1];
    SSf[2] = (double)SS[2];
    SSf[3] = (double)SS[3];
    double prop[2];
    random->bivariatenormal(mmf, SSf, prop);
    double mprob = metroprob((double)prop[0], (double)prop[1], yy, a0, a, a1,
                             b0, b, b1, sz0,sz1,
                             mm, SS);
    mprob = 1.0;
    if(random->uniform() < mprob){
        abnew[0] = (double)prop[0];
        abnew[1] = (double)prop[1];
    }else{
        abnew[0] = a;
        abnew[1] = b;
    }
}

void SVJ2F::simulatestatesadaption(int k)
{
    for(int j=0;j<k;j++){
            double a0 = 0.0;
            double a1 = 0.0;
            double b0 = 0.0;
            double b1 = 0.0;
            double yy = 0.0;
            double sz0 = 0.0;
            double sz1 = 0.0;
            this->alphafirst = this->sigma1*random->normal();
            this->betafirst = this->sigma2*random->normal();
            for(int i=0;i<this->n;i++){
                if(i !=0 && i != (this->n - 1)){
                    a0 = this->alpha[i-1];
                    a1 = this->alpha[i+1];
                    sz0 = (double)this->S[i-1]*Z[i-1];
                    sz1 = (double)this->S[i]*Z[i];
                    b0 = this->beta[i-1];
                    b1 = this->beta[i+1];
                    yy = this->y[i];
                }else if(i==0){
                    a0 = this->alphafirst;
                    a1 = this->alpha[i+1];
                    sz0 = 0.0;
                    sz1 = (double)this->S[i]*Z[i];
                    b0 = this->betafirst;
                    b1 = this->beta[i+1];
                    yy = this->y[i];
                }else{
                    a0 = this->alpha[i-1];
                    a1 = this->alphalast;
                    sz0 = (double)this->S[i-1]*Z[i-1];
                    sz1 = (double)this->S[i]*Z[i];
                    b0 = this->beta[i-1];
                    b1 = this->betalast;
                    yy = this->y[i];
                }
                double abnew[2];
                this->singlestatesimulateadaptation(yy, a0, this->alpha[i], a1, b0, this->beta[i], b1, sz0, sz1, abnew);
                this->alpha[i] = abnew[0];
                this->beta[i] = abnew[1];
            }
            this->alphalast = this->phi1*this->alpha[this->n - 2] +
            this->sigma1*this->random->normal();
            this->betalast = this->phi2*this->beta[this->n - 2] +
            this->sigma2*this->random->normal();
    }
}

void SVJ2F::simulateparameters()
{
    AR1JModel *ar1j = new AR1JModel(this->beta,this->n,this->phi2,this->sigma2,this->muj,this->sigmaj,this->theta,this->random);
    //
    ar1j->setupdatemu(this->updatephi2);
    ar1j->setmudiffuse(this->phi2diffuse);
    ar1j->setmupriortype(this->phi2priortype);
    ar1j->setmuprior(this->phi2prior);
    ar1j->setmubetaprior(this->phi2betaprior);
    //
    ar1j->setupdatesigma(this->updatesigma2);
    ar1j->setsigmadiffuse(this->sigma2diffuse);
    ar1j->setsigmapriortype(this->sigma2priortype);
    ar1j->setsigmaprior(this->sigma2prior);
    ar1j->setsigmalognormalprior(this->sigma2lognormalprior);
    ar1j->setsigmainvgaussianprior(this->sigma2invgaussianprior);
    //
    ar1j->setupdatemuj(this->updatemuj);
    ar1j->setmujdiffuse(this->mujdiffuse);
    ar1j->setmujpriortype(this->mujpriortype);
    ar1j->setmujprior(this->mujprior);
    ar1j->setmujbetaprior(this->mujbetaprior);
    //
    ar1j->setupdatesigmaj(this->updatesigmaj);
    ar1j->setsigmajdiffuse(this->sigmajdiffuse);
    ar1j->setsigmajpriortype(this->sigmajpriortype);
    ar1j->setsigmajprior(this->sigmajprior);
    ar1j->setsigmajlognormalprior(this->sigmajlognormalprior);
    ar1j->setsigmajinvgaussianprior(this->sigmajinvgaussianprior);
    //
    ar1j->setupdatetheta(this->updatetheta);
    ar1j->setthetaprior(this->thetaprior);
    //
    ar1j->simulateparameters();
    //
    this->phi2   =  ar1j->getmu();
    this->sigma2 =  ar1j->getsigma();
    this->muj    =  ar1j->getmuj();
    this->sigmaj =  ar1j->getsigmaj();
    this->theta  =  ar1j->gettheta();
    //
    for(int i=0;i<(this->n-1);i++){
        this->S[i] = ar1j->getS(i);
        this->Z[i] = ar1j->getZ(i);
    }
    this->S[this->n - 1] = 0;
    this->Z[this->n - 1] = 0.0;
    //
    delete ar1j;
    //
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
    if(this->updatesigmav){
        double *err = (double*)malloc(this->n*sizeof(double));
        for(int i=0;i<this->n;i++) err[i] = this->y[i]/exp(0.5*this->alpha[i] + 0.5*this->beta[i]);
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
    //
}


double SVJ2F::getmuj()
{
    return this->muj;
}

double SVJ2F::getsigmaj()
{
    return this->sigmaj;
}

double SVJ2F::gettheta()
{
    return this->theta;
}


void SVJ2F::setprior(struct PriorStructSVJ2F prior)
{
    this->sigmavpriortype = prior.sigmavpriortype;
    this->phi1priortype = prior.phi1priortype;
    this->sigma1priortype = prior.sigma1priortype;
    this->phi2priortype = prior.phi2priortype;
    this->sigma2priortype = prior.sigma2priortype;
    this->mujpriortype = prior.mujpriortype;
    this->sigmajpriortype = prior.sigmajpriortype;
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
       //
       //
          if(this->phi2priortype == 0){
              this->phi2prior[0] = prior.phi2prior[0];
              this->phi2prior[1] = prior.phi2prior[1];
          }else{
              this->phi2betaprior[0] = prior.phi2prior[0];
              this->phi2betaprior[1] = prior.phi2prior[1];
          }
          //
          if(this->sigma2priortype == 0){
              this->sigma2prior[0] = prior.sigma2prior[0];
              this->sigma2prior[1] = prior.sigma2prior[1];
          }else if(this->sigma2priortype == 1){
              this->sigma2lognormalprior[0] = prior.sigma2prior[0];
              this->sigma2lognormalprior[1] = prior.sigma2prior[1];
          }else{
              this->sigma2invgaussianprior[0] = prior.sigma2prior[0];
              this->sigma2invgaussianprior[1] = prior.sigma2prior[1];
          }
    
    if(this->mujpriortype == 0){
        this->mujprior[0] = prior.mujprior[0];
        this->mujprior[1] = prior.mujprior[1];
    }else{
        this->mujbetaprior[0] = prior.mujprior[0];
        this->mujbetaprior[1] = prior.mujprior[1];
    }
    //
    if(this->sigmajpriortype == 0){
        this->sigmajprior[0] = prior.sigmajprior[0];
        this->sigmajprior[1] = prior.sigmajprior[1];
    }else if(this->sigmajpriortype == 1){
        this->sigmajlognormalprior[0] = prior.sigmajprior[0];
        this->sigmajlognormalprior[1] = prior.sigmajprior[1];
    }else{
        this->sigmajinvgaussianprior[0] = prior.sigmajprior[0];
        this->sigmajinvgaussianprior[1] = prior.sigmajprior[1];
    }
    //
    this->thetaprior[0] = prior.thetaprior[0];
    this->thetaprior[1] = prior.thetaprior[1];
}
