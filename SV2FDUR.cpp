//
//  SV2FDUR.cpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "SV2FDUR.hpp"


void svdurfirstsecondderiv(double sigmav,double phi1,double sigma1,
                           double phi2,double sigma2,
                           double lambda,double gam,
                           double y,double d,double a,double b,
                           double a0,double a1,double b0,double b1,
                           double fd[2],double sdinv[4])
{
    double sigmavsq = sigmav*sigmav;
    double sigma1sq = sigma1*sigma1;
    double sigma2sq = sigma2*sigma2;
    //
    double t1 = log(gam) + log(d);
    double t2 = exp(-0.5*a)*lambda;
    double t3 = -0.5*t2*t1 + 0.5*t2*digamma(t2);
    double t4 = -0.25*t2*t1 + 0.25*t2*digamma(t2) + 0.25*trigamma(t2)*t2*t2;
    double dfa = t3 - 0.5 + 0.5*y*y*exp(-a - b)/sigmavsq +
    (a1 - phi1*a)*phi1/sigma1sq - (a - phi1*a0)/sigma1sq;
    double dfb = -0.5 + 0.5*y*y*exp(-a - b)/sigmavsq +
    (b1 - phi2*b)*phi2/sigma2sq - (b - phi2*b0)/sigma2sq;
    //
    //
    double dfaa = 0.5*y*y*exp(-a - b)/sigmavsq + (1.0 + phi1*phi1)/sigma1sq + t4;
    double dfab = 0.5*y*y*exp(-a - b)/sigmavsq;
    double dfbb = 0.5*y*y*exp(-a - b)/sigmavsq + (1.0 + phi2*phi2)/sigma2sq;
    //
    //
    fd[0] = dfa;
    fd[1] = dfb;
    //
    double ddd = dfaa*dfbb - dfab*dfab;
    sdinv[0] = dfbb/ddd;
    sdinv[1] = -1.0*dfab/ddd;
    sdinv[2] = sdinv[1];
    sdinv[3] = dfaa/ddd;
}

void svdurnewton(double sigmav,double phi1,double sigma1,
                 double phi2,double sigma2,double lambda,double gam,
                 double y,double d,double a0,double a1,double b0,double b1,
                 double sol[2],double sdinv[4])
{
    double fd[2];
    sol[0] = 0.5*(a0 + a1);
    sol[1] = 0.5*(b0 + b1);
    svdurfirstsecondderiv(sigmav, phi1, sigma1,phi2,sigma2,lambda,gam,
                          y,d,sol[0],sol[1],a0,a1,b0,b1,fd,sdinv);
    int niter = 0;
    while(norm(fd) > 0.000001 && niter < 30){
        svdurfirstsecondderiv(sigmav, phi1, sigma1,phi2,sigma2,lambda,gam,
                              y,d,sol[0],sol[1],a0,a1,b0,b1,fd,sdinv);
        sol[0] = sol[0] + sdinv[0]*fd[0] + sdinv[1]*fd[1];
        sol[1] = sol[1] + sdinv[2]*fd[0] + sdinv[3]*fd[1];
        niter++;
    }
}


double SV2FDUR::logpdfbivariatenormal(double x[2],double mm[2],double vcv[4])
{
    double rho = vcv[1]/sqrt(vcv[0]*vcv[3]);
    double mu1 = mm[0];
    double mu2 = mm[1];
    double s1sq = vcv[0];
    double s2sq = vcv[3];
    double err1 = x[0] - mu1;
    double err2 = x[1] - mu2;
    double zz = err1*err1/s1sq + err2*err2/s2sq - 2.0*rho*err1*err2/sqrt(s1sq*s2sq);
    return -1.0*log(2.0*3.1415926535897932385*sqrt(s1sq)*sqrt(s1sq)*sqrt(1.0 - rho*rho)) - 0.5*zz/(1.0 - rho*rho);
}


double SV2FDUR::logpdfsv2fdur(double yy,double d,double a0,double a,double a1,double b0,double b,double b1)
{
    double tt1 = this->lambda*exp(-0.5*a);
    double tt2 = tt1*log(this->gam*d) - lgamma(tt1);
    double t1 = -0.5*a - 0.5*b - 0.5*yy*yy/(this->sigmav*this->sigmav*exp(a + b));
    double t2 = -0.5*(a1 - this->phi1*a)*(a1 - this->phi1*a)/(this->sigma1*this->sigma1);
    double t3 = -0.5*(a - this->phi1*a0)*(a - this->phi1*a0)/(this->sigma1*this->sigma1);
    double t4 = -0.5*(b1 - this->phi2*b)*(b1 - this->phi2*b)/(this->sigma2*this->sigma2);
    double t5 = -0.5*(b - this->phi2*b0)*(b - this->phi2*b0)/(this->sigma2*this->sigma2);
    return tt2 + t1 + t2 + t3 + t4 + t5;
}

double SV2FDUR::metroprob(double anew,double bnew, double y,double d,
                 double a0,double a,double a1,
                 double b0,double b,double b1,double mm[2],double vcv[4])
{
    double t1 = logpdfsv2fdur(y,d,a0,anew,a1,b0,bnew,b1);
    double t2 = logpdfsv2fdur(y,d,a0,a,a1,b0,b,b1);
    double prop0[2] = {a,b};
    double prop1[2] = {anew,bnew};
    double t3 = logpdfbivariatenormal(prop0, mm, vcv);
    double t4 = logpdfbivariatenormal(prop1, mm, vcv);
    double tt = t1 - t2 + t3 - t4;
    if(tt < -15.0) tt = -15.0;
    if(tt > 0.0) tt = 0.0;
    return exp(tt);
}

void SV2FDUR::singlestatesimulateAdaptation(double yy,double dd,double a0,double a,double a1,
                                  double b0,double b,double b1,
                                  double abnew[2])
{
    double mm[2];
    double SS[4];
    svdurnewton(this->sigmav,this->phi1,this->sigma1,
                this->phi2,this->sigma2,this->lambda,this->gam,yy,dd,a0,a1,b0,b1,mm,SS);
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
    double mprob = 1.0;
    if(random->uniform() < mprob){
        abnew[0] = (double)prop[0];
        abnew[1] = (double)prop[1];
    }else{
        abnew[0] = a;
        abnew[1] = b;
    }
}


void SV2FDUR::singlestatesimulate(double yy,double dd,double a0,double a,double a1,
                         double b0,double b,double b1,
                         double abnew[2])
{
    double mm[2];
    double SS[4];
    svdurnewton(this->sigmav,this->phi1,this->sigma1,
                this->phi2,this->sigma2,lambda,gam,yy,dd,a0,a1,b0,b1,mm,SS);
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
    double mprob = metroprob((double)prop[0], (double)prop[1],
                             yy, dd, a0, a, a1, b0, b, b1,mm, SS);
    if(random->uniform() < mprob){
        abnew[0] = (double)prop[0];
        abnew[1] = (double)prop[1];
    }else{
        abnew[0] = a;
        abnew[1] = b;
    }
}


void SV2FDUR::simulatestates()
{
    double a0 = 0.0;
    double a1 = 0.0;
    double b0 = 0.0;
    double b1 = 0.0;
    double yy = 0.0;
    double dd = 0.0;
    this->alphafirst = this->sigma1*random->normal();
    this->betafirst = this->sigma2*random->normal();
    for(int i=0;i<this->n;i++){
        if(i !=0 && i != (this->n - 1)){
            a0 = this->alpha[i-1];
            a1 = this->alpha[i+1];
            b0 = this->beta[i-1];
            b1 = this->beta[i+1];
            yy = this->y[i];
            dd = this->dur[i];
        }else if(i==0){
            a0 = this->alphafirst;
            a1 = this->alpha[i+1];
            b0 = this->betafirst;
            b1 = this->beta[i+1];
            yy = this->y[i];
            dd = this->dur[i];
        }else{
            a0 = this->alpha[i-1];
            a1 = this->alphalast;
            b0 = this->beta[i-1];
            b1 = this->betalast;
            yy = this->y[i];
            dd = this->dur[i];
        }
        double abnew[2];
        singlestatesimulate(yy, dd, a0,
                            this->alpha[i], a1, b0, this->beta[i], b1, abnew);
        this->alpha[i] = abnew[0];
        this->beta[i] = abnew[1];
    }
    this->alphalast = this->phi1*this->alpha[this->n - 2] +
    this->sigma1*this->random->normal();
    this->betalast = this->phi2*this->beta[this->n - 2] +
    this->sigma2*this->random->normal();
}


void SV2FDUR::simulatestatesAdaptation(int k)
{
    for(int jj=0;jj<k;jj++){
        double a0 = 0.0;
        double a1 = 0.0;
        double b0 = 0.0;
        double b1 = 0.0;
        double yy = 0.0;
        double dd = 0.0;
        this->alphafirst = this->sigma1*random->normal();
        this->betafirst = this->sigma2*random->normal();
        for(int i=0;i<this->n;i++){
            if(i !=0 && i != (this->n - 1)){
                a0 = this->alpha[i-1];
                a1 = this->alpha[i+1];
                b0 = this->beta[i-1];
                b1 = this->beta[i+1];
                yy = this->y[i];
                dd = this->dur[i];
            }else if(i==0){
                a0 = this->alphafirst;
                a1 = this->alpha[i+1];
                b0 = this->betafirst;
                b1 = this->beta[i+1];
                yy = this->y[i];
                dd = this->dur[i];
            }else{
                a0 = this->alpha[i-1];
                a1 = this->alphalast;
                b0 = this->beta[i-1];
                b1 = this->betalast;
                yy = this->y[i];
                dd = this->dur[i];
            }
            double abnew[2];
            singlestatesimulateAdaptation(yy, dd, a0,
                                          this->alpha[i], a1, b0, this->beta[i], b1,abnew);
            this->alpha[i] = abnew[0];
            this->beta[i] = abnew[1];
        }
        this->alphalast = this->phi1*this->alpha[this->n - 2] +
        this->sigma1*this->random->normal();
        this->betalast = this->phi2*this->beta[this->n - 2] +
        this->sigma2*this->random->normal();
    }
}

void SV2FDUR::simulateparameters()
{
    SV2F::simulateparameters();
    BayesianGammaDuration *bdur =
    new BayesianGammaDuration(this->dur,this->alpha,
                              this->n,this->lambda,this->gam,this->random);
    bdur->simulateparameters();
    this->lambda = bdur->getlambda();
    this->gam = bdur->getbeta();
    delete bdur;
}


double SV2FDUR::getlambda()
{
    return this->lambda;
}

double SV2FDUR::getgam()
{
    return this->gam;
}
