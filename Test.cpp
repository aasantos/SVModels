//
//  Test.cpp
//  Nov20
//
//  Created by António Alberto Santos on 20/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "Test.hpp"

void estimateSVDUR()
{
    int n;
    double *y = readArrayDouble("y5aapl.txt", &n);
    double *dur = readArrayDouble("d5aapl.txt", &n);
    
    double sigmav = 0.1;
    double phi1 = 0.95;
    double sigma1 = 0.2;
    double lambda = 5.0;
    int niter = 10000;
    int nwarmup = 10000;
    int thin = 5;
    //
    Random *random = new Random(720);
    SVDUR *obj = new SVDUR(y,n,sigmav,phi1,sigma1,random,dur,lambda);
    obj->simulatestatesadaptation(250);
    
    for(int i=0;i<nwarmup;i++){
        if(i%100 == 0) printf("Warmup iteration: %d/%d\n",i,nwarmup);
        obj->simulatestates();
        obj->simulateparameters();
        printf("%.4f %.4f %.4f %.4f\n",
        obj->getsigmav(),obj->getphi1(),obj->getsigma1(),obj->getlambda());
    }
    //
    FILE *fp;
    fp = fopen("svdur5aapl.txt","wa");
    fprintf(fp,"sigmav phi1 sigma1 lambda\n");
    for(int i=0;i<niter;i++){
        if(i%100 == 0) printf("Main iteration: %d/%d\n",i,niter);
        for(int k=0;k<thin;k++){
            obj->simulatestates();
            obj->simulateparameters();
        }
        fprintf(fp,"%.4f %.4f %.4f %.4f\n",
                obj->getsigmav(),obj->getphi1(),obj->getsigma1(),obj->getlambda());
        }
    delete obj;
    fclose(fp);
    
    delete random;
    free(y);
    free(dur);
}

void estimateSV2FDUR()
{
    int n;
    double *y = readArrayDouble("y5aapl.txt", &n);
    double *dur = readArrayDouble("d5aapl.txt", &n);
    //
    double sigmav = 0.1;
    double phi1 = 0.98;
    double sigma1 = 0.2;
    double phi2 = 0.9;
    double sigma2 = 0.2;
    double lambda = 5.0;
    //
    int niter = 10000;
    int nwarmup = 10000;
    int thin = 5;
    //
    Random *random = new Random(720);
    SV2FDUR *obj = new SV2FDUR(y,n,sigmav,phi1,sigma1,random,phi2,sigma2,dur,lambda);
    obj->simulatestatesadaptation(250);
    
    for(int i=0;i<nwarmup;i++){
        if(i%100 == 0) printf("Warmup iteration: %d/%d\n",i,nwarmup);
        obj->simulatestates();
        obj->simulateparameters();
        printf("%.4f %.4f %.4f %.4f %.4f %.4f\n",
        obj->getsigmav(),obj->getphi1(),obj->getsigma1(),
               obj->getphi2(),obj->getsigma2(),obj->getlambda());
    }
    //
    FILE *fp;
    fp = fopen("sv2fdur5aapl.txt","wa");
    fprintf(fp,"sigmav phi1 sigma1 phi2 sigma2 lambda\n");
    for(int i=0;i<niter;i++){
        if(i%100 == 0) printf("Main iteration: %d/%d\n",i,niter);
        for(int k=0;k<thin;k++){
            obj->simulatestates();
            obj->simulateparameters();
        }
         fprintf(fp,"%.4f %.4f %.4f %.4f %.4f %.4f\n",
               obj->getsigmav(),obj->getphi1(),obj->getsigma1(),
                      obj->getphi2(),obj->getsigma2(),obj->getlambda());
        }
    delete obj;
    fclose(fp);
    
    delete random;
    free(y);
    free(dur);
}

