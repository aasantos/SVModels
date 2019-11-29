#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "func.hpp"
#include "Random.hpp"
#include "Stats.hpp"
#include "NormalBayesian.hpp"
#include "NormalMixture.hpp"
#include "RegModel.hpp"
#include "AR1Model.hpp"
#include "AR1JModel.hpp"
#include "DLM.hpp"
#include "DLM2F.hpp"
#include "DLM2FB.hpp"
//
#include "SV.hpp"
#include "SVJ.hpp"
#include "SV2F.hpp"
#include "SVJ2F.hpp"
//
#include "SVDUR.hpp"
#include "SVJDUR.hpp"
#include "SV2FDUR.hpp"
#include "SVJ2FDUR.hpp"
//
#include "GammaDuration.hpp"


int main(int argc,const char *argv[])
{
    int n;
    double *y = readArrayDouble(argv[1], &n);
    //
    double sigmav = sqrt(Stats(y,n).variance());
    double phi1 = 0.95;
    double sigma1 = 0.1;
    double phi2 = 0.9;
    double sigma2 = 0.2;
    double muj = 0.5;
    double sigmaj = 1.0;
    double theta = 0.01;
    //
    int niter = atoi(argv[2]);
    int nwarmup = atoi(argv[3]);
    int thin = atoi(argv[4]);
    //
    Random *random = new Random(720);
    SVJ2F *obj = new SVJ2F(y,n,sigmav,phi1,sigma1,random,phi2,sigma2,muj,sigmaj,theta);
    obj->simulatestatesadaptation(250);
    //
    for(int i=0;i<nwarmup;i++){
        if(i%100 == 0) printf("Warmup iteration: %d/%d\n",i,nwarmup);
        obj->simulatestates();
        obj->simulateparameters();
    }
    //
    FILE *fp;
    fp = fopen(argv[5],"wa");
    fprintf(fp,"sigmav phi1 sigma1 phi2 sigma2 muj sigmaj theta\n");
    for(int i=0;i<niter;i++){
        if(i%100 == 0) printf("Main iteration: %d/%d\n",i,niter);
        for(int k=0;k<thin;k++){
            obj->simulatestates();
            obj->simulateparameters();
        }
        fprintf(fp,"%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
        obj->getsigmav(),obj->getphi1(),obj->getsigma1(),obj->getphi2(),obj->getsigma2(),obj->getmuj(),obj->getsigmaj(),obj->gettheta());
        }
    delete obj;
    fclose(fp);
    //
    delete random;
    free(y);
}
