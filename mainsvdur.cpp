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
    double *dur = readArrayDouble(argv[2],&n);
    
    double sigmav = sqrt(Stats(y,n).variance());
    double phi1 = 0.95;
    double sigma1 = 0.2;
    double lambda = Stats(dur,n).mean() + 1.0;
    //
    int niter = atoi(argv[3]);
    int nwarmup = atoi(argv[4]);
    int thin = atoi(argv[5]);
    //
    Random *random = new Random(720);
    SVDUR *obj = new SVDUR(y,n,sigmav,phi1,sigma1,random,dur,lambda);
    obj->simulatestatesadaptation(250);
    
    for(int i=0;i<nwarmup;i++){
        if(i%100 == 0) printf("Warmup iteration: %d/%d\n",i,nwarmup);
        obj->simulatestates();
        obj->simulateparameters();
    }
    //
    FILE *fp;
    fp = fopen(argv[6],"wa");
    fprintf(fp,"sigmav phi1 sigma1 lambda\n");
    for(int i=0;i<niter;i++){
        if(i%100 == 0) printf("Main iteration: %d/%d\n",i,niter);
        for(int k=0;k<thin;k++){
            obj->simulatestates();
            obj->simulateparameters();
        }
        fprintf(fp,"%.4f %.4f %.4f %4f\n",
                obj->getsigmav(),obj->getphi1(),obj->getsigma1(),obj->getlambda());
        }
    delete obj;
    fclose(fp);
    
    delete random;
    free(y);
    free(dur);
	return 0;
}
