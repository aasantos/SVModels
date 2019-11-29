//
//  GammaDuration.hpp
//  Nov20
//
//  Created by António Alberto Santos on 20/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef GammaDuration_hpp
#define GammaDuration_hpp

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "func.hpp"
#include "Random.hpp"
#include "Stats.hpp"

class GammaDuration{
    
protected:
    int n;
    double *dur;
    double *alpha;
    double lambda;
    double lambdaprior[2];
    //
    double lambdainf;
    double lambdasup;
    Random *random;

public:
    GammaDuration(double *dur,double *alpha,int n,double lambda,Random *random)
    {
        this->dur  = dur;
        this->alpha = alpha;
        this->n = n;
        this->lambda = lambda;
        this->lambdaprior[0] = Stats(this->dur, this->n).mean();
        this->lambdaprior[1] = 1.0;
        this->lambdainf = 1.0;
        this->lambdasup = 100.0;
        this->random = random;
    }
    
    double lognormalpdf(double x,double m,double s);
    double loglikdur(double ll);
    double loglikdurprior(double ll,double a,double b);
    double metroproblambda(double llnew,double ll,double a,double b,double mm,double ss);
    void lambdameanvar(double l0,double out[2]);
    double dfa(double x);
    double ddfa(double x);
    double newtownlambda(double l0);
    void simulatelambda();
    double getlambda();
    void setlambdalim(double ll[2]);
    void setlambdaprior(double ll[2]);
    
};


#endif /* GammaDuration_hpp */
