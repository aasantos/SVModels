#ifndef BayesianGammaDuration_hpp
#define BayesianGammaDuration_hpp


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "Utils.hpp"
#include "Random.hpp"


class BayesianGammaDuration{
//
protected:
	int n;
	double *dur;
	double *alpha;
	double lambda;
	double beta;
	//
	Random *random;

public:
	BayesianGammaDuration(double *dur,double *alpha,int n,double lambda,double beta,Random *random)
	{
		this->dur  = dur;
		this->alpha = alpha;
		this->n = n;
		this->lambda = lambda;
		this->beta = beta;
		this->random = random;
		//
	}


    double lognormalpdf(double x,double m,double s);
    void simulatebeta();
    double dfa(double x);
    double ddfa(double x);
    double newtownlambda(double l0);
    void lambdameanvar(double l0,double out[2]);
    void simulatelambda();
    void simulateparameters();
    double getlambda();
    double getbeta();
};


#endif
