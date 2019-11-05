#ifndef Stats_hpp
#define Stats_hpp


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>



class Stats{
    double *x;
    int n;
    
public:
    
    Stats(){}
    
    ~Stats()
    {
    }
    
    Stats(double *x,int nobs){
        this->n = nobs;
        this->x = x;
    }
    
    double sum()
    {
        double result = 0.0;
        for(int i=0;i<n;i++){
            result += x[i];
        }
        return result;
    }
    
    
    double sumsq()
    {
        double result = 0.0;
        for(int i=0;i<n;i++){
            result += x[i]*x[i];
        }
        return result;
    }
    
    double sumsqminusmu(double mu)
    {
        double result = 0.0;
        for(int i=0;i<n;i++){
            result += (x[i] - mu)*(x[i]-mu);
        }
        return result;
    }
    
    //calculates the mean
    double mean()
    {
        double mvalue = 0.0;
        for(int i=0;i<n;i++){
            mvalue += x[i];
        }
        return mvalue/(double)n;
    }
    
    //calculates the variance
    double variance()
    {
        double mm = mean();
        double msqvalue = 0.0;
        for(int i=0;i<n;i++){
            msqvalue += x[i]*x[i];
        }
        msqvalue /= (double)n;
        return msqvalue - mm*mm;
    }
    
    double kurtosis()
    {
        double mm = mean();
        double vv = variance();
        double value = 0.0;
        for(int i=0;i<n;i++){
            double tt = (x[i] - mm);
            value += tt*tt*tt*tt;
        }
        value /= (double)n;
        return value/(vv*vv);
    }
    
};



#endif
