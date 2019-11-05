#ifndef Random_hpp
#define Random_hpp


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define PI  3.141592653589793

#define  A1  (-3.969683028665376e+01)
#define  A2   2.209460984245205e+02
#define  A3  (-2.759285104469687e+02)
#define  A4   1.383577518672690e+02
#define  A5  (-3.066479806614716e+01)
#define  A6   2.506628277459239e+00

#define  B1  (-5.447609879822406e+01)
#define  B2   1.615858368580409e+02
#define  B3  (-1.556989798598866e+02)
#define  B4   6.680131188771972e+01
#define  B5  (-1.328068155288572e+01)

#define  C1  (-7.784894002430293e-03)
#define  C2  (-3.223964580411365e-01)
#define  C3  (-2.400758277161838e+00)
#define  C4  (-2.549732539343734e+00)
#define  C5   4.374664141464968e+00
#define  C6   2.938163982698783e+00

#define  D1   7.784695709041462e-03
#define  D2   3.224671290700398e-01
#define  D3   2.445134137142996e+00
#define  D4   3.754408661907416e+00

#define P_LOW   0.02425
/* P_high = 1 - p_low*/
#define P_HIGH  0.97575




//interface for an object used to generate random numbers
class Random{
public:
    
    Random(){
        
        this->x = 3023082353252;
        this->y = 10970392323412324;
        this->z1 = 2370343136432;
        this->c1 = 2429979799696;
        this->z2 = 74624524523452;
        this->c2 = 436452452245736;
    }
    
    
    Random(int idxx)
    {
        unsigned long long idx = (unsigned long long)idxx;
        unsigned long long seed[6] =
        {(idx + 546)*79794,(idx+890)*564534,
            (idx+98)*580345,(idx+9)*45,
            (idx+23)*897,(idx+ 67)*212879987};
        this->x = seed[0];
        this->y = seed[1];
        this->z1 = seed[2];
        this->c1 = seed[3];
        this->z2 = seed[4];
        this->c2 = seed[5];
    }
    
     ~Random(){};
    
    unsigned long long JLKISS64()
    {
        unsigned long long t = 0;
        x = 1490024343005336237ULL * x + 123456789;
        y ^= y << 21; y ^= y >> 17; y ^= y << 30;
        /* Do not set y=0! */
        t = 4294584393ULL * z1 + c1;
        c1 = t >> 32;
        z1 = t;
        t = 4246477509ULL * z2 + c2;
        c2 = t >> 32;
        z2 = t;
        return x + y + z1 + ((unsigned long long)z2 << 32);
        /* Return 64-bit result */
    }
    
    
    double uniform()
    {
        unsigned long long div = 18446744073709551615ULL; //2^64 -1
        return (double)(JLKISS64()/(1.0*div));
    }
    
    double normal()
    {
        double p =  uniform();
        double x = 0.0;
        double q = 0.0;
        double r = 0.0;
        if ((0 < p )  && (p < P_LOW)){
            q = sqrt(-2*log(p));
            x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
        }
        else{
            if ((P_LOW <= p) && (p <= P_HIGH)){
                q = p - 0.5;
                r = q*q;
                x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
            }
            else{
                if ((P_HIGH < p)&&(p < 1)){
                    q = sqrt(-2*log(1-p));
                    x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
                }
            }
        }
        return x;
    }
    
    
    double exponential(double lambda)
    {
        double u = uniform();
        return -1.0*logf(u)/lambda;
    }
    
    double randGammabest(double a)
    {
        int flag = 0;
        double b = 0.0;
        double c = 0.0;
        double u = 0.0;
        double v = 0.0;
        double w = 0.0;
        double y = 0.0;
        double x = 0.0;
        double z = 0.0;
        double test1 = 0.0;
        double test2 = 0.0;
        double test3 = 0.0;
        double temp = 0.0;
        b = a - 1.0;
        c = (12.0*a - 3.0)/4.0;
        flag = 0;
        while(flag == 0){
            u = uniform();
            v = uniform();
            w = u*(1.0-u);
            y = sqrtf(c/w)*(u-0.5);
            x = b + y;
            if(x>0){
                z = 64.0*v*v*w*w*w;
                test1 = 1.0 - 2*y*y/x;
                test2 = 2.0*(b*logf(x/b)-y);
                test3 = logf(z);
                if(z <= test1 || test2 >= test3){
                    temp = x;
                    flag = 1;
                }
            }
        }
        return temp;
    }
    
    double randGammaless(double a)
    {
        double temp = 0.0;
        double x = 0.0;
        double u = uniform();
        temp = randGammabest((a + 1.0));
        x = temp*powf(u,1.0/a);
        return x;
    }
    
    double gamma(double a,double b)
    {
        double g = 0.0;
        if(a<1.0){
            g = (1/b)*randGammaless(a);
        }
        else if(a==1.0){
            g = (1/b)*exponential(a);
        }
        else{
            g = (1/b)*randGammabest(a);
        }
        return g;
    }
    
    double beta(double a,double b)
    {
        double x = gamma(a,1.0);
        double y = gamma(b,1.0);
        double result = x/(x + y);
        return result;
    }
    
    double t(int df)
    {
        double u = normal();
        double v = gamma(0.5*(double)df,0.5);
        return u/sqrt(v/(double)df);
    }
    
    int bernoulli(double p)
    {
        if(uniform() < p){
            return 1;
        }else{
            return 0;
        }
    }
    
    
    void bivariatenormal(double m[2],double S[4],double x[2])
    {
        double a = S[0];
        double b = S[1];
        double c = S[3];
        //
        double t[2];
        t[0] = normal();
        t[1] = normal();
        x[0] = sqrt(a)*t[0];
        double tt = sqrt((a*c - b*b)/a);
        //stupid error stupid stupid
        //x[1] = (b/sqrt(a))*t[0] + sqrt(tt)*t[1];
        x[1] = (b/sqrt(a))*t[0] + tt*t[1];
        x[0] += m[0];
        x[1] += m[1];
    }
    
    
private:
    unsigned long long x;
    unsigned long long y;
    unsigned long long z1;
    unsigned long long c1;
    unsigned long long z2;
    unsigned long long c2;
};


#endif
