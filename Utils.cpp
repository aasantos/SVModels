//
//  Utils.cpp
//  October2019
//
//  Created by António Alberto Santos on 15/10/2019.
//  Copyright © 2019 Antonio Santos. All rights reserved.
//

#include "Utils.hpp"

float* readArrayFloat(const char* file,int *n)
{
    FILE *fp;
    fp = fopen(file, "r");
    char buff[64];
    int nrow = -1;
    while (!feof(fp)) {
        fscanf(fp, "%s",buff);
        nrow++;
    }
    *n = nrow;
    rewind(fp);
    float *result = (float*)malloc(nrow*sizeof(float));
    for(int i=0;i<nrow;i++){
        fscanf(fp, "%s",buff);
        result[i] = atof(buff);
    }
    fclose(fp);
    return result;
}


double* readArrayDouble(const char* file,int *n)
{
    FILE *fp;
    fp = fopen(file, "r");
    char buff[64];
    int nrow = -1;
    while (!feof(fp)) {
        fscanf(fp, "%s",buff);
        nrow++;
    }
    *n = nrow;
    rewind(fp);
    double *result = (double*)malloc(nrow*sizeof(double));
    for(int i=0;i<nrow;i++){
        fscanf(fp, "%s",buff);
        result[i] = atof(buff);
    }
    fclose(fp);
    return result;
}


void writeArrayDouble(double *x,const char* file,int n)
{
    FILE *fp;
    fp = fopen(file, "wa");
    for(int i=0;i<n;i++) fprintf(fp,"%.6f\n",x[i]);
    fclose(fp);
}

double digamma(double x)
{
    double r, f, t;
    r = 0;
    while (x<=5)
    { r -= 1/x;
        x += 1;
    }
    f = 1/(x*x);
    t = f*(-1/12.0 + f*(1/120.0 +
        f*(-1/252.0 + f*(1/240.0 + f*(-1/132.0
        + f*(691/32760.0 + f*(-1/12.0 + f*3617/8160.0)))))));
    return r + log(x) - 0.5/x + t;
}

double trigamma(double x)
{
    double h = 0.000001;
    double f1 = digamma(x + h);
    double f0 = digamma(x - h);
    return (f1 - f0)/(2.0*h);
}


double norm(double x[2])
{
    return sqrt(x[0]*x[0] + x[1]*x[1]);
}

/*
double digamma(double x)
{
    double xx = x;
    if(xx < 0.000001) xx = 0.000001;
    return gsl_sf_psi(xx);
}


double trigamma(double x)
{
    double xx = x;
    if(xx < 0.000001) xx = 0.000001;
    return gsl_sf_psi_1(xx);
}
*/

 
 
double logpdfbivariatenormal(double x[2],double mm[2],double vcv[4])
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


