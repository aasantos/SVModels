//
//  func.cpp
//  Nov20
//
//  Created by António Alberto Santos on 20/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#include "func.hpp"

/*
double digamma(double x)
{
    return gsl_sf_psi(x);
}

double trigamma(double x)
{
    return gsl_sf_psi_1(x);
}
*/

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



float* readArrayFloat(const char* file,int *n)
{
    FILE *fp;
    fp = fopen(file, "r");
    if(!fp){
        printf("File not found .... \n");
        return NULL;
    }
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
    if(!fp){
        printf("File not found .... \n");
        return NULL;
    }
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


