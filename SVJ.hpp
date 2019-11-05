//
//  SVJ.hpp
//  November2019
//
//  Created by António Alberto Santos on 03/11/2019.
//  Copyright © 2019 António Alberto Santos. All rights reserved.
//

#ifndef SVJ_hpp
#define SVJ_hpp

#include <stdio.h>
#include "SV.hpp"

class SVJ: public SV{
protected:
    double muj;
    double sigmaj;
    double theta;
    double *Z;
    int *J;

public:
    SVJ(double *y,double sigmav,double phi1,double sigma1,Random *random,double muj,double sigmaj,double theta):SV(y,n,sigmav,phi1,sigma1,random)
    {
        this->muj = muj;
        this->sigmaj = sigmaj;
        this->theta = theta;
        this->Z = (double*)calloc(n,sizeof(double));
        this->J = (int*)calloc(n,sizeof(int));
    }
    
    ~SVJ()
    {
        if(Z){
            free(Z);
        }
        if(J){
            free(J);
        }
    }
    
    
};

#endif /* SVJ_hpp */
