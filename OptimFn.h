//
//  OptimFn.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/28/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#ifndef hawkesInCpp_OptimFn_h
#define hawkesInCpp_OptimFn_h
class OptimFn{
public:
    virtual double getFunVal(int n, double* Bvec){
        cout << "Entered Base Class \n";
        return 0;
    }
};

class OptimFn_Parabola : public OptimFn{
public:
    double getFunVal(int n, double* Bvec){
        cout << "Entered Parabola Function \n";
        return ((Bvec[0]-3)*(Bvec[0]-3) + (Bvec[1]-6)*(Bvec[1]-6))+10000;
    };
};

class OptimFn{
public:
    virtual double getFunVal(int n, double* Bvec){
        cout << "Entered Base Class \n";
        return 0;
    }
};

class OptimFn_Parabola : public OptimFn{
public:
    double getFunVal(int n, double* Bvec){
        cout << "Entered Parabola Function \n";
        return ((Bvec[0]-3)*(Bvec[0]-3) + (Bvec[1]-6)*(Bvec[1]-6))+10000;
    };
};

#endif
