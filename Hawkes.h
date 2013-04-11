//
//  Hawkes.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/12/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#ifndef hawkesInCpp_Hawkes_h
#define hawkesInCpp_Hawkes_h

#include "SubSpikeStruct.h"
#include "Matrix.h"
#include "OptimFn.h"

class Hawkes : public OptimFn {
    
    
public:
	SubSpikeStruct s;
	Matrix<double> paramMatrix;
	int currentDimension;
    int trace;
    int numTypes;  // M
    
public:
    // constructor
	Hawkes();
	Hawkes(int numTypes);
    // big three
	Hawkes(const Hawkes & src);
	const Hawkes & operator=(const Hawkes & rhs);
	virtual ~Hawkes();
    
    // functionalities
	double computeNegativeLogLikelihood(SubSpikeStruct & ss, Matrix<double> paramMatrixptr, int dimension, int trace);
    Matrix<double> getInitialParams(SubSpikeStruct& ss);
    int typeRequiringUpdateExcept(SubSpikeStruct& ss, Matrix<double> paramMatrix, MyVector<int> ignoreTheseTypes);
    
	// setters
	void setData(SubSpikeStruct s){
		this->s = s;
	}
	void setParams(Matrix<double> paramMatrix){
		this->paramMatrix = paramMatrix;
	}
	void setCurrentDimension(int dimension){
		this->currentDimension = dimension;
	}
    void setTrace(int trace){
        this->trace = trace;
    }
    
    //optimization
    double getFunVal(int n, double* Bvec){
        if(trace) cout << "Entered Hawkes Function \n";
        for(int i = 0; i < n; i++){
            paramMatrix.Set_Element(currentDimension, i, Bvec[i]);
        }
        
        return computeNegativeLogLikelihood(s, paramMatrix, currentDimension,trace)+10000;
    };
    

    
private:
	// admin helpers
	void init(int dim);
    void copy(const Hawkes & copyObj);
	void clear();
    
};

#endif
