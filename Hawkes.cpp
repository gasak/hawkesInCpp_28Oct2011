//
//  Hawkes.cpp
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/12/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#include "Hawkes.h"
#include "math.h"
#include "stdio.h"
#include <iostream>
#include "MyVector.h"

using namespace std;



/*********************************************************************************************/

/* Return the following functions that are specific to the Hessian model.
 *
 * f = negative log likelihood function. the function that has to be minimized.
 * G = gradient of f w.r.t. parameters
 * H = Hessian of f w.r.t parameters
 *
 */

Matrix<double> Hawkes::getInitialParams(SubSpikeStruct& ss){
	int numTypesA = ss.getNumTypes(); int paramDimension = 2*numTypesA+1;
    Matrix<double> initParams(numTypesA,paramDimension);
    
	
    for(int type = 0; type < numTypesA; type++){
        // lambdaBar_T,m
        int sizeOfType = ss.individualSizes[type];
        double lambdaBar = 0.0;
        if(ss.individualSizes[type] >= 2){
            double totalTime = (ss.individualLists[type])[sizeOfType-1]->getTime() - (ss.individualLists[type])[0]->getTime();
            lambdaBar = ((double) sizeOfType-1)/((double) totalTime);
        }
        else if(ss.individualSizes[type] == 1){
            double totalTime = (ss.individualLists[type])[sizeOfType-1]->getTime();
            lambdaBar = 1/((double) totalTime);
        }
        else{
            lambdaBar = 0.0;
        }
        initParams.Set_Element(type, 0, lambdaBar);
    
        // gamma's
        for(int i = 0; i < ss.numTypes ; i++){
            // compute gamma_m',m
            initParams.Set_Element(type, i+1, 0.0);
        }
    
        // alpha's
        for(int i = 0; i < ss.numTypes ; i++){
            // compute alpha_m',m
            initParams.Set_Element(type, numTypesA+i+1, 0.0);
        }
    }
    return initParams;
}

int Hawkes::typeRequiringUpdateExcept(SubSpikeStruct &ss, Matrix<double> paramMatrix, MyVector<int> ignoreTheseTypes){
    MyVector<int> typesToConsider(ss.getNumTypes()-ignoreTheseTypes.GetLength());
    MyVector<int> allTypes(ss.getNumTypes());
    for(int i = 0; i < allTypes.length; i++){
        allTypes.Set_Element(i, i);
    }
    typesToConsider = allTypes.AnotB(ignoreTheseTypes);
    //typesToConsider.GetTransposed().print();
    
    MyVector<double> rho(typesToConsider.GetLength());
    double T = ss.entireList[ss.getSize()-1]->getTime();
    for(int j=0; j < typesToConsider.GetLength(); j++){
        int temp = typesToConsider.Get_Element(j);
        double lambda = paramMatrix.Get_Element(temp, 0);
        int tempsize = ss.individualSizes[temp];
        if (tempsize >= 2) {
            double tlast = (ss.individualLists[temp])[tempsize-1]->getTime();
            rho.Set_Element(j, exp(lambda*(T-tlast)));
        }
        else{
            rho.Set_Element(j, 0);
        }
    }
    double sum = rho.sum();
    rho = rho*(1/sum);
    //cout << "printing rho" << endl;
    //rho.GetTransposed().print();
    double val = (double)rand()/(double)RAND_MAX; //cout << "Vals is " << val << endl;
    int flag = 1, ret = -1, j=0; double cumsum = 0;
    do {
        cumsum = cumsum+rho.Get_Element(j);
        if( val <= cumsum){
            ret = j; flag = 0; 
        }
        j++;
    } while (flag == 1);
    ret = typesToConsider.Get_Element(ret);
    
    /*
    // choose highest rho.
    int highestElem = -1; double highVal = 0;
    for(int i = 0 ; i < rho.GetLength(); i++){
        if(rho.Get_Element(i) > highVal){
            highVal = rho.Get_Element(i);
            highestElem = i;
        }
    }
    ret = typesToConsider.Get_Element(highestElem);
    */
    
    return ret;
}

double Hawkes::computeNegativeLogLikelihood(SubSpikeStruct &ss, Matrix<double> paramMatrixptr, int dimension, int trace){
    double alpha = 1e+9; // penalty
    double logLik;
    
    // compute Cond Int at spikes.
    vector<int> sizes = ss.individualSizes; int size = ss.getSize();
    int spikes_size = sizes[dimension], numTypes = ss.getNumTypes();
    
    // are all Params Positive?
    bool condition_allparamspositive = true;
    for(int iter = 0 ; iter < (2*numTypes+1); iter++){
        if(paramMatrixptr.Get_Element(dimension, iter) < 0){
            condition_allparamspositive = false;
        }
    }
    if(!condition_allparamspositive){
        if(trace) cout << "Negative Parameter - Penalty added.\n";
        logLik = -alpha;
    }
    else{
    // f = sum(log(condInt.At.Spikes)) - int_{t_0}^{t_1} cond.Int(t) dt
    
    
        MyVector<double> condIntensities(spikes_size);
        MyVector<double> timingsForEvaluation = ss.getTimingsVector_ForType(dimension);
        for(int j = 0; j < spikes_size; j++){
            double lambda_m = paramMatrixptr(dimension,0);
            double condint = lambda_m;
            for(int i=0; i < numTypes; i++){
                double alpha_mprime = paramMatrixptr(dimension,numTypes+i+1); //cout << "alpha" << alpha_mprime <<endl;
                double gamma_mprime = paramMatrixptr(dimension,i+1); //cout << "gamma" << gamma_mprime << endl;
                MyVector<double> timingsMprime = ss.getTimingsVector_ForType(i);
                
                MyVector<double> timingsLessThanEvalTime = timingsMprime[timingsMprime < timingsForEvaluation(j)];
                const double tval = timingsForEvaluation(j);
                //cout << "timing for eval = " << timingsForEvaluation(j) << " and timingsLessThanEvalTime = " << endl;
                //timingsLessThanEvalTime.GetTransposed().print();
                //(timingsLessThanEvalTime-tval).GetTransposed().print();
                Matrix<double> temp1 =(timingsLessThanEvalTime-tval)*alpha_mprime + log(gamma_mprime);
                //cout << "temp1 is " << endl;
                //temp1.GetTransposed().print();
                MyVector<double> temp(temp1);
                condint += temp.expo().sum();
                //cout << temp.expo().sum() << endl;

            }
            condIntensities.Set_Element(j, condint);
        }



        // compute Integral First To Last A
        const double t0 = ss.entireList[0]->getTime(), t1=ss.entireList[size-1]->getTime();
        double lambda_m = paramMatrixptr(dimension,0);
        double integral_condint = lambda_m*(t1-t0);

        for(int i=0; i < numTypes; i++){
            double alpha_mprime = paramMatrixptr(dimension,numTypes+i+1);
            double gamma_mprime = paramMatrixptr(dimension,i+1);
            MyVector<double> timingsMprime = ss.getTimingsVector_ForType(i);
            
            MyVector<double> timingsLessThanEvalTime = timingsMprime[timingsMprime < t1];
            MyVector<double> temp((timingsLessThanEvalTime-t1)*alpha_mprime);
            integral_condint += (sizes[i]-(temp.expo().sum()))*gamma_mprime/alpha_mprime;
        }
        logLik = condIntensities.loga().sum() - integral_condint;
        if(trace == 1){
            cout << "conditional Intensities at spikes are " << endl;
            condIntensities.GetTransposed().print();
            cout << integral_condint << " is the integral" << endl;
            cout << "negative LogLikelihood is " << -1*logLik << endl;
        }

        double barrier = 0; 
        double mu = 0.5;
        for(int param = 0 ; param < (2*numTypes)+1; param++){
            barrier += -10* log(mu*paramMatrixptr.Get_Element(dimension, param));
        }
    }
    return -1*logLik;
    
}



/************************************** ADMIN STUFF - MEMORY MANAGEMENT *************************************************/
/*
 *  init(int dim)
 *  helper file for constructor
 */
void Hawkes::init(int dim){
	this->numTypes = dim; // M = dim
    paramMatrix.Set_Size(dim, 2*dim+1);
}

/*
 *  constructer for hawkes
 *  initiates - dimensions, param matrix.
 */

Hawkes::Hawkes(int dim){
	init(dim);
}

Hawkes::Hawkes(){
    
}

Hawkes::Hawkes(const Hawkes & copyObj){
	copy(copyObj);
}

const Hawkes & Hawkes::operator=(const Hawkes & rhs){
	if(this != &rhs)
	{
		clear();
		copy(rhs);
	}
	return *this;
}

/*
 * destructor
 */

Hawkes::~Hawkes() {
	clear();
}

void Hawkes::copy(const Hawkes & copyObj){
	this->s = copyObj.s;
	this->paramMatrix = copyObj.paramMatrix;
	this->currentDimension = copyObj.currentDimension;
	this->numTypes = copyObj.numTypes;
    this->trace = copyObj.trace;
}

void Hawkes::clear(){
}
