//
//  main.cpp
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#include <iostream>
#include "Matrix.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <string.h>
#include <list>
#include "stdio.h"
#include <stdlib.h>
#include "spikeTrain.h"
#include "SubSpikeStruct.h"
#include "csvManagement.h"
#include "Matrix.h"
#include "Hawkes.h"
#include "MyVector.h"
#include "limits"
#include "OptimFn.h"
#include "Optimizer.h"

void testMatrixFile(){
    Matrix<int> a;
    a.Set_Size(2,3);
    a.Set_Element(1, 1, 2);
    cout << a.Get_Element(1, 1);
}

void testVector(){
    MyVector<double> x(5);
    x.Set_Element(0,1.25);
    x.Set_Element(1,2);
    x.Set_Element(2,3);
    x.Set_Element(3,4);
    x.Set_Element(4,5);
    
    x.print();
    cout << "check substraction." << endl;
    (x - 3).print();
    x.print();
    
    cout << "check multiplication." << endl;
    ((x-3)*2).print();
    cout << exp(x(0)) << endl;
    cout << "check exponent." << endl;
    MyVector<double> z = x.expo();
    z.print();
    cout << "check sum" << endl;
    cout << z.sum() << endl;
    cout << "check less than equal to" << endl;
    MyVector<bool> w = (x<= 3);
    (x[w]).print();
    
}

void testSpikeTrain2(){
    
	csvManagement c;    
	const char * fileName = "/Users/sivaiitm/Documents/xcodeworkspace/hawkesInCpp/hawkesInCpp/aapl_cpp_test.csv";
    
    // read data from the file <fileName>. 
	SubSpikeStruct ss = c.csvToSpikeTrain(fileName);
    
	//ss.printSpikeTrain();
    //ss.getSpikesBetween(0, 30).printSpikeTrain();
    ss.printSummary();
    
    //vector<string> wantedList; wantedList.push_back("A"); wantedList.push_back("B");
    //SubSpikeStruct sss = ss.filterEventsOfType(wantedList);
    //sss.printSummary();
    
    cout << "*****************************\n";
    cout << "Testing Model - Hawkes\n";
    
    int numTypes = ss.getNumTypes(); 
    
    Hawkes h(ss.getNumTypes());
	Matrix<double> paramMatrixPtr1M(numTypes,(2*numTypes+1));
    
	for(int i=0; i<numTypes; i++){
		for(int j = 0; j<(2*numTypes+1); j++){
			if(j == 0){
				paramMatrixPtr1M.Set_Element(i,j,1.0);
			}
			else{
				paramMatrixPtr1M.Set_Element(i,j,0.2);
			}
		}
	}
    cout << h.computeNegativeLogLikelihood(ss, paramMatrixPtr1M, 1,1) << endl;
    
    cout << "**************************************************************************\n";
    cout << "Start Optimization - Single Step - MLE Over Entire DataSet\n";
    
    cout << "Initial Parameter Matrix - (lambda,gamma_{i,1},gamma_{i,2},alpha_{i,1},alpha_{i,2}) = " << endl;
    paramMatrixPtr1M.print();
    cout << "Initial NegLogLikelihood = (";
    for(int typeiter = 0; typeiter < ss.getNumTypes(); typeiter++){
        cout << h.computeNegativeLogLikelihood(ss, paramMatrixPtr1M, typeiter, 0) << "," ;
    }
    cout << ")" << endl;
    
    Optimizer o;
    // optimization parameters
    int n = 2*ss.getNumTypes()+1;
    double abstol = -INFINITY; // for non-neg functions. general function, put -Inf
    double intol = 1e-8;  // if fun reduces by less than eps, stop.
    double alpha = 1;
    double bet = 0.5;
    double gamm = 2;
    int trace = 0;
    int* fncount = new int;
    int maxit = 5000; //max iterations
    
    for(int currentDimension = 0; currentDimension < ss.getNumTypes(); currentDimension++){
        cout << "Optimizing over Type " << ss.eventTypes[currentDimension] << endl;
        //int currentDimension = 1;
        h.setData(ss);
        h.setParams(paramMatrixPtr1M);
        h.setCurrentDimension(currentDimension);
        h.setTrace(0); //Don't trace
        
        OptimFn * fptr = &h;
        double* X; X = new double[n]; //?
        double* Fmin = new double; //final min value
        int* fail = new int;    // minimization success
        double* Bvec; Bvec = new double[n]; // initial parameters
        for(int i = 0; i < n; i++){
            Bvec[i] = paramMatrixPtr1M.Get_Element(currentDimension, i);
        }
        
        //cout << "n is " << n << " Bvec[0] = " << Bvec[0] << " Bvec[1] = " << Bvec[1] << endl;
        o.nmmin(n, Bvec,X, Fmin, fptr, fail, abstol, intol, 
                alpha,  bet,  gamm,  trace,
                fncount,  maxit);
        
        for(int i = 0; i < n; i++){
            paramMatrixPtr1M.Set_Element(currentDimension, i, X[i]);
        }
        cout << *fncount << " number of iterations to coverge. Function Minimum is " << *Fmin-10000 << " observed at " << endl;
        paramMatrixPtr1M.print();
    }
    
    cout << "Final NegLogLikelihood = (";
    for(int typeiter = 0; typeiter < ss.getNumTypes(); typeiter++){
        cout << h.computeNegativeLogLikelihood(ss, paramMatrixPtr1M, typeiter, 0) << "," ;
    }
    cout << ")" << endl;
    
    cout << "Done." << endl;
    
}

/*
void testHawkesMLEstimateUpdating(){
    const char * fileName = "bigtest.csv";
    
	int numTypes = 12;
    
	string* eventTypesPtr;
	string eventTypes[12] = {"A","B","C","D","J","K","O","P","Q","X","Y","Z"};
	eventTypesPtr = eventTypes;
    
    string* eventTypesPtr1;
	string eventTypes1[2] = {"A","B"};
	eventTypesPtr1 = eventTypes1;
    
    // read data from the file <fileName>. Specify which events are of interest (Other event types filtered out). 
    csvManagement c;
    spikeTrain s = c.getSpikeTrain(fileName,numTypes,eventTypesPtr);
    
	//s.setEventTypesPtr(eventTypesPtr);
	s.createMainStruct();
    SubSpikeStruct ss;
	ss = s.getMainStruct();
    ss.printSummary();
    Hawkes h(s.getNumTypes());
   

    cout << "**************************************************************************\n";
    cout << "Start Optimization\n";
    
    // counters
    MyVector<int> counter_encountered(numTypes);
    MyVector<int> counter_randomlychosen(numTypes);
    MyVector<int> summary_initialEvents(numTypes);
    MyVector<int> summary_endEvents(numTypes);
    MyVector<double> initialLikelihood(numTypes);
    MyVector<double> finalLikelihood(numTypes);
    
    // get first 30 (out of 44) events.
    int windowsize = 300;
    SubSpikeStruct sOpt = ss.getSpikesBetween(0, windowsize-1);
    Matrix<double> paramMatrixPtr1M = h.getInitialParams(sOpt);
    for(int i=0; i<numTypes; i++){
		for(int j = 0; j<(2*numTypes+1); j++){
			if(j != 0){
				paramMatrixPtr1M.Set_Element(i,j,0.2);
			}
		}
	}
    //paramMatrixPtr1M.print();
    
    cout << "Initial Parameter Matrix - (lambda,gamma_{i,1},gamma_{i,2},alpha_{i,1},alpha_{i,2}) = " << endl;
    paramMatrixPtr1M.print();
    cout << "Initial NegLogLikelihood = (";
    for(int typeiter = 0; typeiter < ss.getNumTypes(); typeiter++){
        double tempval = h.computeNegativeLogLikelihood(ss, paramMatrixPtr1M, typeiter, 0);
        cout << tempval << "," ;
        initialLikelihood.Set_Element(typeiter, tempval);
    }
    cout << ")" << endl;
    
    Optimizer o;
    // optimization parameters
    int n = 2*s.getNumTypes()+1;
    double abstol = -INFINITY; // for non-neg functions. general function, put -Inf
    double intol = 1e-8;  // if fun reduces by less than eps, stop.
    double alpha = 1;
    double bet = 0.5;
    double gamm = 2;
    int trace = 0;
    int* fncount = new int;
    int maxit = 5000; //max iterations
    
    for(int movingWindow = 0; movingWindow < 131; movingWindow++){
        cout << "=====window " << movingWindow+0 << "-" << movingWindow+windowsize-1 << endl;
        sOpt = ss.getSpikesBetween(movingWindow+0, movingWindow+windowsize-1);
        
        if(movingWindow == 0){
            for(int tempiter = 0; tempiter < numTypes; tempiter++){
                summary_initialEvents.Set_Element(tempiter, sOpt.individualSizes[tempiter]);
            }
        }

        int currentDimension = sOpt.getLastEventTypeAsInteger();
        cout << "Encountered event is of type " << sOpt.getLastEventType() << "-" << currentDimension << endl; 
        double* X; X = new double[n];
        double* Fmin = new double; //final min value
        int* fail = new int;    // minimization success
        if(movingWindow == 0) 
            maxit = 200;
        else
            maxit = 100; //max iterations
        o.performNMmin(paramMatrixPtr1M, currentDimension, h, sOpt, n, X, abstol, intol, alpha, bet, gamm, trace, fncount, maxit, Fmin, fail );
        
        for(int i = 0; i < n; i++){
            paramMatrixPtr1M.Set_Element(currentDimension, i, X[i]);
        }
        
        counter_encountered.Set_Element(currentDimension, counter_encountered.Get_Element(currentDimension)+1);
        cout << *fncount << " number of iterations to coverge. Function Minimum is " << *Fmin-10000 << " observed at " << endl;
        //paramMatrixPtr1M.print();
            
        MyVector<int> ignoreThisType(1); ignoreThisType.Set_Element(0, currentDimension);
        int alsoUpdate = h.typeRequiringUpdateExcept(sOpt, paramMatrixPtr1M, ignoreThisType);
        cout << "also updating event " << eventTypesPtr[alsoUpdate] << " - " << alsoUpdate << endl;
        
        o.performNMmin(paramMatrixPtr1M, alsoUpdate, h, sOpt, n, X, abstol, intol, alpha, bet, gamm, trace, fncount, maxit, Fmin, fail );
        
        for(int i = 0; i < n; i++){
            paramMatrixPtr1M.Set_Element(alsoUpdate, i, X[i]);
        }
        counter_randomlychosen.Set_Element(alsoUpdate, counter_randomlychosen.Get_Element(alsoUpdate)+1);
        cout << *fncount << " number of iterations to coverge. Function Minimum is " << *Fmin-10000 << " observed at " << endl;
        //paramMatrixPtr1M.print();
             
    }
    
    cout << "Final Parameter Matrix - (lambda,gamma_{i,1},gamma_{i,2},alpha_{i,1},alpha_{i,2}) = " << endl;
    paramMatrixPtr1M.print();
    cout << "Final NegLogLikelihood = (";
    for(int typeiter = 0; typeiter < ss.getNumTypes(); typeiter++){
        double tempval = h.computeNegativeLogLikelihood(ss, paramMatrixPtr1M, typeiter, 0);
        cout << tempval << "," ;
        finalLikelihood.Set_Element(typeiter, tempval);
    }
    cout << ")" << endl;
    
    cout << "Summary of parameter updates:" << endl;
    cout << "Type - InitEvents - EndEvents - Encountered - Randomly Updated - TotalUpdates - Init. Likelihood - FinaL Likelihood" << endl;
    for(int i=0; i < numTypes; i++){
        cout << eventTypesPtr[i] << " - " << summary_initialEvents.Get_Element(i) << " - " << summary_initialEvents.Get_Element(i)+counter_encountered.Get_Element(i) << " - " << counter_encountered.Get_Element(i) << " - " << counter_randomlychosen.Get_Element(i) << " - "<< counter_encountered.Get_Element(i)+counter_randomlychosen.Get_Element(i) << " - " << initialLikelihood.Get_Element(i)<< " - " << finalLikelihood.Get_Element(i) << endl;
    }
    cout << "Total - " << summary_initialEvents.sum() << " - " << summary_initialEvents.sum()+counter_encountered.sum() << " - " <<counter_encountered.sum() << " - " << counter_randomlychosen.sum() << " - "<< counter_encountered.sum()+counter_randomlychosen.sum() << " - "<< initialLikelihood.sum()<< " - " << finalLikelihood.sum() <<endl;
    
    cout << "Done." << endl;
    

    
}
*/
void testNumericLimits(){
    cout << "Minimum value for int: " << numeric_limits<double>::min() << endl;
    cout << "Maximum value for int: " << numeric_limits<double>::max() << endl;
    cout << "int is signed: " << numeric_limits<double>::is_signed << endl;
    cout << "Non-sign bits in int: " << numeric_limits<double>::digits << endl;
    cout << "int has infinity: " << numeric_limits<double>::has_infinity << endl;
}




void testOptimization(){
        
    OptimFn_Parabola fminfn;
    OptimFn * fptr = &fminfn;
    Optimizer o;
    
    int n = 2;
    double* Bvec; Bvec = new double[2]; Bvec[0]=10;Bvec[1]=9; // params
    double* X; X = new double[2]; //?
    double* Fmin = new double; //final min value
    int* fail = new int;    // minimization success
    double abstol = -INFINITY; // ?
    double intol = 1e-8;  // ?
    double alpha = 1;
    double bet = 0.5;
    double gamm = 2;
    int trace = 1;
    int* fncount = new int;
    int maxit = 1000; //max iterations
     cout << "n is " << n << " Bvec[0] = " << Bvec[0] << " Bvec[1] = " << Bvec[1] << endl;
    o.nmmin(n, Bvec,X, Fmin, fptr, fail, abstol, intol, 
                alpha,  bet,  gamm,  trace,
                fncount,  maxit);
    
    for(int i = 0; i < n; i++){
        cout << Bvec[i] << "," << X[i] << endl;
    }
 
}

int main (int argc, const char * argv[])
{
    //testVector();
    //testMatrixFile();
    //testNumericLimits();
    //testSpikeTrain();
    testSpikeTrain2();
    //testHawkesMLEstimateUpdating();
    //testOptimization();
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

