//
//  SubSpikeStruct.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#ifndef hawkesInCpp_SubSpikeStruct_h
#define hawkesInCpp_SubSpikeStruct_h

#include <list>
#include <vector>
#include "MyVector.h"
using namespace std;

class spike {
    
	float time;  // hh*3600+mm*60+ss+ff    
	string type; // Type (m)
    double amp; // magnitude (x)
    
	int UIN;  // SpikeID
    
public:
    string getType() const{
        return  type;
    }
    void setType(string type){
        this->type = type;
    }
    double getTime() const{
        return time;
    }
    void setTime(double time){
        this->time = time;
    }
    double getAmp() const{
        return amp;
    }
    void setAmp(double amp){
        this->amp = amp;
    }
    // other functions
    void printSpike(){
        std::cout << time << "," << type << endl;
    }
    
    const int getUIN(){
        return UIN;
    }
    void setUIN(int UIN){
        this->UIN = UIN;
    }
    
    
};

class SubSpikeStruct {
    
public:
    vector<spike> actualData;
	vector<spike*> entireList;
	vector<vector<spike*> > individualLists;
    
	int size;
	vector<int> individualSizes;
	vector<string> eventTypes;
	int numTypes;
    
	int startRefIndex;
	int endRefIndex;
    
public:
    SubSpikeStruct();
    SubSpikeStruct(int numTypes, string* eventTypesPtr);
    
    // Get and Set
    int getSize();
    void setSize(int size);
    vector<int> getIndividualSizes();
    MyVector<double> getTimingsVector_All();
    MyVector<double> getTimingsVector_ForType(int type);
    int getSizeOfType(int type);
    const int getNumTypes();
    void setNumTypes(int numTypes);
    void setEvents(vector<string> eventTypes);
    
    //
    SubSpikeStruct getSpikesBetween(int startEvent, int endEvent);
    SubSpikeStruct filterEventsOfType(vector<string> wantedTypes);
    string getLastEventType() const;int getLastEventTypeAsInteger();
    int getIndexOfType(string type);
    
    //printing
     void writeData2CSV(const char * fileName); // to file
     void printSpikeTrain(); // all outputs
     void printSummary();
};




#endif
