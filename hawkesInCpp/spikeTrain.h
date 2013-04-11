//
//  spikeTrain.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#ifndef hawkesInCpp_spikeTrain_h
#define hawkesInCpp_spikeTrain_h

#include "spike.h"
#include "SubSpikeStruct.h"
#include <list>
#include <vector>
using namespace std;

class spikeTrain {
private:
	spike* spikeTrainArrayPtr; // actual data
    
	// params
	int size;
	int numTypes;
	string* eventTypesPtr;
	SubSpikeStruct mainStruct;
    
public:
	spikeTrain();
	spikeTrain(int size, int numTypes);
	spikeTrain(int size, int numTypes, string* eventTypesPtr);
	spikeTrain(spikeTrain const & copyObj);
	spikeTrain const & operator=(spikeTrain const & rhs);
	virtual ~spikeTrain();
    
    
	// things to be set while reading data.
	void setSpikeAt(int i,spike s);
    
	int getIndexOfType(const string type);
	void createMainStruct();
    
	SubSpikeStruct getLastKSpikes(int K);
    
    
	// getters and setters
	int getSize();
	const int getNumTypes();
	void setSize(int size);
	void setNumTypes(int numTypes);
	void setEventTypesPtr(string* ptr);
	vector<spike*> getListOfEventsOfType(int type);
	vector<spike*> getListOfEventsOfType(string type);
	SubSpikeStruct getMainStruct();
    
	// printing
	void printSpikeTrain();
	void printTypes();
	void printListOfEventsOfType(const string type);
	void printListOfSpikePtrs(vector<spike*> listOfSpikePtrs);
    
private:
	void copy(spikeTrain const & copyObj);
	void clear();
    
    
	//vector<spike*>* getEntireList();
    
	// other functions *************************************************************************
	// getting events of particular type m'
	//vector<spike*>* getListsOfAllEventTypes();
	//vector<spike*> getSpikesBetween(double time0, double time1, const int typeIndex);
    
    
    
};

#endif
