//
//  SubSpikeStruct.cpp
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include "SubSpikeStruct.h"
#include "stdio.h"
#include "iostream.h"


SubSpikeStruct::SubSpikeStruct() {
}

SubSpikeStruct::SubSpikeStruct(int numTypes, string* eventTypesPtr){
	this->numTypes = numTypes;
	for(int i=0; i < numTypes; i++){
		eventTypes.push_back(eventTypesPtr[i]);
	}
}

// Get and Set
int SubSpikeStruct::getSize(){
	return size;
}
void SubSpikeStruct::setSize(int size){
	this->size = size;
}
int SubSpikeStruct::getSizeOfType(int type){
	return this->individualSizes[type];
}
const int SubSpikeStruct::getNumTypes(){
	return numTypes;
}
void SubSpikeStruct::setNumTypes(int numTypes){
	this->numTypes = numTypes;
}

SubSpikeStruct SubSpikeStruct::getSpikesBetween(int startEvent, int endEvent){
    SubSpikeStruct ss; ss.setNumTypes(numTypes);ss.setEvents(eventTypes);
	if(startEvent < 0){
        cout << "starting number is negative\n";
        startEvent = 0;
    }
    if(endEvent >= size){
        cout << "ending number is greater than size of data\n";
        endEvent = size-1;
    }
    
    ss.size = 0; ss.individualSizes.resize(numTypes); ss.individualLists.resize(numTypes);
	for(int i=0; i < numTypes ; i++){
		ss.individualSizes[i] = 0;
	}
	for(int i=startEvent; i < endEvent+1 ; i++){
		ss.entireList.push_back(entireList[i]);
		ss.size++;
		string eventType = entireList[i]->getType();
		int typeOfIthSpike = getIndexOfType(eventType);
		if(typeOfIthSpike == -1){
			cout << "New Type Detected. Event type is " << eventType << endl;
		}
		else{
			ss.individualLists[typeOfIthSpike].push_back(entireList[i]);
			ss.individualSizes[typeOfIthSpike]++;
		}
	}
        
	return ss;
}

SubSpikeStruct SubSpikeStruct::filterEventsOfType(vector<string> wantedTypes){
    int wantedTypesSize = wantedTypes.size();
    SubSpikeStruct ss; ss.setNumTypes(wantedTypesSize);ss.setEvents(wantedTypes);
	
    ss.size = 0; ss.individualSizes.resize(wantedTypesSize); ss.individualLists.resize(wantedTypesSize);
	for(int i=0; i < wantedTypesSize ; i++){
        string wantedType = wantedTypes[i]; int wantedTypeIndex = this->getIndexOfType(wantedType);
		ss.individualSizes[i] = this->individualSizes[wantedTypeIndex];
        ss.individualLists[i] = this->individualLists[wantedTypeIndex];
	}
    
	for(int i=0; i < this->size ; i++){
		spike* sptr = entireList[i]; string typeOfsptr = sptr->getType();
        bool isWanted = false;
        for(int j = 0 ; j < wantedTypesSize; j++){
            if(typeOfsptr == wantedTypes[j]){
                isWanted = true; 
            }
        }
        if(isWanted){
            ss.entireList.push_back(sptr);
            ss.size++;
        }
	}
        
	return ss;
}


int SubSpikeStruct::getLastEventTypeAsInteger(){
    return getIndexOfType(this->entireList[size-1]->getType());
}
string SubSpikeStruct::getLastEventType() const{
    return this->entireList[size-1]->getType();
}

void SubSpikeStruct::setEvents(vector<string> eventTypes){
    this->eventTypes = eventTypes;
}

int SubSpikeStruct::getIndexOfType(string type){
	// input - (string) type of event.
	// output - index of type.

	for(int i=0; i<numTypes; i++){
		if(eventTypes[i] == type){
			return i;
		}
	}
    cout << numTypes << endl;
    for(int i=0; i<numTypes; i++){
        cout << eventTypes[i] << ",";
    }cout << "\n";
    cout << type << endl;
	cout << "mentioned event type not in the event list." << endl;
	return -1;
    
}

void SubSpikeStruct::writeData2CSV(const char *fileName){
    ofstream myfile;
    myfile.open(fileName);
    
    for(int i=0; i < size; i++){
        spike* s = this->entireList[i];
        double time= s->getTime();
        string type = s->getType();
        double magnitude = s->getAmp();
        
        myfile << time << "," << type << "," << magnitude << "," << endl;
    }
    myfile.close();    
}

MyVector<double> SubSpikeStruct::getTimingsVector_All(){
    int size = this->getSize();
    MyVector<double> ret(size);
    for(int i = 0; i < size; i++){
        ret.Set_Element(i, this->entireList[i]->getTime());
    }
    return ret;
}

MyVector<double> SubSpikeStruct::getTimingsVector_ForType(int type){
    if(type >= numTypes)
        cout << "invalid numType" << endl;
    vector<int> sizes = individualSizes;
    int size = sizes[type];
    MyVector<double> ret(size);
    for(int i = 0; i < size; i++){
        ret.Set_Element(i, this->individualLists[type][i]->getTime());
    }
    return ret;
}

// printings
void SubSpikeStruct::printSpikeTrain(){
    for(int i=0; i< size; i++){
		cout << "t = " << this->entireList[i]->getTime() << ", type = " << this->entireList[i]->getType() << endl;
	}
	cout << endl;
}

void SubSpikeStruct::printSummary(){
    cout << "number of types = " << numTypes << endl;
    for(int i=0;  i< numTypes; i++){
        cout << eventTypes[i] << ",";
    }
    cout << endl;
    cout << "size of data = " << size << ". Individual sizes are"<<endl;
    for(int i=0;  i< numTypes; i++){
        cout << individualSizes[i] << ",";
    }
    cout << endl;
}