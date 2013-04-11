//
//  spikeTrain.cpp
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#include "spikeTrain.h"
#include "stdio.h"
#include <iostream>

/*
 *  Big Three
 */
/*************************************************************************************************/
// copy constructor, destructor and operator=
// clear and copy.


void spikeTrain::copy(spikeTrain const & copyObj ){
	size = copyObj.size;
	numTypes = copyObj.numTypes;
    
	if(copyObj.spikeTrainArrayPtr != NULL){
		this->spikeTrainArrayPtr = new spike[size];
		for(int i=0; i < size ; i++){
			this->spikeTrainArrayPtr[i] = copyObj.spikeTrainArrayPtr[i];
		}
        //		cout << "spike Train Pointer copy is completed." << endl;
	}
	else{
		spikeTrainArrayPtr = NULL;
        //		cout << "spike Train Pointer is null." << endl;
	}
    
	if(copyObj.eventTypesPtr!=NULL){
		this->eventTypesPtr = new string[numTypes];
		for(int i=0; i < numTypes ; i++){
			this->eventTypesPtr[i] = copyObj.eventTypesPtr[i];
		}
        //		cout << "event Types Pointer copy is completed." << endl;
	}
	else{
		eventTypesPtr = NULL;
        //		cout << "event Types Pointer is null." << endl;
	}
}

void spikeTrain::clear(){
    //	cout << "Entered." << endl;
	if(spikeTrainArrayPtr!= NULL){
        //		cout << "Clearing spike Train Pointer Array." << endl;
		delete [] spikeTrainArrayPtr;
        //		cout << "spike Train Pointer Array cleared." << endl;
	}
	else{
        //		cout << "SPike Train pointer is null." << endl;
	}
    
	if(eventTypesPtr != NULL){
        //		cout << "Clearing et Pointer Array." << endl;
		delete [] eventTypesPtr;
        //		cout << "event Types pointer cleared." << endl;
	}
	else{
        //		cout << "etp is null." << endl;
	}
    //	cout << "Leaving" << endl;
}

spikeTrain::spikeTrain() {
    
}

spikeTrain::spikeTrain(spikeTrain const & copyObj){
	copy(copyObj);
}

spikeTrain const & spikeTrain::operator=(spikeTrain const & rhs)
{
	if(this != &rhs)
	{
		clear();
		copy(rhs);
	}
	return *this;
}

spikeTrain::~spikeTrain() {
	clear();
}

spikeTrain::spikeTrain(int size, int numTypes) {
	this->spikeTrainArrayPtr = new spike[size];
	this->numTypes = numTypes;
	this->eventTypesPtr = new string[numTypes];
	mainStruct.create(numTypes);
}

spikeTrain::spikeTrain(int size, int numTypes, string* eventTypesPtr) {
	this->spikeTrainArrayPtr = new spike[size];
	this->numTypes = numTypes;
	this->eventTypesPtr = new string[numTypes];
	mainStruct.create(numTypes);
	for(int i=0; i < numTypes; i++){
        this->eventTypesPtr[i] = eventTypesPtr[i];
        this->mainStruct.eventTypesPtr[i] = eventTypesPtr[i];
	}
}




/**********************************************************************************/

void spikeTrain::createMainStruct(){
	mainStruct.size = 0;
	for(int i=0; i < numTypes ; i++){
		mainStruct.individualSizes[i] = 0;
	}
    
	for(int i=0; i < size ; i++){
		mainStruct.entireList.push_back(&spikeTrainArrayPtr[i]);
		mainStruct.size++;
		string eventType = spikeTrainArrayPtr[i].getType();
		int typeOfIthSpike = getIndexOfType(eventType);
		if(typeOfIthSpike == -1){
			cout << "New Type Detected. Event type is " << eventType << endl;
		}
		else{
			mainStruct.individualLists[typeOfIthSpike].push_back(&spikeTrainArrayPtr[i]);
			mainStruct.individualSizes[typeOfIthSpike]++;
		}
	}
    
	mainStruct.startRefIndex = 0;
	mainStruct.endRefIndex = size-1;
    
}

vector<spike*> spikeTrain::getListOfEventsOfType(string type){
	return getListOfEventsOfType(getIndexOfType(type));
}

vector<spike*> spikeTrain::getListOfEventsOfType(int type){
	return mainStruct.individualLists[type];
}

void spikeTrain::setSpikeAt(int i,spike s){
	this->spikeTrainArrayPtr[i] = s;
	this->spikeTrainArrayPtr[i].setUIN(i);
}

void spikeTrain::setSize(int size){
	this->size = size;
}

int spikeTrain::getSize(){
	return size;
}

void spikeTrain::setEventTypesPtr(string* ptr){
	for(int i=0; i < numTypes; i++){
		this->eventTypesPtr[i] = ptr[i];
		this->mainStruct.eventTypesPtr[i] = ptr[i];
	}
}

const int spikeTrain::getNumTypes(){
	return numTypes;
}

SubSpikeStruct spikeTrain::getMainStruct(){
	return mainStruct;
}

void spikeTrain::printSpikeTrain(){
	cout << "Priting starts, size of spike train is " << size << endl;
	for(int i=0; i< size; i++){
		cout << "t = " << spikeTrainArrayPtr[i].getTime() << ", type = " << spikeTrainArrayPtr[i].getType() << endl;
	}
    
	cout << "printing done." << endl;
}

void spikeTrain::printTypes(){
	for(int i=0; i<numTypes; i++){
		cout << eventTypesPtr[i] << ",";
	}
	cout << endl;
}

int spikeTrain::getIndexOfType(const string type){
	// input - (string) type of event.
	// output - index of type.
    
	for(int i=0; i<numTypes; i++){
		if(eventTypesPtr[i] == type){
			return i;
		}
	}
    
	cout << "mentioned event type not in the event list." << endl;
	return -1;
    
}


void spikeTrain::printListOfEventsOfType(const string type){
	vector<spike*> listOfSpikePtrs = getListOfEventsOfType(type);
	cout << "Printing Events of Type " << type << endl;
	printListOfSpikePtrs(listOfSpikePtrs);
}

void spikeTrain::printListOfSpikePtrs(vector<spike*> listOfSpikePtrs){
	int size = listOfSpikePtrs.size();
	for(int i=0; i< size; i++){
		listOfSpikePtrs[i]->printSpike();
	}
}


