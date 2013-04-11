//
//  csvManagement.cpp
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/11/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#include <iostream>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string.h>
#include <list>
#include "stdio.h"
#include <stdlib.h>
#include "spikeTrain.h"
#include "csvManagement.h"

using namespace std;

double csvManagement::string_to_double( const std::string& s ){
	std::istringstream i(s);
	double x;
	if (!(i >> x))
		return 0;
	return x;
}

SubSpikeStruct csvManagement::csvToSpikeTrain(const char * fileName, int numTypes, vector<string> eventTypes ){
 	ifstream inFile(fileName);
    
    SubSpikeStruct ss; ss.setNumTypes(numTypes);ss.setEvents(eventTypes);
	ss.size = 0;
    
    string line;
	int linenum = 0;
    
	getline(inFile,line); // remove header.
    //cout << line << " - header " << endl;
    
	spike referenceSpike; // spike at t=- (initial time)
    
	while (getline (inFile, line))
	{
		istringstream linestream(line);
		linenum++;
        // cout << "LineNum:" << linenum << " text: " << line << endl;
		string item,hour,minute,second;
        double frac,hour_d,minute_d,seconds_d;
		spike s;
        getline (linestream, hour, ':' );  hour_d = string_to_double(hour);//s.setHour(string_to_double(hour));
		getline (linestream, minute, ':' );minute_d = string_to_double(minute);//s.setMinute(string_to_double(minute));
		getline (linestream, second, ',' );seconds_d = string_to_double(second); //s.setSeconds(string_to_double(second));
		getline (linestream, item, ',' ); frac = string_to_double(item);//s.setFractionOfSec(string_to_double(item));
		double time = hour_d*3600+minute_d*60+seconds_d+frac;
        s.setTime(time);
        
		getline (linestream, item, ',' );
        item = item.substr(0,1);
        s.setType(item);
		
		if(linenum==1){
			referenceSpike = s;
		}
        
		bool typeOfInterest = false;
        
		for(int i=0; i < numTypes; i++){
            if( item == eventTypes[i]){
				typeOfInterest = true;
            }
            
		}
        
		if(typeOfInterest){
			ss.actualData.push_back(s);
            ss.entireList.push_back(&s);
            ss.size++;
            string eventType = s.getType();
            int typeOfIthSpike = ss.getIndexOfType(eventType);
            if(typeOfIthSpike == -1){
                cout << "New Type Detected. Event type is " << eventType << endl;
            }
            else{
                ss.individualLists[typeOfIthSpike].push_back(&s);
                ss.individualSizes[typeOfIthSpike]++;
            }
        }
	}

    ss.startRefIndex = 0;
    ss.endRefIndex = ss.size-1;
    
    return ss;
}

SubSpikeStruct csvManagement::csvToSpikeTrain(const char * fileName){
    SubSpikeStruct ss;
    ifstream inFile(fileName);
    
	ss.size = 0;
    ss.numTypes = 0;
    
    string line;
	int linenum = 0;
    
	getline(inFile,line); // remove header.
    //cout << line << " - header " << endl;
    
	while (getline (inFile, line))
	{
		istringstream linestream(line);
		linenum++;
        // cout << "LineNum:" << linenum << " text: " << line << endl;
		string item,hour,minute,second;
        double frac,hour_d,minute_d,seconds_d;
		spike s;
        getline (linestream, hour, ':' );  hour_d = string_to_double(hour);//s.setHour(string_to_double(hour));
		getline (linestream, minute, ':' );minute_d = string_to_double(minute);//s.setMinute(string_to_double(minute));
		getline (linestream, second, ',' );seconds_d = string_to_double(second); //s.setSeconds(string_to_double(second));
		getline (linestream, item, ',' ); frac = string_to_double(item);//s.setFractionOfSec(string_to_double(item));
		double time = hour_d*3600+minute_d*60+seconds_d+frac;
        s.setTime(time);
        
		getline (linestream, item, ',' );
        item = item.substr(0,1);
        s.setType(item);
		
        // if ss.eventTypes contains item, then get its Index and push the spike up that type.
        
        bool eventTypesContainsElement = false; int typeIndex = -1;
        
        for(int tempit = 0; tempit < ss.eventTypes.size(); tempit++){
            if (item == ss.eventTypes[tempit]) {
                eventTypesContainsElement = true; typeIndex = tempit;
            }
        }
        
        if(!eventTypesContainsElement){
            ss.numTypes = ss.numTypes+1; ss.eventTypes.push_back(item); typeIndex = ss.numTypes;
            ss.individualSizes.push_back(0); 
            vector<spike* > vs; ss.individualLists.push_back(vs);
        }
        ss.actualData.push_back(s);
        ss.entireList.push_back(&s);
        ss.size++;
        ss.individualLists[typeIndex].push_back(&s);
        ss.individualSizes[typeIndex]++;
        
	}
    
    ss.startRefIndex = 0;
    ss.endRefIndex = ss.size-1;
    
    return ss;
}
