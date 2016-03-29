/*
 * textParams.h
 *
 *  Created on: Dec 8, 2015
 *      Author: greggfranklin
 */

#ifndef TEXTPARAMS_H_
#define TEXTPARAMS_H_

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;
#include "TString.h"

#define BUFFER_LEN 5000
#define LABEL_LEN 50
#define VALUE_LEN 100

class textParams {
public:
	int findFlag(const char* flag); //see if flag exists
	int getString(char* value, const int maxsize, const char* flag);  //return string
	long getInt(const char* flag);
	float getFloat(const char* flag);
	int getRunNumber(){return runNumber;};
	void updateParamBuffer(const char* fileName,int runNumber);
	void printParamBuffer();
	void printParamList();
	textParams();
	virtual ~textParams();

private:
	int runNumber;
	char paramBuffer[BUFFER_LEN];
	void append(char* output, int outputLen, const char* input);
	char getToken(char output[], int* indexOut, int*indexIn,
			int outputSize, char inputString[]);
	void removeRepeats(char buffer[]);  //remove redundant flags
};

#endif /* TEXTPARAMS_H_ */
