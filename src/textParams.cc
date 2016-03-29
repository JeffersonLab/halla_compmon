/*
 * textParams.cpp
 *
 *  Created on: Dec 8, 2015
 *      Author: greggfranklin
 */

#include "textParams.h"

int textParams::getString(char *value, const int maxSize, const char *flag)
{
	int rtn=findFlag(flag);
	if(rtn>0){
		int indexBuffer=rtn; //first char. of value
		int indexOut=0;
		while(indexBuffer<BUFFER_LEN&&indexOut<maxSize){
			value[indexOut]=paramBuffer[indexBuffer++];
			if(value[indexOut]==',')break;
			indexOut++;
		}
		value[indexOut]='\0';
	} else {
	  printf("WARNING: Missing parameter flag: %s\n",flag);
	}
	return rtn;
}

int textParams::findFlag(const char *flag)
{
	// looks up flag in parameter input string buffer
	// rtn= 0  if flag not found
	// rtn= -1 if flag found, but not followed by "=:
	// rtn>0  if flag has a value (rtn=index of first charcter of value)
	int rtn;
	char* pos=strstr(paramBuffer,flag);
	if(pos==NULL){
		rtn=0;
	}else {
		rtn=-1;     //We've found it, does it also have a value?
		int indexFlag=pos-paramBuffer;
		for(int i=indexFlag; i<BUFFER_LEN; i++){
			if(paramBuffer[i]=='='){
				rtn=i+1;  //start of Value token
				break;
			}
			if(paramBuffer[i]==',') break;
		}
	}
	return rtn;
}

char textParams:: getToken(char output[], int* indexOut, int* indexIn,
		int outputSize, char inputString[]){
	int len1=strlen(inputString);
	char rtnchr='\0';
	//find deliminating character
	for(int i=*indexIn; i<len1; i++){
		if(strchr("=;\0",inputString[i])!=NULL){
			rtnchr=inputString[i];
			if(rtnchr==';')rtnchr='\0';
			len1=i;
			break;
		}
	}
	int lastNonBlank=*indexIn;
	for(int i=*indexIn; i<len1; i++){
		if(inputString[i]!=' ') lastNonBlank=i;
	}
	int firstNonBlank=*indexIn;
	for(int i=*indexIn; i<len1; i++){
		if(inputString[i]!=' ')break;
		firstNonBlank++;
	}
	int j=*indexOut;
	for(int i=firstNonBlank; i<lastNonBlank+1; i++){
		output[j++]=inputString[i];
		if(j==outputSize-1) break;
	}
	output[j]='\0';
	*indexIn=len1+1;
	*indexOut=j;
	return rtnchr;
}
void textParams:: append(char output[], int maxSize, const char* inputString){
	int pntOut=strlen(output);
	int len1=strlen(inputString);
	for(int i=0; i<len1; i++){
		if(pntOut>maxSize-2) break;
		output[pntOut++]=inputString[i];
	}
	output[pntOut]='\0';
	return;
}
void textParams::removeRepeats(char buffer[]){
	char label[50];
	int indexIn=0;
	int indexOut=0;
	char* pnt;

	while (buffer[indexIn]!='\0'){
		label[0]=',';
		int indexLabel=1;
		int indexTmp=indexIn;
		for(int i=1;i<49;i++){
			if(buffer[indexTmp]=='\0' ) break;
			label[indexLabel++]=buffer[indexTmp++];
			if(label[indexLabel-1]=='=') break;
		}
		label[indexLabel]='\0';
		pnt=strstr(&buffer[indexIn],label);
		if(pnt!=NULL) {
			while (buffer[indexIn]!='\0'){
				if(buffer[indexIn++]==',') break;
			}
		}else{
			while (buffer[indexIn]!='\0'){
				buffer[indexOut++]=buffer[indexIn++];
				if(buffer[indexIn-1]==',') break;
			}
		}
	}
	buffer[indexOut]='\0';
	return;
}
long textParams::getInt(const char* flag){
	char valueString[VALUE_LEN];
	int foundIt=getString(valueString, VALUE_LEN, flag);
	long rtn=0;
	TString tmp=TString(valueString);
	if(foundIt>0){
		try{
		  rtn=tmp.Atoi();
		} catch (std::exception const &exc) {
			printf("Parameter Error: cannot convert value to integer: %s=%s\n",flag,valueString);
			rtn=0;
		}
	}
	return rtn;
}

float textParams::getFloat(const char* flag){
	char valueString[VALUE_LEN];
	int foundIt=getString(valueString, VALUE_LEN, flag);
	float rtn=0;
	TString tmp=TString(valueString);
	if(foundIt>0){
		try{
			rtn =tmp.Atof();
		} catch (std::exception const &exc) {
			printf("Parameter Error: cannot convert value to flaot: %s=%s\n",flag,valueString);
			rtn=0;
		}
	}
	return rtn;
}
void textParams::printParamList(){
	// print out flag=value format line by line
	for(int i=1;i<5000;i++){
		if(paramBuffer[i]=='\0') break;
		if(paramBuffer[i]==','){
			cout<<endl;
		}else{
			cout<<paramBuffer[i];
		};
	}
	cout<<endl;
	return;
}
void textParams::printParamBuffer(){
	// dump out list of parameters as single line
	printf("%s\n",paramBuffer);
	cout<<endl;
	return;
}
void textParams::updateParamBuffer(const char* fileName,int runNumberIn){
	printf("Parameter Update from file: %s run %d\n",fileName,runNumberIn);
	const char COMMENT_CHAR = ';';
	const int BUFF_IN_LEN=200;
	runNumber=runNumberIn;
	ifstream *fd=new ifstream();
	fd->open(fileName);
	if( !fd->is_open() ){
		printf("Failed to open parameter file\n");
	} else {
		char inputLine[BUFF_IN_LEN];
		char inputLabel[LABEL_LEN];
		char inputValue[VALUE_LEN];
		char* pntNext;
		int runTransition=0;
		paramBuffer[0]=','; //start with a deliminator for simplicity
		paramBuffer[1]='\0';
		while( !fd->eof() ){
			fd->getline(inputLine,BUFF_IN_LEN-1);
			
			pntNext=strchr(inputLine,COMMENT_CHAR);
			if(pntNext!=NULL)*pntNext='\0';
			int indexIn=0;
			int indexOut=0;
			int EOL=strlen(inputLine);
			char delim1, delim2;
			while(indexIn<EOL) {
				indexOut=0;
				indexIn=0;
				delim1=getToken(inputLabel,&indexOut,&indexIn,LABEL_LEN,inputLine);
				indexOut=0;
				if(delim1=='='){
					delim2=getToken(inputValue,&indexOut,&indexIn,VALUE_LEN,inputLine);
				} else {
					inputValue[0]='\0';
				}
				//Now decide if we put them into output buffer
				if(strcmp(inputLabel,"TRANSITION")==0){
					try{
					  runTransition=TString(inputValue).Atoi();
					} catch (std::exception const &exc) {
						printf("Cannot interpret TRANSITION=%s\n",inputValue);
					}
				} else if(runNumber>=runTransition) {
					append(paramBuffer,BUFFER_LEN,inputLabel);
					if(indexOut>0){
						append(paramBuffer,BUFFER_LEN,"=");
						append(paramBuffer,BUFFER_LEN,inputValue);
					}
					append(paramBuffer,BUFFER_LEN,",");
				}
			}
		}
		removeRepeats(paramBuffer);
	}
	return;
}
textParams::textParams() {
	// TODO Auto-generated constructor stub

}

textParams::~textParams() {
	// TODO Auto-generated destructor stub
}
