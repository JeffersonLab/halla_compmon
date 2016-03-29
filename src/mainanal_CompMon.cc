// mainanal_CompMon.cc
//
//  analysis routine for fadccoda
//  Simplified version designed for monitoring Compton fadc data
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <bitset>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include "root.h"
#include "TString.h"
#include "THaEpics.h"
#include "THaCodaFile.h"
#include "THaCodaData.h"
#include "THaEpics.h"
#include "TTree.h"
#include "bankstructure.h"         //coda bank structure
//COMPTON fadc classes

#include "comptonStatus.h"         //tracks helicity, laser state, etc.
#include "textParams.h"
#include "fadcdata.h"         //unpacks fadc dat
#include "vmeauxdata.h"     //unpacks VME auxillary data
#include "fadcTriggered.h"    //analyzes FADC triggered data
#include "fadcAccums.h"      //analyzes FADC accumulator data

int main(int argc, char** argv)
{
  //+++++++++++++++++++++++++++++++
  //  initialize instances of classes, etc.
  //
  //Hall A general CODA and EPICS classes
  THaCodaFile *codaData;    //the Compton data
  THaEpics* epicsHandler=new THaEpics();  //handles epics data
  textParams* theTextParams=new textParams();
  fadcdata* theFADCdata=new fadcdata(theTextParams);  // unpacks FADC data
  //tracks helicty, laser state, etc.
  comptonStatus* theComptonStatus=new comptonStatus(theTextParams);
  //auxillary data from VME crate
  vmeauxdata* theVMEauxdata=new vmeauxdata(theTextParams);
//triggered and summed fadc data
  fadcTriggered* theTriggered=new fadcTriggered(theTextParams,theComptonStatus);
  //accumulator data
  fadcAccums* theAccums=new fadcAccums(theTextParams,theComptonStatus);

  codaData=new THaCodaFile();
  //
  // Setup intput and output files.  
  //
  TString fileName ;
  int run;
  // 
  // Open CODA data file
  if (argc>1)
  {
    run = atoi(argv[1]);
    //    fileName=Form("/home/franklin/HallA/newCompton/data/run%d.dat",run);
    fileName=Form("lnkCompMon.input");  //let script setup link to input file
    if (codaData->codaOpen(fileName, "r") != 0)
      {
	cout << "\nERROR: Cannot open CODA file\n";
	return -1;  
      } 	// end can't open CODA file
    cout<<"input file: "<<fileName<<endl;
  }else{
    cout<< "Error:  Called with no run number arguement"<<endl;
  }
  int maxEvents;
  if(argc>2){
    maxEvents=atoi(argv[2]);
  }else{
    maxEvents=-1;  //flag for no maxEvents count
  }
  if(maxEvents>0) printf(" WARNING maxEvents set to %d\n",maxEvents);
  //
  // Open ROOT outputfile  (keep it simple for now)
  //  char* outfilename =Form("fadc.root");
  char* outfilename =Form("lnkCompMon.output");
  cout<<"output file: "<<outfilename<<endl;

  TFile* rootfile = new TFile(outfilename,"RECREATE","fadc data");
  rootfile->SetCompressionLevel(0);
  //
  // setup parameters  (edit fadcparams.h to change values)
  // then setup histograms and trees
  //
  theTextParams->updateParamBuffer("compmon.params",run);
  theTextParams->printParamList();

  theComptonStatus->DefineTree(); //run-wise tree
  theComptonStatus->newRun();  //initialize status-tracking 

  theVMEauxdata->newRun();

  theTriggered->DefineHistos();   
  theTriggered->DefineTriggeredTree();
  theTriggered->newRun();

  theAccums->DefineHistos();
  theAccums->DefineTree();
  theAccums->newRun();
  // initialize for new run
  int status;
  bool laserTransition;
  int eventType;
  uint32_t* pROCdata;    //will get pointer to ROC data for each event
  int counter=0;     //counts all CODA events read
  int accumulatorCounter=0;  //counts physics events (accumulator readouts)
  int epicsCounter=0;  //counts EPICS events
  int usereventCounter=0;  //count special USER eents
  //  int fadcROCnum=6;   //search for ROC number 6 data 
  int fadcROCnum=0;   //search for ROC number (zero-> use first ROC found) 

  int verbose=0;
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // "first pass" analyzes the FADC data, etc.
  //
  //Loop through CODA data
  // codaData->codaClose();
  // codaData->codaOpen(fileName, "r");
  printf("BEGIN PASS ONE\n");
  int calorimeterChan=theTextParams->getInt("channel_calorimeter_PMT"); //fadc chan used for PMT
  while (eventType>=-2 && eventType!=20  &&
	 (counter<maxEvents||maxEvents<0)){
    counter++;
    status=codaData->codaRead();  //read a CODA event buffer in
    if(status<0) break;
    //codaData->evDump(1);
    eventType=codaData->UnpackEventHeader(verbose);
    //1 -> standard physics event
    //16 synch event, 17 PreStart, 18 Go, 19 Pause, 20 End
    // 131 epics event?
    // 121  USER event used for CONFIGUATION info
    if(eventType<0) {
      printf("event %d eventcode %x", counter, eventType);
    } else if(eventType==1) {
      pROCdata=codaData->getROCdata(fadcROCnum,0); //assume any ROC is THE ROC for now
      if(pROCdata!=NULL){
	theFADCdata->newMPS();     //clear counters, ,etc.
	status=theFADCdata->UnpackAllSubbanks(codaData,theVMEauxdata);
	//	theVMEauxdata->DebugDump(0);
	int codaEventNumber=codaData->GetEventNumber();
	//
	// finished upacking raw data,  now start analysis 
	// determine laserState, helicity, beamON, etc. with theComptonStatus
	//
	laserTransition=theComptonStatus->newMPS(codaEventNumber,
						 theFADCdata,theVMEauxdata);
	//ComptonStatus DebugDump	
	// mode =1  dump line every helicity window
	// mode =2  dump transitions only
	theComptonStatus->DebugDump(2);	 //dump some status info
	//
	// DUMP a few events for checking scalers and bits
	if(codaEventNumber>=100 && codaEventNumber<=105){
	  printf("Dump Scalers codaEvent=%d    ",codaEventNumber);
	  theVMEauxdata->DebugDump(0);  //dump TIR, scalers, etc.
	}
	if(laserTransition){
	  status=theAccums->DoLaserTransition(0); //do completed laser state
	  //(also call this at end of run)
	}
      	//if(theFADCdata->IsSumsValid(calorimeterChan)){   //trigged sums
	  status=theTriggered->DoSummedPulses(theVMEauxdata,theFADCdata);
	  //}
	if(theFADCdata->IsSamplesValid(calorimeterChan)){   //waveform samples
	  status=theTriggered->DoSampledWaveforms(theFADCdata);
	}
	if(theFADCdata->IsAccumValid(calorimeterChan)){  //FADC accumulator data
	  status=theAccums->DoAccums(theVMEauxdata,theFADCdata);
	  status=theAccums->BuildQuad();
	  accumulatorCounter++;
	}
	
      }
    } else if(eventType==131){  //EPICS
      epicsCounter++;
      //      theEPICSdata->Unpack(epicsHandler,codaData->getEvBuffer());
      theComptonStatus->UnpackEpics(epicsHandler,codaData->getEvBuffer());
      if(epicsCounter%1000==5) {
	printf("EPICS Event %d\n", epicsCounter);
	theComptonStatus->EpicsDebugDump();
      }
    } else if(eventType==121){
      printf("CONFIG event encountered\n");
      usereventCounter++;
      uint32_t *buffp=codaData->getEvBuffer();
      char* flagLineRaw=(char*)(buffp+3);
      char flagLine[2000];
      int len=(int) buffp[2];
      printf("Length of Flags String:%d\n",len);
      int j;
      printf("  ");
      for (int i=0; i<len+4; i++){
	j= i + 3 -2*(i%4);
	flagLine[i]=flagLineRaw[j];
	printf("%c",flagLine[i]);
	if(flagLine[i]==',') printf("\n  ");
	if(flagLine[i]=='\0') break;
      }
      printf("\n");
      flagLine[len]=NULL;
      } else if(eventType==17){
      } else if(eventType==18){
      } else if(eventType==20){
      }else{
	printf(" Unrecognized eventType=%d\n",eventType);
      }
    if(counter%1000==0){
      status=theTriggered->DoNormalizedHistos();
    }
    if(counter%10000==0){
      printf(" counter %d eventType %d \n",counter,eventType);
    }
    if(eventType==20) {
      break;
    }
  }	// end of first-pass analysis
  printf("Events Encounted:    %10d\n",counter);
  printf(" EPICS events;       %10d\n",epicsCounter);
  printf(" user(flags) events: %10d\n",usereventCounter);
  status=theAccums->DoLaserTransition(1);//wrap up last laser-period if needed
  rootfile->Write();
  codaData->codaClose();
  rootfile->Close();
  return 0;
}
