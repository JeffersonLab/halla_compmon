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
  // Parse command line options
  std::string configFile("compmon.config");
  int maxEvents = -1;
  int run = -1;
  bool paramsFound = false;
  bool paramsError = false;
  for(int i = 1; i < argc && !paramsError; i++) {
    std::string paramName;
    std::string paramValue;
    std::string arg(argv[i]);
    size_t pos1 = arg.find_first_of("-");
    if(pos1 != std::string::npos && pos1 == 0) { // found a -
      size_t length = arg.length();
      if(length > 1 ) { // Valid parameter
        paramsFound = true;
        paramName=argv[i][1];
        size_t pos2 = arg.find_first_of("=");
        if(pos2 != std::string::npos && pos2==2) { // found an equals
          if(length <=3 ) {
            paramsError=true;
            std::cerr << "Invalid command line option: '"
              << argv[i] << "'" << std::endl;
          } else {
            paramValue = arg.substr(3,length);
          }
        } else if(length > 2) { // No equal sign used, still valid
          paramValue = arg.substr(2,length);
        } else { // Next parameter must be the actual value
          if(i+1 < argc) { // Ensure there are enough parameters
            i++;
            paramValue = argv[i];
          } else {
            std::cerr << "Missing argument to command: '" << argv[i] << "'"
              << std::endl;
            paramsError = true;
          }
        }
      } else {
        paramsError = true;
        std::cerr << "Invalid command line option: '"
          << argv[i] << "'" << std::endl;
      }
    } else { // Must be the default run number
      paramName = "r";
      paramValue = argv[i];
      paramsFound = true;
    }

    // Finally, get the required and optional parameters
    if(paramsError)
      continue;
    if(paramName.compare("r") == 0) {
      run = atoi(paramValue.c_str());
    } else if (paramName.compare("c") == 0) {
      configFile = paramValue;
    } else if (paramName.compare("n") == 0) {
      maxEvents = atoi(paramValue.c_str());
    } else {
      paramsError = true;
      std::cerr << "Invalid command line option: '"
        << argv[i] << "'" << std::endl;
    }
  }

  if(run < 0) {
    paramsError = true;
    std::cerr << "Invalid run number specified." << std::endl;
  }

  if(!paramsFound || paramsError) {
    std::cerr << "Usage: " << std::endl
      << "  " << argv[0] << " [options] runnumber " << std::endl
      << "    where options are: " << std::endl
      << "      -c config_file" << std::endl
      << "      -n maxEvents"  << std::endl;
    return -1;
  }

  // Open up the config file
  TString dataPath = "/data/cmu";
  TString dataPrefix = "FadcCalo_";
  TString dataPostfix = ".dat";
  TString outPath = "/data/cmuwork";
  TString outPrefix ="compmon_";
  TString paramsFileName = Form("%s/config/compmon.params", getenv("COMPMON_DIR"));
  // First, open up the config file and change any hard coded setting here
  fstream inConfig;
  inConfig.open(configFile,std::ios::in);
  std::string tmpConfigLine;
  while( inConfig >> tmpConfigLine && !inConfig.eof() ) {
    size_t pos1 = tmpConfigLine.find_first_of("=");
    if(pos1 == std::string::npos) {
      std::cerr << "Invalid line in config file: " << tmpConfigLine
        << std::endl;
      return -2;
    }
    std::string varName = tmpConfigLine.substr(0,pos1);
    std::string varValue = tmpConfigLine.substr(pos1+1,tmpConfigLine.length());

    if(varName.compare("dataPath")==0) {
      dataPath = varValue;
    } else if ( varName.compare("dataPrefix") == 0) {
      dataPrefix = varValue;
    } else if ( varName.compare("dataPostfix") == 0 ) {
      dataPostfix = varValue;
    } else if ( varName.compare("outPath") == 0 ) {
      outPath = varValue;
    } else if ( varName.compare("outPrefix") == 0 ) {
      outPrefix = varValue;
    } else if ( varName.compare("paramsFileName") == 0) {
      paramsFileName = varValue;
    }
  }
  TString dataFileName = TString::Format("%s/%s%d%s",dataPath.Data(),
      dataPrefix.Data(), run,dataPostfix.Data());
  TString outFileName = TString::Format("%s/%s%d.root",outPath.Data(),
      outPrefix.Data(),run);
  TString outHistosFileName = TString::Format("%s/%shistos_%d.root",outPath.Data(),
      outPrefix.Data(),run);

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
  int crlVersion = 1; /* The CRL version. Read from CODA file, otherwise
                         asume Spring 2016 running. */
  //
  // Setup intput and output files.  
  //
  /* Old method
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

  // TODO: Must fix this. For now, will hardcode a path because really long
  // (multi day) runs tend to break this simply lnkCompMon.output routine
  //char* outfilename =Form("lnkCompMon.output");
  char* outfilename =Form("/home/cornejo/scratch/compton/rootfiles/jc2_%d.root",run);
  cout<<"output file: "<<outfilename<<endl;
  */

  if(codaData->codaOpen(dataFileName, "r") != 0) {
    std::cerr << "\nERROR: Cannot open CODA file: "
      << dataFileName.Data() << std::endl;
    return -1;  
  }


  TFile* rootfilehistos = new TFile(outHistosFileName,"RECREATE","fadc histo data");
  TFile* rootfile = new TFile(outFileName,"RECREATE","fadc data");
  rootfile->SetCompressionLevel(0);
  rootfile->cd();
  //
  // setup parameters  (edit fadcparams.h to change values)
  // then setup histograms and trees
  //
  theTextParams->updateParamBuffer(paramsFileName.Data(),run);
  theTextParams->printParamList();

  theComptonStatus->DefineTrees(); //run-wise and epics-wise trees
  theComptonStatus->newRun();  //initialize status-tracking 

  theVMEauxdata->newRun();

  rootfilehistos->cd();
  theTriggered->DefineHistos();   
  theAccums->DefineHistos();
  rootfile->cd();
  theTriggered->DefineTriggeredTree();
  theTriggered->newRun();

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
	theComptonStatus->DebugDump(-1);	 //dump some status info
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
	  //status=theAccums->BuildQuad();
	  // Where a multiplet can be anything (pair, quad, octet, etc...)
	  status=theAccums->BuildMultiplet();
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
      /* Check to see if the version information got stored. If so, the
       * readout will depend on the version.
       * For Spring 2016 running, version information was not set.
       * For Fall 2016 running, version information is set and will
       * follow the new readout, accordingly.
       */
      j = 3+(len+3)/4; // Skip to the end of the flags input
      if( (buffp[j] >> 16 ) == 0xFADC ) {
        crlVersion = buffp[j] & 0xFFFF;
        theFADCdata->SetCRLVersion(crlVersion);
        printf("\n\nFound CRL verison in CODA file: CRL Version %d\n",crlVersion);
        j++;

        // If we found versioning info, then that means it's at least CRL
        // version 3 or newer. So we now look to see if we implemented
        // the new waveformReadout
        if( ( buffp[j] >> 16 ) == 0xFBE3 ) {
          printf("NewWaveformReadout is %s\n",(buffp[j]&0xFFFF ? "ENABLED" :
                "DISABLED" ));
          theFADCdata->SetWaveformReadout(buffp[j] & 0xFFFF);
        }
        if(crlVersion >= 5) {
          j++;
          if( ( buffp[j] >> 16 ) == 0xFBE4 ) {
            printf("MM is %s\n",(buffp[j]&0xFFFF ? "ENABLED" :
                "DISABLED" ));
            theFADCdata->SetMMEnabled(buffp[j] & 0xFFFF);
          }
        } else {
          theFADCdata->SetMMEnabled(2); // Set it to two so it won't break compatibility with old data
        }
        printf("\n");
      }
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
  // Apparently if it creates a sub file, we get a crash. So this is kind of a cludge to fix that
  rootfile = theAccums->GetTree()->GetCurrentFile();
  rootfile->Write();
  rootfile->Close();
  rootfilehistos->Write();
  rootfilehistos->Close();
  codaData->codaClose();
  return 0;
}
