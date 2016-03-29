/////////////////////////////////////////////////////////////////////
//
//  THaCodaFile
//  File of CODA data
//
//  CODA data file, and facilities to open, close, read,
//  write, filter CODA data to disk files, etc.
//  A lot of this relies on the "evio.cpp" code from
//  DAQ group which is a C++ rendition of evio.c that
//  we have used for years, but here are some useful
//  added features.
//
//  author  Robert Michaels (rom@jlab.org)
//
//  updated 2/1/2015 Gregg Franklin
/////////////////////////////////////////////////////////////////////

#include "THaCodaFile.h"

#ifndef STANDALONE
ClassImp(THaCodaFile)
#endif

//Constructors 

  THaCodaFile::THaCodaFile() {       // do nothing (must open file separately)
       ffirst = 0;
       init(" no name ");
  }
  THaCodaFile::THaCodaFile(TString fname) {
       ffirst = 0;
       init(fname);
       int status = codaOpen(fname.Data(),"r");       // read only 
       staterr("open",status);
  }
  THaCodaFile::THaCodaFile(TString fname, TString readwrite) {
       ffirst = 0;
       init(fname);
       int status = codaOpen(fname.Data(),readwrite.Data());  // pass read or write flag
       staterr("open",status);
  }

//Destructor
  THaCodaFile::~THaCodaFile () {
       int status = codaClose();
       staterr("close",status);
  };       

  int THaCodaFile::codaOpen(TString fname) {  
    cout<<"DEBUG: codaOpen called"<<endl;
    char inputfile[80];
    char rw[10];
    rw[0]='r';
    strcpy(inputfile,fname);
       int status = evOpen(inputfile,rw,&handle);
       staterr("open",status);
       return status;
  };

  int THaCodaFile::codaOpen(TString fname, TString readwrite) {  
    char inputfile[80];
    char rw[80];
    strcpy(rw,readwrite);
    strcpy(inputfile,fname);
      int status = evOpen(inputfile,rw,&handle);
      staterr("open",status);
      return status;
  };


  int THaCodaFile::codaClose() {
// Close the file. Do nothing if file not opened.
    if( handle ) {
      int status = evClose(handle);
      handle = 0;
      return status;
    }
    return CODA_OK;
  }

  int THaCodaFile::codaRead() {
// codaRead: Reads data from file, stored in evbuffer.
// Must be called once per event.
    int status;
    if ( handle ) {
       status = evRead(handle, evbuffer, MAXEVLEN);
       staterr("read",status);
       if (status != S_SUCCESS) {
  	  if (status == EOF) return status;  // ok, end of file
          status = CODA_ERROR;
       }
    } else {
      if(CODA_VERBOSE) {
         cout << "codaRead ERROR: tried to access a file with handle = 0" << endl;
         cout << "You need to call codaOpen(filename)" << endl;
         cout << "or use the constructor with (filename) arg" << endl;
      }
      status = CODA_ERROR;
    }
    return status;
  };



  int THaCodaFile::codaWrite(uint32_t *evbuffer) {
// codaWrite: Writes data to file
     int status;
     if ( handle ) {
       status = evWrite(handle, evbuffer);
       staterr("write",status);
     } else {
       cout << "codaWrite ERROR: tried to access file with handle = 0" << endl;
       status = CODA_ERROR;
     }
     return status;
   };



  uint32_t* THaCodaFile::getEvBuffer() {
// Here's how to get raw event buffer, evbuffer, after codaRead call
      return evbuffer;
  }

void THaCodaFile::evDump(int verbose=0){
  // HEX DUMP
  int NumWrds=(int) evbuffer[0];
  printf("evDump: Number of word in event %d\n",NumWrds);
  for (int i=0; i<NumWrds; i++){
    if(i%10==0) printf("\n %6d",i);
    printf("%11X ",evbuffer[i]);
  }
  printf("\n*******\n");
  return;
}
void THaCodaFile::evDumpStructure(int verbose=0){
  // DUMP of header information, including ROC and subbanks for event type 1
  int NumWrds=(int) evbuffer[0];
  EventType=( (evbuffer[1]>>16)&0xffff);
  EventDataType=( (evbuffer[1]>>8)&0xff);
  EventDataNum=(evbuffer[1]&0xff);
  printf(" evDumpStructure: NumWrds %5d ", NumWrds);
  printf("     Event Type=%d EventDataNum=%d\n",EventType,EventDataNum);
  if(EventType==1){
    //read out subbanks
    //first subbank is standard CODA 4 word subbank
    int Index=2;
    // for(int i=0; i<60; i++){
    //   printf(" 0x%8x ", evbuffer[2+i] );
    //   if(i%10==0) printf("\n");
    // }
    // printf("\n");
    int bankWords=evbuffer[Index];  //# of datawords 
    int Header=evbuffer[Index+1]; 
    int EventNumber=evbuffer[Index+2];
    //    int Data=&ROCData[Index+2];  //pointer to data in subbank
    int Tag=((Header>>16)&0xffff);
    printf("EventID Bank: Tag 0x%x NumWrds %d EventNumber %d\n",
	   Tag, bankWords, EventNumber);
    Index+=5;  //skip over header subbank
    // start loop on ROC data banks
    while(Index<NumWrds){
      ROCWords=(int) evbuffer[Index];
      ROCnum= (int) (evbuffer[Index+1]>>16)&0x1f;
      printf("ROC Bank: ROCnum= %d, Num Wrds= %d\n",ROCnum,ROCWords);
      int subIndex=Index+2;
      int SubbankWords;
      int SubbankType;
      int SubbankFlag;
      while (subIndex<ROCWords){
	SubbankWords=evbuffer[subIndex]-1;  //# of datawords 
	Header=evbuffer[subIndex+1]-1; 
	EventNumber=evbuffer[subIndex+2];
	//    int Data=&ROCData[Index+2];  //pointer to data in subbank
	SubbankType=((Header>>16)&0xffff);
	SubbankFlag=(Header&0xffff);
	printf("  Subbank: Tag 0x%4x NumWrds %5d Flag 0x%x\n",
	       SubbankType,SubbankWords,SubbankFlag);
	if(SubbankWords<1) break;
	subIndex+=SubbankWords+2;
      } 
      if(bankWords<1) break;
      Index+=ROCWords;
    }
    printf("\n*******\n");
  }
  return;
}

  int THaCodaFile::UnpackEventHeader(int verbose){
  // Can be called after codaRead
  // unpacks events record and returns EventType
    // verbose=0 for minimal output (event start/etc.)
    // verbose=-1 no output
    // verbse>0 to dump info for all events
    EventType=( (evbuffer[1]>>16)&0xffff);
    EventDataType=( (evbuffer[1]>>8)&0xff);
    EventDataNum=(evbuffer[1]&0xff); //don't know what this is (should be 0?
    if(verbose>0) printf("Event Type %d\n",EventType);
    if(EventType==SYNC_EVTYPE){
     if(verbose>=0) printf("SYNC event encountered\n");
    }else if(EventType==PRESTART_EVTYPE){
      StartTime=evbuffer[2];
      RunNumber=evbuffer[3];
      if(verbose>=0) {
	printf("PRESTART event encountered\n");
	printf("Run Number %d\n",RunNumber);
      }
      RunOnTime=0;         //sum total time (seconds) run is enabled
      GoTime=0;            //time of latest GO event
    }else if(EventType==GO_EVTYPE){
      printf("GO event encountered\n");
      //      GoTime=evbuffer[1];
      GoTime=evbuffer[2];  //fix?  619/15 gbf
    }else if(EventType==PAUSE_EVTYPE){
      printf("PAUSE event encountered\n");
      RunOnTime+=evbuffer[2]-GoTime;
      GoTime=0;
    }else if(EventType==END_EVTYPE){
      RunOnTime+=evbuffer[2]-GoTime;
      GoTime=0;
      EndTime=evbuffer[2];
      printf("END event encountered\n");
      printf("Total Run Time %d seconds. Run Enabled Time %d seconds\n",
	     EndTime-StartTime,RunOnTime);
    } else if(EventType>0 && EventType<=MAX_PHYS_EVTYPE){
      //read out event ID bank for physics events
      EventNumber=evbuffer[4];
      EventClass=evbuffer[5];
      EventStatus=evbuffer[6];
      if(verbose>1) printf(" Event Number %d Class %d Status %d\n",
			 EventNumber,EventClass,EventStatus);
    }else{
      EventNumber=-1;
      EventClass=-1;
      EventStatus=-1;
    }
   if(verbose>0){
     printf("Event Number      %d\n",EventNumber);
     if(verbose>1){
       for(int i=0; i<verbose; i++){
	 printf(" Buffer[%d]=%d\n",i,evbuffer[i]);
       }
     }
   }
   return EventType;
  }

uint32_t* THaCodaFile::getROCdata(int SeekROCnum=0, int verbose=0){
  // looks for first set of data from ROC number SeekROCnum
  //  finds with first ROC data encountered if SeekROCnum<=0)
  // sets pointers to  first bank
  //  LOOKS LIKE IT ONLY WORKS FOR ONE ROC ??
  int bankIndex=7;         //pointer to word count in first ROC bank
  ROCWords=(int) evbuffer[bankIndex];
  ROCnum= (int) (evbuffer[bankIndex+1]>>16)&0x1f;
  if(verbose>1) printf("found ROC# %d \n",ROCnum);
  if(SeekROCnum==ROCnum || SeekROCnum==0){
    ROCData= &evbuffer[bankIndex];
    CurrentSubbank.BankIndex=-2;       //flag for reseting subbank
    if(verbose>0) printf(" ROC %3d   Num ROC Data Wrds= %d\n",ROCnum,ROCWords);
    return (uint32_t*) &evbuffer[bankIndex];   //return pointer to ROC data
  }else{
    if(verbose>0) printf(" ROC %d not found in event\n",SeekROCnum);
    return NULL;
  }
}

bankstructure THaCodaFile::getNextSubbank(int verbose){
  if(verbose>1) printf("getNextSubbank: CurrentSubbank.BankIndes=%d\n",
		       CurrentSubbank.BankIndex);
  if(CurrentSubbank.BankIndex==-1){
    return CurrentSubbank;                //already have all the subbanks
  } else if (CurrentSubbank.BankIndex<-1){
    CurrentSubbank.BankIndex=2;
  } else {
    CurrentSubbank.BankIndex+= (int) ROCData[CurrentSubbank.BankIndex]+1; //use subbank word count
  }
  if(CurrentSubbank.BankIndex>ROCWords-2){
    CurrentSubbank.BankIndex=-1;   //flag for no more subbanks in this ROC event
    CurrentSubbank.DataWords=0;
    CurrentSubbank.Header=0;
    CurrentSubbank.Data=NULL;
    CurrentSubbank.Tag=-1;
    return CurrentSubbank;
  }
  
  int Index=CurrentSubbank.BankIndex;
  CurrentSubbank.DataWords=ROCData[Index]-1;  //# of datawords 
  CurrentSubbank.Header=ROCData[Index+1]-1;  //# of datawords 
  CurrentSubbank.Data=&ROCData[Index+2];  //pointer to data in subbank
  CurrentSubbank.Tag=((CurrentSubbank.Header>>16)&0xffff);

  if(verbose>0){
    printf("Subbank Dump: Tag= %10x  NumWrds= %d\n",
	   CurrentSubbank.Tag,CurrentSubbank.DataWords);
    if(verbose>1){
      for (int i=0; i<min(verbose,CurrentSubbank.DataWords); i++){
	int j=Index+i;
	printf(" %5d      %10d   %8x \n",j,ROCData[j],ROCData[j]);
      }
    }
  }
  return CurrentSubbank;
}

  int THaCodaFile::filterToFile(TString output_file) {
// A call to filterToFile filters from present file to output_file
// using filter criteria defined by evtypes, evlist, and max_to_filt 
// which are loaded by public methods of this class.  If no conditions 
// were loaded, it makes a copy of the input file (i.e. no filtering).

       int i;
       if(output_file == filename) {
	 if(CODA_VERBOSE) {
           cout << "filterToFile: ERROR: ";
           cout << "Input and output files cannot be same " << endl;
           cout << "This is to protect you against overwriting data" << endl;
         }
         return CODA_ERROR;
       }
       FILE *fp;
       if ((fp = fopen(output_file.Data(),"r")) != NULL) {
          if(CODA_VERBOSE) {
  	    cout << "filterToFile:  ERROR:  ";
            cout << "Output file `" << output_file << "' exists " << endl;
            cout << "You must remove it by hand first. " << endl;
            cout << "This forces you to think and not overwrite data." << endl;
	  }
          fclose(fp);
          return CODA_ERROR;
       }
       THaCodaFile* fout = new THaCodaFile(output_file.Data(),"w"); 
       int nfilt = 0;

       while (codaRead() == S_SUCCESS) {
           uint32_t* rawbuff = getEvBuffer();
           int evtype = rawbuff[1]>>16;
           int evnum = rawbuff[4];
           int oktofilt = 1;
           if (CODA_DEBUG) { 
	     cout << "Input evtype " << dec << evtype;
             cout << "  evnum " << evnum << endl; 
             cout << "max_to_filt = " << max_to_filt << endl;
             cout << "evtype size = " << evtypes[0] << endl;
             cout << "evlist size = " << evlist[0] << endl;
	   }
           if ( evtypes[0] > 0 ) {
               oktofilt = 0;
               for (i=1; i<=evtypes[0]; i++) {               
                   if (evtype == evtypes[i]) {
                       oktofilt = 1;
                       goto Cont1;
	 	   }
	       }
	   }
Cont1:
           if ( evlist[0] > 0 ) {
               oktofilt = 0;
               for (i=1; i<=evlist[0]; i++) {               
                   if (evnum == evlist[i]) {
                       oktofilt = 1;
                       goto Cont2;
		   }
	       }
	   }
Cont2:
	   if (oktofilt) {
             nfilt++;
             if (CODA_DEBUG) {
	       cout << "Filtering event, nfilt " << dec << nfilt << endl;
	     }
	     int status = fout->codaWrite(getEvBuffer());
             if (status != S_SUCCESS) {
	       if (CODA_VERBOSE) {
		 cout << "Error in filterToFile ! " << endl;
                 cout << "codaWrite returned status " << status << endl;
	       }
               goto Finish;
	     }
             if (max_to_filt > 0) {
    	        if (nfilt == max_to_filt) {
                  goto Finish;
	        }
	     }
	   }
       }
Finish:
       delete fout;
       return S_SUCCESS;
  };



  void THaCodaFile::addEvTypeFilt(int evtype_to_filt)
// Function to set up filtering by event type
  {
     int i;
     initFilter();
     if (evtypes[0] >= maxftype-1) {
        TArrayI temp = evtypes;
        maxftype = maxftype + 100;
        evtypes.Set(maxftype);
        for (i=0; i<=temp[0]; i++) evtypes[i]=temp[i];
        temp.~TArrayI();
     }
     evtypes[0] = evtypes[0] + 1;  // 0th element = num elements in list
     int n = evtypes[0];
     evtypes[n] = evtype_to_filt;
     return;
  };


  void THaCodaFile::addEvListFilt(int event_num_to_filt)
// Function to set up filtering by list of event numbers
  {
     int i;
     initFilter();
     if (evlist[0] >= maxflist-1) {
        TArrayI temp = evlist;
        maxflist = maxflist + 100;
        evlist.Set(maxflist);
        for (i=0; i<=temp[0]; i++) evlist[i]=temp[i];
        temp.~TArrayI();
     }
     evlist[0] = evlist[0] + 1;  // 0th element = num elements in list
     int n = evlist[0];
     evlist[n] = event_num_to_filt;
     return;
  };

  void THaCodaFile::setMaxEvFilt(int max_event)
// Function to set up the max number of events to filter
  {
     max_to_filt = max_event;
     return;
  };



void THaCodaFile::staterr(TString tried_to, int status) {
// staterr gives the non-expert user a reasonable clue
// of what the status returns from evio mean.
// Note: severe errors can cause job to exit(0)
// and the user has to pay attention to why.
    if (status == S_SUCCESS) return;  // everything is fine.
    if (tried_to == "open") {
       cout << "THaCodaFile: ERROR opening file = " << filename << endl;
       cout << "Most likely errors are: " << endl;
       cout << "   1.  You mistyped the name of file ?" << endl;
       cout << "   2.  The file has length zero ? " << endl;
       cout << "Job aborting !! " << endl;
       exit(0);
    }
    switch (status) {
      case S_EVFILE_TRUNC :
	 cout << "THaCodaFile ERROR:  Truncated event on file read" << endl;
         cout << "Evbuffer size is too small.  Job aborted." << endl;
         exit(0);  //  If this ever happens, recompile with MAXEVLEN
	           //  bigger, and mutter under your breath at the author.    
      case S_EVFILE_BADBLOCK : 
        cout << "Bad block number encountered " << endl;
        break;
      case S_EVFILE_BADHANDLE :
        cout << "Bad handle (file/stream not open) " << endl;
        break;
      case S_EVFILE_ALLOCFAIL : 
        cout << "Failed to allocate event I/O" << endl;
        break;
      case S_EVFILE_BADFILE :
        cout << "File format error" << endl;
        break;
      case S_EVFILE_UNKOPTION :
        cout << "Unknown option specified" << endl;
        break;
      case S_EVFILE_UNXPTDEOF :
        cout << "Unexpected end of file while reading event" << endl;
        break;
      case EOF: 
        if(CODA_VERBOSE) {
          cout << "Normal end of file " << filename << " encountered" << endl;
	}
        break;
      default:
        cout << "Error status  0x" << hex << status << endl;
      }
  };

  void THaCodaFile::init(TString fname) {
    handle = 0;
    filename = fname;
  };

  void THaCodaFile::initFilter() {
    if (!ffirst) {
       ffirst = 1;
       maxflist = 100;
       maxftype = 100;
       evlist.Set(maxflist);
       evtypes.Set(maxftype);
       evlist[0] = 0;
       evtypes[0] = 0;
    }
  };












