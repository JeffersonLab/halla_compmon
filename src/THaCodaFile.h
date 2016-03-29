#ifndef THaCodaFile_h
#define THaCodaFile_h

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
/////////////////////////////////////////////////////////////////////

#include "Rtypes.h"
#include "THaCodaData.h"
#include <stdlib.h>
#include "evio.h"
#include "TString.h"
#include "TArrayI.h"
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "bankstructure.h"

class THaCodaFile : public THaCodaData 
{

public:

  THaCodaFile();


  THaCodaFile(TString filename);
  THaCodaFile(TString filename, TString rw);
  ~THaCodaFile();
  int codaOpen(TString filename);
  int codaOpen(TString filename, TString rw);
  int codaClose();
  int codaRead(); 
  int codaWrite(uint32_t* evbuffer);
  uint32_t *getEvBuffer();     

  int UnpackEventHeader(int verbose);
  uint32_t* getROCdata(int SeekRocNum, int verbose); 
  void evDump(int verbose);
  void evDumpStructure(int verbose);  //dumps structure of Event Type 1
  bankstructure getNextSubbank(int verbose);
  int GetEventNumber(){return EventNumber;};

  int filterToFile(TString output_file);     // filter to an output file
  void addEvTypeFilt(int evtype_to_filt);    // add an event type to list
  void addEvListFilt(int event_to_filt);     // add an event num to list
  void setMaxEvFilt(int max_event);          // max num events to filter

private:

     static const int MAXROC = 32;  
     static const int MAXSLOT = 27;  

     static const int MAX_PHYS_EVTYPE  = 14;  // Types up to this are physics
     static const int SYNC_EVTYPE      = 16;
     static const int PRESTART_EVTYPE  = 17;
     static const int GO_EVTYPE        = 18;
     static const int PAUSE_EVTYPE     = 19;
     static const int END_EVTYPE       = 20;
     static const int TS_PRESCALE_EVTYPE  = 120;
     static const int EPICS_EVTYPE     = 131;
     static const int PRESCALE_EVTYPE  = 133;
     static const int DETMAP_FILE      = 135;
     static const int TRIGGER_FILE     = 136;
     static const int SCALER_EVTYPE    = 140;

  THaCodaFile(const THaCodaFile &fn);
  THaCodaFile& operator=(const THaCodaFile &fn);
  void init(TString fname);
  void initFilter();
  void staterr(TString tried_to, int status);  // Can cause job to exit(0)
  int ffirst;
  int max_to_filt,handle;
  int maxflist,maxftype;
  TArrayI evlist, evtypes;

  //  int EventBuffer[MAX_EVENT_SIZE];
  int PhysRecPoint;
  int PhysRecWords;
  int PhysRecHeader[8];
  int EventWordCount;
  int EventType;
  int EventDataType;
  int EventDataNum;
  int EventNumber;
  int EventClass;
  int EventStatus;
  int ROCHeader[2];
  int ROCWords;
  int RunNumber;
  int StartTime;
  int EndTime;
  int RunOnTime;
  int GoTime;
  int ROCnum;
  uint32_t* ROCData;
  bankstructure CurrentSubbank; //current subbank


#ifndef STANDALONE
  ClassDef(THaCodaFile,0)   //  File of CODA data
#endif

};

#endif






