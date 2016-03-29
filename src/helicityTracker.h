//
//	helicityTracker Class
//      keeps track of helicity
//      create 10/07/2015  G.B. Franklin
//
#include <stdlib.h>
#include <stdint.h>
#include <TROOT.h>
#include "TTree.h"
#include "TMath.h"

#ifndef helicityTracker_h
#define helicityTracker_h

#define HELPLUS 1
#define HELMINUS 0

#ifndef helicityStatus_h
#define helicityStatus_h

struct helicityStatus
{
  int helicityState;
  int currentHelicityBit;
  int predictedCurrentBit;
  int indexPair;
  int indexQuartet;
  int indexOctet;
  int countPairs;
  int countQuartets;
  int countOctets;
  int countHWs;
  int countSynchedHWs;
  bool helicityValid;
  bool helicitySynchValid;
  bool helicityPredictionValid;
  bool helicityBitInAgreement;
};
#endif

class helicityTracker {
 public:
  helicityTracker();
  void SetHelicityMode(int helicityDelay, int helicityShiftBits,
		       int helicityStructure, int helicityBitFlip);
  int newRun();  //called at start of run to initizlize status variables  
  int newHelicityWindow(int currentHelicityBit,int helicitySetSyncBit, int mpsScaler,
			int verbose); //call once per helicity window
  //
  // newHelicityWindow returns 0 or 1 for valid helicity, -1 for invalid
  // usecall with  HelicitySetSynchBit=-1 for automatic synching

  // use follow routines to retriev additional info form statusHW structure
  helicityStatus* GetHelicityStatus(){return statusHW;}; //point to helicity status info
  void DebugDump(int mode);
  int GetHelicityState();  //returns -1 if not valid, otherwise 0 or 1 
  int GetCountHelicityWindows(){return statusHW->countHWs;};
  int GetCountConsecutiveGoodWindows(){return statusHW->countSynchedHWs;};
  int GetCountPairs(){return statusHW->countPairs;};
  int GetIndexPair(){return statusHW->indexPair;}; //index within a quartet
  int GetCountQuartets(){return statusHW->countQuartets;};
  int GetIndexQuartets(){return statusHW->indexQuartet;}; //index within a quartet
  int GetCountOctets(){return statusHW->countOctets;};
  int GetIndexOctet(){return statusHW->indexOctet;}; //index within a quartet
  bool IsHelicityValid(){return statusHW->helicityValid;};
  bool IsHelicityPredictionValid(){return statusHW->helicityPredictionValid;};
  bool IsHelicityBitInAgreement(){return statusHW->helicityBitInAgreement;};
  bool IsHelicitySynchValid(){return statusHW->helicitySynchValid;};

  //******************************
 private:
  bool newHelicitySet(int currentBit);
  int updateShiftRegister(int hRead,int registerType);
  //
  // internal variables
  helicityStatus* statusHW;  //will point to helicity window info
  //
  //helicity parameters filled in at newRun
  int helicityStructure; //#of windows in each set  (4 for "quads", 8 for Octs)
  int helicityShiftBits;  //# bits in helicity generating shift register
  int helicityDelay;  //filled in at start of run using comptomParams
  int helicityPredictorBitFlip;
  int helicityRegisterType;  
  
  //internal helcity tracking stuff
  int patternIndex;  //0 to 3 if doing quartets, 0 to 7 if doing octets
  uint32_t helBitHistory;  //keep track of actual bits read
  uint32_t helPredHistory;  //actual helicity history starting 16 MPS in the future
  uint32_t fgShreg;     //value for helicity sequence algorithm
  int fgNShreg;    //count since fgShreg was reset
  //current helicity set patterns for Bit Readout and for Actual Helicity
  //(same thing if not running with helicty bit delayed)
  uint32_t currentBitInHistory;  //bit corresponding to current Bit
  bool helicitySetTransitionFound;
  //debugging and resynching info
  int currentMPSscaler;
  int previousMPSscaler;
};
#endif


