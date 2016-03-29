//*****************************
//*  helicityTracker.cc
//*  created 10/01/2015 G.B. Franklin
//*  Designed for keeping track of Compton helicity setsa
//*
//*****************************
#include <iostream>
using namespace std;
#include "helicityTracker.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
helicityTracker::helicityTracker(){
  //constructor
  helicityShiftBits=0;  //use as flag for unitialized helicitdy parameters
  statusHW=new helicityStatus;  //helicityStatus structure defined
  return;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void helicityTracker::SetHelicityMode(int helicityDelayIn,
				      int helicityShiftBitsIn,
				 int helicityStructureIn, int helicityBitFlipIn){
  helicityDelay=helicityDelayIn;        //0 for no delay in helicity bit, etc.
  helicityShiftBits=helicityShiftBitsIn;//number of bits in shift register
  helicityStructure=helicityStructureIn;  // 2,4,or 8 helicity windows in a set
  helicityPredictorBitFlip=helicityBitFlipIn;   //special flip bit in random registers
  //                 (appears to be needed with some Helicity Control Boards)
  return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int helicityTracker::newRun(){
  //Init all counters and status
  statusHW->helicityState=-1;
  statusHW->indexPair=-1;  //where are we within a helicity pair?
  statusHW->indexQuartet=-1; //etc.
  statusHW->indexOctet=-1;
  statusHW->countPairs=0; //zero pair counter.
  statusHW->countQuartets=0; //etc.
  statusHW->countOctets=0;
  statusHW->countHWs=0;  //increment for eahch newHelicityWindow call
  statusHW->countSynchedHWs=0;  //HWs since helicty sets are synched
  statusHW->helicityValid=false;
  statusHW->helicitySynchValid=false;
  statusHW->helicityPredictionValid=false;
  statusHW->helicityBitInAgreement=false;
  currentMPSscaler=-1;
  previousMPSscaler=-1;
  //
  //stuff for synching predictions of helicity
  patternIndex=-1;     //used by helicty decoder
  fgShreg=0;   //helcity register decoding
  fgNShreg=0;   //counts since fgShreg reset
  helPredHistory=0;  //shift register for actual helicity bit read
  helBitHistory=0;   //history of helicity bits actualy read
  currentBitInHistory=0x80;  //Bit inhelPredHistory holding info for helicity
  currentBitInHistory= currentBitInHistory<<helicityDelay;
  if(helicityShiftBits==30){
    helicityRegisterType=0;
  }else if (helicityShiftBits==24){
    helicityRegisterType=1;
  }else{
    helicityRegisterType=-1;
    printf("FATAL ERROR:  helicityShiftBits=%d Unknown\n",helicityShiftBits);
  }
  //
  //
  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
int helicityTracker::GetHelicityState(){
  if(statusHW->helicityValid) {
    return statusHW->helicityState;
  }else{
    return -1;
  }
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void helicityTracker::DebugDump(int mode){
  // Dump Status info for debugging
  //  mode==0   lots of info
  if(mode==1){
    printf("helicityTracker Dump \n");
  }
  return ;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
int helicityTracker::newHelicityWindow(int currentBit,int helicitySetSynchBit,
				       int mpsScaler,int verbose){
  //update helicity state
  // call once for each new helicity window
  // currentBit and helicitySetSynchBit should be signals from helicity board
  // OR if helicitySetSynchBit=-1, the first newHelicityWindow 
  // calls used are to find transition between two setsof helicity window
  // WARNING helcitySetSynchBit not yet implemented  (always does an autosynch)
  //      (patternIndex set to nonzero) and helicitySetTransitionFound set to true
  // next calls used to get to start of next helicity window
  //    helicitySynchValid set to true at start of first synched helicity set
  // one helcitySynchvalid==true, newHelicitySet is called at start of every helicity set
  // next 30 helicity sets  (240 windows if running octs) used to fill shiftregister
  //     helicityValid set to true once shiftregister if filled and the future can be predited
  // 
  // uses first events to  train shift register,
  // then predicts future helicity bit patterns if in delayed reporting mode
  int setCheck;
  bool validPredictionRegister;
  uint32_t tmp;
  //
  // Update history of helicity predictins and actual currentHelicitbit
  helPredHistory=helPredHistory<<1;            //history of predictions
  helBitHistory= ( helBitHistory<<1 | currentBit); //history of currentBit 
  statusHW->countHWs++;  //count calls to newHelicityWindow
  statusHW->currentHelicityBit=currentBit;  //same has helicity if running with no delayed reporting
  //
  //  If not already synched, see if we are at a transitin between two
  // helicity sets by looking at recent history of helicity bit
  //
  previousMPSscaler=currentMPSscaler;
  currentMPSscaler=mpsScaler;
  if( !helicitySetTransitionFound){
    if(helicityStructure==8&&statusHW->countHWs>8){
      tmp=helBitHistory& 0xFF;  //look at lower 8 bits of inut bit history
      if(tmp==0xaa || tmp==0x66) {
	//OctetA to Octet A transition or OctetB to OctetB transition
	patternIndex=2;  //halfway into an octet
	helicitySetTransitionFound=true;
	printf("Synched to octet sets at Helicity Window %d \n",
	       statusHW->countHWs);
      }
    }else if(helicityStructure==4&&statusHW->countHWs>4){
      tmp=helBitHistory& 0xF;  //look at lower 4 bits of history
      if(tmp==0xa || tmp==0x5) {
	//QuadA to Quad A transition or QuadB to QuadB transition
	patternIndex=0;  //end of a quad
	helicitySetTransitionFound=true;
	printf("Synched to quartet sets at Helicity Window %d pattern 0x%8x\n",
	       statusHW->countHWs,helBitHistory);
      }
    }else if(helicityStructure==2&&statusHW->countHWs>2){
      tmp=helBitHistory& 0x3;  //look at lower 2 bits of history
      if(tmp==0x0 || tmp==0x3) {
	//PairA to PairA transition or PairB to PairB transition
	patternIndex=1;  //halfway into an pair
	helicitySetTransitionFound=true;
	printf("Synched to pair sets at Helicity Window %d \n",
	       statusHW->countHWs);
      }
    }
    statusHW->helicitySynchValid=false;//begining of a helicity set not yet found
    statusHW->helicityPredictionValid=false;//not yet able to make predictions
  }
  if(patternIndex>=0) {
    patternIndex++;
    if(patternIndex>=helicityStructure) patternIndex=0;  //cycle pattern index
    if(patternIndex==0){
      statusHW->helicitySynchValid=true;   //found start of a helicity set
      validPredictionRegister=newHelicitySet(currentBit);  //update for start of next set
      if(validPredictionRegister){
	statusHW->helicityPredictionValid=true;
      }else{
	statusHW->helicityPredictionValid=false;
      }
    }
  }
  //
  // now wrap up bookeeping for each new helicity window
  if(statusHW->helicitySynchValid){
    statusHW->countSynchedHWs++;
  }else{
    statusHW->countSynchedHWs=-1;
  }

  if(statusHW->helicityPredictionValid){
    statusHW->helicityState= (helPredHistory & 0x80) !=0;   //predicted helicity
    setCheck= (helPredHistory & currentBitInHistory) !=0; //predicted bit
    statusHW->predictedCurrentBit=setCheck;    //predicted current bit
    statusHW->helicityBitInAgreement= (setCheck==currentBit); //do they agree?
    statusHW->helicityValid=true;  //calling routine delas with setCceck!=currentBit
  } else if (helicityDelay==0) {
    statusHW->helicityState=currentBit;  //override prediction if no delay
    statusHW->helicityValid=true;
  }else{
    statusHW->helicityState=-1;
    statusHW->helicityValid=false;
  }
  //now compute index for pairs, quartets, and octets
  if(statusHW->helicityValid && statusHW->helicitySynchValid){
    int tmp=statusHW->countSynchedHWs;
    statusHW->indexOctet =tmp & 0x7;
    statusHW->indexQuartet = tmp & 0x3;
    statusHW->indexPair = tmp & 0x1;
    if(statusHW->indexOctet==0) statusHW->countOctets++;
    if(statusHW->indexQuartet==0) statusHW->countQuartets++;
    if(statusHW->indexPair==0) statusHW->countPairs++;
  }else{
    statusHW->indexOctet = -1;  //flag for not valid Octet index (not synched)
    statusHW->indexQuartet = -1;
    statusHW->indexPair = -1;
  }
  if(verbose==2){
    if(statusHW->helicityPredictionValid){
      if(patternIndex==0) printf("\n");
      printf(" %8d %4d Quad %d currentBit /predBit  Helicity  %d / %d    %d    0x%06x",
	     mpsScaler,patternIndex, statusHW->indexQuartet,currentBit,
	     setCheck,statusHW->helicityState,helPredHistory);
      if(setCheck!=currentBit){
	printf("   Error\n");
      }else{
	printf("\n");
      }

    }else {
      if(patternIndex==0) printf("\n");
      printf(" %8d %4d Quad %d currentBit /predBit  Helicity  %d / NA    %d    0x%06x\n",
	     mpsScaler,patternIndex, statusHW->indexQuartet,currentBit,
	     statusHW->helicityState,helPredHistory&0x3FF);
    }
  }else if(verbose==1){
    if(statusHW->helicityPredictionValid){
      if(setCheck!=currentBit) {
	printf(" %8d %4d Quad %d currentBit /predBit  Helicity  %d / %d    %d    0x%06x  Helicity Pattern Error \n",
	       mpsScaler,patternIndex, statusHW->indexQuartet,currentBit,
	       setCheck,statusHW->helicityState,helPredHistory&0x3FF);
      }
    }
  }
  return statusHW->helicityState;
} 


//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
bool helicityTracker::newHelicitySet(int currentBit){
  //  called at start of a set of helicity windows
  // Update shiftregister if at the start of a quad
  // if fgNShreg<30, we are still training the shift regiseter
  //  when fgNShreg==30, we update shiftregister to get prediction'
  // of first bit of future helicity patternsa
  //  afgShreg used by "updateShiftRegister" to generatre psuedo-randoms
  //  fgNShreg used to count training bits (need 30 to fill register
  //
  int predictedHelicitySet;
  unsigned int helPattern;
  if(fgNShreg<30){     //are we still training the shift register?
    helPredHistory=0;  //no valid history at this time
    predictedHelicitySet =
      updateShiftRegister(currentBit,helicityRegisterType); //yes, feed it the bit
    helPredHistory=predictedHelicitySet;   //first bit of predictionhistory
    fgNShreg++;
    if(fgNShreg==30) {
      //finished training.  jump forward in in predictions and
      //also record bit within shiftregister for "currentBit" preddiciton"
      if(currentBit==1){
	helPattern=0x96;
      }else{
	helPattern=0x69;
      }
      helPredHistory=helPattern;
      if(helicityDelay>0){
	for(int i=0; i<helicityDelay/helicityStructure; i++){
	  predictedHelicitySet = updateShiftRegister(2,helicityRegisterType);
	  if(predictedHelicitySet==1){
	    helPattern=0x96;
	  }else{
	    helPattern=0x69;
	  }
	  helPredHistory=helPredHistory<<helicityStructure;
	  helPredHistory=(helPredHistory&0xFFFFFF00)|helPattern;
	}
      }
      printf("Helicity Predictions established at Helicity Window %d\n",
	     statusHW->countHWs);
    }
  }else{
    predictedHelicitySet = updateShiftRegister(2,helicityRegisterType); //no, let it decide
    if(predictedHelicitySet==1){
      helPattern=0x96;
    }else{
      helPattern=0x69;
    }
    helPredHistory=(helPredHistory&0xFFFFFF00) + helPattern;   //add next set of predictions 
  }
  if(fgNShreg>=30){
    return true;
  }else{
    return false;
  }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
int helicityTracker::updateShiftRegister(int hRead,int registerType){
  // updates helcity-tracking shifdt-register used for pseudo-randoms
  //  input hRead= current helicty bit read for training registrer
  //  input hRead= 2 to return prediction of next bit
  uint32_t newbit;
  uint32_t inputBit;
  if(helicityPredictorBitFlip==1){
    inputBit=1-hRead;
  }else{
    inputBit=hRead;
  }
  if(registerType==0){
    //30 bit register
    uint32_t bit7    = (fgShreg & 0x00000040) != 0;
    uint32_t bit28   = (fgShreg & 0x08000000) != 0;
    uint32_t bit29   = (fgShreg & 0x10000000) != 0;
    uint32_t bit30   = (fgShreg & 0x20000000) != 0;
    newbit = (((bit30 ^ bit29) ^ bit28) ^ bit7) & 0x1;
    fgShreg = ( (hRead == 2 ? newbit : inputBit) | (fgShreg << 1 )) & 0x3FFFFFFF;
  }else{
    //24 bit register
    const uint32_t IB1 = 0x1;           // Bit 1 mask
    const uint32_t IB3 = 0x4;           // Bit 3 mask
    const uint32_t IB4 = 0x8;           // Bit 4 mask
    const uint32_t IB24 = 0x800000;         // Bit 24 mask
    const uint32_t MASK = IB1+IB3+IB4+IB24; // 100000000000000000001101
    newbit = (fgShreg & IB24) ? 1 : 0;
    if ((hRead == 2 ? newbit : inputBit) == 1)
      fgShreg = ((fgShreg ^ MASK) << 1) | IB1;
    else
      fgShreg <<= 1;
  }
  if(helicityPredictorBitFlip==1) newbit=1-newbit;
  return newbit;
}
