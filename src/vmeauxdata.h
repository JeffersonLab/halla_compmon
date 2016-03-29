#include <stdio.h>
#include <TROOT.h>
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "bankstructure.h"
#include "textParams.h"

#ifndef vmeaux_h
#define vmeaux_h

// Add scaler data?
class vmeauxdata {
  private:
  textParams* theParams; // pointer to parameter handling class
  // header data
  int numAuxWords;
    // TIR data
  int helBit;   //helicity bit
  int dithering; //dithering on bit
    int qrt;
    int cavpowBit;  //cavity power Bit (but also look at real-time cav power
    int cavpolBit;
    int evlen;
    int mpsSignal;

    // HAPPEX Timing Board data
    int rampDelay;
    int intTime;
    int DACsetting;   //DAC setting for pulser
    //
    unsigned int scalers[17];
    unsigned int IPscalers[17];

    //channel mapping for scalers
    int chan_runBCM;
    int chan_runClock;
    int chan_runL1A;
    int chan_runMPS;
    int chan_runTriggers;
    int chan_ipBCM;
    int chan_ipClock;
    int chan_ipBPM2Aym;
    int chan_ipBPM2Ayp;
    int chan_ipBPM2Axm;
    int chan_ipBPM2Axp;
    int chan_ipBPM2Bym;
    int chan_ipBPM2Byp;
    int chan_ipBPM2Bxm;
    int chan_ipBPM2Bxp;
    int chan_ipPowLeft;
    int chan_ipPowRight;
    int chan_ipCavPower;
    int chan_ipTriggers;
    int chan_ipVFinger;
    int chan_ipHFinger;
    //allow for bit flip mistakes in dAQ
    int bitflip_cavpol;
    int bitflip_helicity;

  public:

    vmeauxdata();		        	// Constructor
    vmeauxdata(textParams* theParams);  // Constructor
    void newRun();

    //Bits
    int GetMPSSignal()		{return mpsSignal;};
    int GetHelicityBit()	{return helBit;};
    int GetDithering()  {return dithering;};
    int GetCavityPowerBit()	{return cavpowBit;};
    int GetCavityPolBit()	{return cavpolBit;};

    //Timing Board Info
    int GetRampDelay()  	{return rampDelay;};
    int GetIntTime()    	{return intTime;};
    int GetDACSetting()         {return DACsetting;};
    //
    // Run Period Scaler
    unsigned int GetScaler(int i){return scalers[i];};

    unsigned int GetBCMScaler()  {return scalers[chan_runBCM];};
    unsigned int GetClockScaler(){return scalers[chan_runClock];};
    unsigned int GetL1A() 	{return scalers[chan_runL1A];};
    unsigned int GetMPSScaler()	{return scalers[chan_runMPS];};
    unsigned int GetPMTScaler() {return scalers[chan_runTriggers];};


    //    unsigned int GetTriggerScaler(){return scalers[6];};
    
    //
    //Integration Period Scalers  (note erros with multiple references
    unsigned int GetIPScaler(int i){return IPscalers[i];};

    unsigned int GetGatedBCM(){return IPscalers[chan_ipBCM];};
    unsigned int GetGatedClock(){return IPscalers[chan_ipClock];};
    unsigned int GetVtoFBPM2AymIPScaler(){return IPscalers[chan_ipBPM2Aym];};
    unsigned int GetVtoFBPM2AypIPScaler(){return IPscalers[chan_ipBPM2Ayp];};
    unsigned int GetVtoFBPM2AxmIPScaler(){return IPscalers[chan_ipBPM2Axm];};
    unsigned int GetVtoFBPM2AxpIPScaler(){return IPscalers[chan_ipBPM2Axp];};
    unsigned int GetVtoFBPM2BymIPScaler(){return IPscalers[chan_ipBPM2Bym];};
    unsigned int GetVtoFBPM2BypIPScaler(){return IPscalers[chan_ipBPM2Byp];};
    unsigned int GetVtoFBPM2BxmIPScaler(){return IPscalers[chan_ipBPM2Bxm];};
    unsigned int GetVtoFBPM2BxpIPScaler(){return IPscalers[chan_ipBPM2Bxp];};

    unsigned int GetVtoFPowLeftIPScaler(){return IPscalers[chan_ipPowLeft];};
    unsigned int GetVtoFPowRightIPScaler(){return IPscalers[chan_ipPowRight];};
    unsigned int GetVtoFPowIPScaler(){return IPscalers[chan_ipCavPower];};
    unsigned int GetTriggerScaler(){return IPscalers[chan_ipTriggers];};
    unsigned int GetTriggerIPScaler(){return IPscalers[chan_ipTriggers];};

    unsigned int GetVerticalFingerIPScaler(){return IPscalers[chan_ipVFinger];};
    unsigned int GetHorizontalFingerIPScaler(){return IPscalers[chan_ipHFinger];};

    //redundant?
    unsigned int GetTest5(){return IPscalers[5];};
    unsigned int GetTest6(){return IPscalers[6];};
    unsigned int GetCavityPowerScaler(){return IPscalers[chan_ipCavPower];};

    int Unpack(bankstructure bank);

    // Which module information to dump?	0 = all
    //	  		    Multiples of 2 = Trigger Interface and mapped scalers
    //			    Multiples of 3 = Timing Board
    //			    Multiples of 5 = Scalers
    // e.g. 6 gives you both trigger interface and timing board
    int DebugDump(int whichmodule=0);
};

#endif 	//vmeaux_h
