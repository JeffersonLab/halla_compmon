//*******************************
//    vmeauxdata.cc
//    created 11/15/2008 D. Parno
//    handles VME data from auxiliary modules:
//		HAPPEX timing board
//		Trigger Interface Register
//		16-channel Caen v560 scaler
//*******************************

#include <iostream>
using namespace std;
#include "vmeauxdata.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
vmeauxdata::vmeauxdata(){               // Constructor
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
vmeauxdata::vmeauxdata(textParams* theParamsIn){   // Constructor
  theParams=theParamsIn;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void vmeauxdata::newRun(){
  helBit = -1;
  mpsSignal = 0;
  qrt = 0;
  cavpowBit = 0;
  cavpolBit = 0;
  evlen = 0;
  rampDelay = 0;
  intTime = 0;
  DACsetting =0;
  for (int i = 0; i < 17; i ++){
    scalers[i] = 0;
    IPscalers[i]=0;
  }
  chan_runBCM=theParams->getInt("chan_runBCM");
  chan_runClock=theParams->getInt("chan_runClock");
  chan_runL1A=theParams->getInt("chan_runL1A" );
  chan_runMPS=theParams->getInt("chan_runMPS" );
  chan_runTriggers=theParams->getInt("chan_runTriggers" );
  chan_ipBCM=theParams->getInt("chan_ipBCM" );
  chan_ipClock=theParams->getInt("chan_ipClock" );
  chan_ipBPM2Aym=theParams->getInt("chan_ipBPM2Aym" );
  chan_ipBPM2Ayp=theParams->getInt("chan_ipBPM2Ayp" );
  chan_ipBPM2Axm=theParams->getInt("chan_ipBPM2Axm" );
  chan_ipBPM2Axp=theParams->getInt("chan_ipBPM2Axp" );

  chan_ipBPM2Bym=theParams->getInt("chan_ipBPM2Bym" );
  chan_ipBPM2Byp=theParams->getInt("chan_ipBPM2Byp" );
  chan_ipBPM2Bxm=theParams->getInt("chan_ipBPM2Bxm" );
  chan_ipBPM2Bxp=theParams->getInt("chan_ipBPM2Bxp" );
  
  chan_ipPowLeft=theParams->getInt("chan_ipPowLeft" );
  chan_ipPowRight =theParams->getInt("chan_ipPowRight" );
  chan_ipCavPower=theParams->getInt("chan_ipCavPower" );
  chan_ipTriggers=theParams->getInt("chan_ipTriggers" );
  chan_ipVFinger=theParams->getInt("chan_ipVFinger" );
  chan_ipHFinger=theParams->getInt("chan_ipHFinger" );

  bitflip_cavpol=theParams->getInt("bitflip_cavpol" );
  bitflip_helicity=theParams->getInt("bitflip_helicity");
  return;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
int vmeauxdata::Unpack(bankstructure bank){
  unsigned int* data = (unsigned int*)bank.Data;

  // Timing Board Data
  // data[0]=1  for version number (ignore for now)
  numAuxWords=(int)data[1];
  rampDelay = (int)data[2];
  intTime = (int)data[3];
  DACsetting= (int)data[4]; //DAC setting used for LED pulser scan
  // TIR data
  int mytirdata = data[5];
  dithering=(mytirdata & 0x100)>>8;     // channel 10? -- DOUBLE CHECK THIS
  helBit = (mytirdata & 0x10)>>4;
  if(bitflip_helicity!=0) helBit=1-helBit;
  cavpolBit = (mytirdata & 0x20)>>5;
  if(bitflip_cavpol!=0) cavpolBit=1-cavpolBit;
  cavpowBit = (mytirdata & 0x40)>>6;
  mpsSignal = (mytirdata & 0x80)>>7;

  int point = numAuxWords+3;
  int point2 = point+data[point-1]+2;

  // Scaler data
  //if ((data[3]&0xcae56000) == 0xcae56000) {
  if(data[point] == 0xcae56000){
    for (int i = 0; i < 16; i ++){
      scalers[i] = data[point+i+1];
    }
  }
  else {              // Something has changed in the subbank
    cout << "\nERROR in vmeauxdata::Unpack: Scaler header does not appear where expected.\n";
  }

  if (data[point2] == 0xdae56000) {
    for (int i = 0; i < 16; i ++){
      IPscalers[i] = data[point2+1+i];
    }
  }
  else {              // Something has changed in the subbank
    cout << "\nERROR in vmeauxdata::Unpack: IP Scaler header does not appear where expected.\n";
  }


  //	DebugDump(2);

  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Which module information to dump?    0 = all
//                                      Multiples of 2 = Trigger Interface
//                                      Multiples of 3 = Timing Board
//                                      Multiples of 5 = Scalers
// e.g. 6 gives you both trigger interface and timing board
int vmeauxdata::DebugDump(int whichmodule){
  cout<< "vmeauxdata DebugDump**************************************"<<endl;
  if (whichmodule%2 == 0){
    cout << "\nTIR data: helicity =            " << helBit;
    cout << "\n          cavity power =        " << cavpowBit;
    cout << "\n          cavity polarization = " << cavpolBit;
    cout << "\n          event length =        " << evlen;
    cout << "\n          mpsSignal bit =       " << mpsSignal << endl;
    cout << "\nMapped Run Period Scaler Channels";
    cout << "\n  MPS    =                        " << GetMPSScaler();
    cout << "\n  BCM    =                        " << GetBCMScaler();
    cout << "\n  Clock   =                       " << GetClockScaler();
    cout << "\n  Cal. PMT =                      " << GetPMTScaler();
    cout << "\nMapped Integration Period Scaler Channels";
    cout << "\n  Clock  =                        "<<GetGatedClock();

    cout << "\n  BCM    =                        " << GetGatedBCM();
    cout << "\n  CavityPower                     " << GetCavityPowerScaler();
    cout << "\n  Triggers=                       " << GetTriggerIPScaler();
    cout <<endl;
    //		if (cavpow>1) cout << "\nHallelujah! Cavity was on!";
  }

  if (whichmodule%3 == 0){
    cout << "\nTiming Board data: ramp delay = " << rampDelay;
    cout << "   integration time = " << intTime << endl;
  }

  if (whichmodule%5 == 0){
    printf("     Run Scaler         IP Scaler \n");
    for (int i = 0; i < 16; i++){
      printf(" %2d: %10d   %2d:  %10d\n",i,scalers[i],i,IPscalers[i]);
    }
  }
  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
