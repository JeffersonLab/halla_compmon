/*
//	fadcTriggered
//      Analysis of triggered fadc dada
*/

#include <stdlib.h>
#include <TROOT.h>
#include "TTree.h"
#include "TMath.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "THaEpics.h"
#include "THaCodaFile.h"
#include "bankstructure.h"
//Compton fadc-specific classes
#include "comptonStatus.h"
#include "textParams.h"
#include "fadcdata.h"
#include "vmeauxdata.h"

#ifndef fadcTriggeredSums_h
#define fadcTriggeredSums_h

#define MAX_DAC_VALUES 20
#define HELPLUS 1
#define HELMINUS 0

// laserstate now defined in class comptonStatus
//enum laserstate {RIGHT, LEFT, OFF, UNKNOWN};
//enum laserstate {RIGHT, LEFT, RIGHTOFF, LEFTOFF, UNKNOWN};

class fadcTriggered {
 public:
  fadcTriggered();
  fadcTriggered(textParams* theParams, comptonStatus* theStatus);
  void newRun();   //init counters at start of a run
  int DefineHistos();
  int DefineTriggeredTree();
  int DoSummedPulses(vmeauxdata* theVMEauxdata,fadcdata *theFADCdata); 
  int DoNormalizedHistos();
  int DoSampledWaveforms(fadcdata *theFADCdata); 
 private:
  textParams* theParams;   //pointer to parameter class
  comptonStatus* theStatus;  //poitner to the status class
  // parmameters from textParams class
  int channel_calorimeter_PMT;
  Float_t ped_value;
  void labelSpinSortedHistos(TH1F* histo); //used to label x-axis of laser-sorted histos
  int NumSamples;  //Filled with number of valid waveform samples (snapshots)
  //
  //Variables output to triggeredWise root tree
  TTree* triggerWiseTree;
  Float_t sumVal;   //fadc triggered data summed over single pulse
  Float_t sumClock;   //fadc triggered data clock time
  Float_t sumPedestal; ///fadc triggered data pedestal measurement (crl after 10/1016)
  Bool_t  sumIsRandom; ///fadc triggered data (is this sum random trigger?)
  //
  // Variables for building pulserWise root tree (MiniMegan pulser info)
  //one entry for each set of 4 MiniMegan settings
  TTree* pulserWiseTree;
  Float_t MMpulse[4]; 
  int MMVarDacIndex;
  int MMVarDacSetting;
  //
  // special tree for snapshots  (sampled waveforms)
  TTree* snapshotsTree;
  int mpsCount;
  int mpsLaserOnCount;
  int mpsLaserOffCount;
  Float_t bcmLaserOnSum;
  Float_t bcmLaserOffSum;
  int helicityState;
  int laserState;
  int snapshotClock;
  Float_t snapshot[1000]; //will use this for a snapshot array
  int randomTime;// 1-> snapsshots israndom time samples
  //
  // Triggered Sums data histograms
  TH1F* hTrig_numSums;  //#triggers actually summed each MPS
  TH1F* hTrig_sums_All;  //histo of all summed pulses
  TH1F* hTrig_sums_laserOn;  //histo of laser-on summed pulses
  TH1F* hTrig_sums_laserOff;  //histo of laser-on summed pulses
  TH1F* hTrig_sums[4];  // pulse sums for Laser Left, Positive Helicity Beam. etc

  // calculated periodically from other histograms 
 //Normalized histos
  TH1F* hNorm_Trigs_Scaler;   //triggers (via scaler) per MPS
  TH1F* hNorm_sums_subtracted; //background subtracted (laserOn-laserOff)
  TH1F* hNorm_sums_asym;   //asymmetry vs energy
  //Triggered Sampled waveform histograms
  TH1F* hTrig_numSamples;
  TH1F* hTrig_wf;       //sum (average) of snapshots
  TH1F* hTrig_wf_Ped;   //gather pedestal information from snapshots
};

#endif
