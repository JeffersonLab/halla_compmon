
/*
//	fadcAccums
//      Analysis of fadc accumulator dada
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
#include "comptonParams.h"
#include "fadcdata.h"
#include "vmeauxdata.h"

#ifndef fadcAccums_h
#define fadcAccums_h

#define MAX_DAC_VALUES 20
#define HELPLUS 1
#define HELMINUS 0

// laserstate now defined in class comptonStatus
//enum laserstate {RIGHT, LEFT, RIGHTOFF, LEFTOFF, UNKNOWN};

class fadcAccums {
 public:
  fadcAccums();
  fadcAccums(comptonParams* theParams);
  void newRun();   //init counters at start of a run
  int DefineHistos();
  int DefineTree(comptonStatus* theComptonStatus);
  int DoAccums(vmeauxdata* theVMEauxdata,
		     fadcdata *theFADCdata, comptonStatus* theComptonStatus); 

  int DoNormalizedHistos();

 private:
  comptonParams* theParams;   //pointer to parameter 
  //
  //Variables output to triggeredWise root tree
  TTree* quadWiseTree;
  //accumulator data
  int nacc[8];
  double accraw[8];
  double accsig[8];
  //
  // Accumulator data histograms
  TH1F* hAcc_acc0_All;  //#triggers actually summed each MPS
  TH1F* hAcc_acc0_laserOn;  //histo of laser-on summed pulses
  TH1F* hAcc_acc0_laserOff;  //histo of laser-off summed pulses
  TH1F* hAcc_acc1_All;  //#triggers actually summed each MPS
  TH1F* hAcc_acc1_laserOn;  //histo of laser-on summed pulses
  TH1F* hAcc_acc1_laserOff;  //histo of laser-off summed pulses

  TH1F* hAcc_acc0[10];  // pulse sums for Laser Left, Positive Helicity Beam. etc
};

#endif
