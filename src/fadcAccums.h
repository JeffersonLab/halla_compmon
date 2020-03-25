/*
//	fadcAccums
//      Analysis of fadc accumulator dada
//      and various MPS-wise and Multipletwise (previously just Quarterwise) data
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

#ifndef fadcAccums_h
#define fadcAccums_h

#define MAX_DAC_VALUES 20
#define NUM_ACCUM_TYPES 8
#define NUM_HISTORY_LASER_PERIODS 10
#define HELPLUS 1
#define HELMINUS 0

// laserstate now defined in class comptonStatus
//enum laserstate {RIGHT, LEFT, RIGHTOFF, LEFTOFF, UNKNOWN};

class fadcAccums {
 public:
  fadcAccums(textParams* theParams, comptonStatus* theStatus);
  void newRun();   //init counters at start of a run
  int DefineHistos();
  int DefineTree();
  TTree *GetTree() { return multipletWiseTree; } // Useful for when ROOT splits files, to get back the current file
  void labelSpinSortedHistos(TH1F* histo);
  int DoAccums(vmeauxdata* theVMEauxdata,fadcdata *theFADCdata); 
  int BuildMultiplet(); 
  int DoLaserTransition(int wrapup);    //called when laser state                                           //transition encountered and at end of run
  int DoNormalizedHistos();
  bool GetSummedAccumulatorRecord(int record,int accumNumber, int* laserState, int*  beamState,
				  int* countMults, double* numerator, double* denominator);
  void StoreAccumulatorInfo(int laserState, int beamState, int countMults,
			   double numerator[NUM_ACCUM_TYPES],
			    double denominator[NUM_ACCUM_TYPES]);

 private:
  char multipletName[200]; //name of multiplet pattern (i.e. Quartet, Octet)
  int helicityStructure; // 2 (pair), 4 (quartet), 8 (octet)
  textParams* theParams;   //pointer to instance of comptonParameter class
  comptonStatus* theStatus;  //pointer to instance of comptonStatus class
  // parameters read by textParams class
  float accumHistoRange;
  float accumHistoDiffMax;
  float accumHistoScaledRange;
  float channel_calorimeter_PMT;
  float ped_value;

  // Define the standard compton variable tree
  comptonHelTree *comptonTree;
  //
  //Variables output to triggeredWise root tree
  TTree* mpsWiseTree;
  //accumulator data
  int nacc[NUM_ACCUM_TYPES];
  int64_t accraw[NUM_ACCUM_TYPES];  //raw acumuilator data
  comptonVariable<double> FIXHELaccsig[NUM_ACCUM_TYPES];   //pedestal corrected accumulator data
  comptonVariable<double> FIXHELnaccsig[NUM_ACCUM_TYPES];   //number of samples in accumulator data
  double mpsPedestal; // Pedestal as determined from triggered+randoms sums.
  double mpsRandomPedestal; // Pedestal as determined from randoms sums.
  double mpsTriggerPedestal; // Pedestal as determined from triggered sums.
  //
  // Multipletwise and Laserwise info
  TTree* multipletWiseTree;
  int lastSubMult;   //track where we are within a helicity multiple set (pair,mult,oct, etc...)
  bool multipletValid;  //used by valid multiplet checking algorithm
  bool multipletStable; //true if laser and beam states stable over multiplet
  int prevLaserStateEncountered;
  int laserStateBeingSummed;
  int beamStateBeingSummed;
  double accSumLastLaserOff[NUM_ACCUM_TYPES];//keep for background subtraction
  double bcmLastLaserOff;  //keep for background subtractio

  //sum periods within a Quartet
  int laserStateThisMult;     //laserState for this Mult 
  int beamStateThisMult;        //beam on, off, or unknown for this Mult
  double accPosMult[NUM_ACCUM_TYPES];  //Sum of Positive helicity bins
  double accNegMult[NUM_ACCUM_TYPES];  //Sum of Negative helicity
  double naccPosMult[NUM_ACCUM_TYPES];  //# samples
  double naccNegMult[NUM_ACCUM_TYPES];  //# samples

  double beamPosMult;  //same for summed beam charge
  double beamNegMult;
  int countPosMinNeg;  //count positive HW windows minus negative HW windos
  //
  //Sum up multiplets within each laser period 
  int countMultsLaserPeriod;
  int firstMPS;
  int multipletHelicity;   //lower bits get helcity pattern
  int multipletReportedHelicity;   // the current lower bits get helcity pattern (should be same as multipletHelicity if not delayed)
  double accPosLaserPeriod[NUM_ACCUM_TYPES]; //sum within a laser period
  double accNegLaserPeriod[NUM_ACCUM_TYPES];
  double beamPosLaserPeriod;
  double beamNegLaserPeriod;
  double epics_multipletBCM;
  //
  //data base for history of last N laser period sums
  int pointerHistory;   //points to most recent mult info in histor arrays
  int countMultsHistory[NUM_HISTORY_LASER_PERIODS];
  int laserStateHistory [NUM_HISTORY_LASER_PERIODS];
  int beamStateHistory[NUM_HISTORY_LASER_PERIODS];
  double accDiffHistory[NUM_ACCUM_TYPES][NUM_HISTORY_LASER_PERIODS]; 
  double accSumHistory[NUM_ACCUM_TYPES][NUM_HISTORY_LASER_PERIODS];
  double beamDiffHistory[NUM_HISTORY_LASER_PERIODS]; 
  double beamSumHistory[NUM_HISTORY_LASER_PERIODS];
  //
  //Multiplet Asymmetries
  //accumulator 0
  
  TH1F* hQ_aligned0_beamOff;  //Acc0 spin aligned beam off
  TH1F* hQ_anti0_beamOff;  //Acc0 spin anti-aligned beam off
  TH1F* hQ_sum0_beamOff;
  TH1F* hQ_diff0_beamOff;
  TH1F* hQ_Asym0_raw_beamOff;

  TH1F* hQ_aligned0_laserOff;
  TH1F* hQ_anti0_laserOff;
  TH1F* hQ_aligned0_laserOn;
  TH1F* hQ_anti0_laserOn;
  TH1F* hQ_sum0_laserLeft;
  TH1F* hQ_sum0_laserOff;
  TH1F* hQ_diff0_laserLeft;
  TH1F* hQ_diff0_laserOff;


  TH1F* hQ_Asym0_raw_laserLeft;
  TH1F* hQ_Asym0_raw_laserOff;
  TH1F* hQ_Signal0;
  TH1F* hQ_Asym0;

  //accumulator 4  (stretched window)
  TH1F* hQ_sum4_beamOff;
  TH1F* hQ_diff4_beamOff;
  TH1F* hQ_sum4_laserLeft;
  TH1F* hQ_sum4_laserOff;
  TH1F* hQ_diff4_laserLeft;
  TH1F* hQ_diff4_laserOff;
  TH1F* hQ_Asym4_raw_laserLeft;
  TH1F* hQ_Asym4_raw_laserOff;
  TH1F* hQ_Signal4;
  TH1F* hQ_Asym4;

//beam
  TH1F* hQ_BCM_laserOff_P;
  TH1F* hQ_BCM_laserOff_N;
  TH1F* hQ_BCM_laserOn_P;
  TH1F* hQ_BCM_laserOn_N;
  TH1F* hQ_BCM_Asym;



  //laser-wise histos
  TH1F* hL_sum0_laserLeft;
  TH1F* hL_diff0_laserLeft;
  TH1F* hL_sum0_laserOff;
  TH1F* hL_diff0_laserOff;
  TH1F* hL_Asym0_raw_laserLeft;
  TH1F* hL_Asym0_raw_laserOff;

  //
  // for a complete (laser-off, laser_on, laser-off) cycle)
  double accAsymRaw[NUM_ACCUM_TYPES];
  double accDilution[NUM_ACCUM_TYPES];
  double accAsym[NUM_ACCUM_TYPES];  
  //
  // MPS-wise histograms
  //
  TH1F* hM_BCM;                 //beam charge
  TH1F* hM_BCM_BeamOn;         //beam charge 
  TH1F* hM_BCM_BeamOff;         //beam charge 
  TH1F* hM_CavPower;
  TH1F* hM_CavPower_LaserOn;
  TH1F* hM_CavPower_LaserOff;
  TH1F* hM_BPMSum;
  TH1F* hM_BPMSum_BeamOn;
  TH1F* hM_BPMSum_BeamOff;

  TH1F* hM_SpinState;            //# MPS events vs Spin State
  TH1F* hM_SpinState_BeamOn;    //Spin State-sorted MPS events for Beam On
  TH1F* hM_BCM_SummedBySpinState;  //BCM sorted by Spin State
  TH1F* hM_Trigs_Accepted;
  TH1F* hM_Trigs_Prescaled;
  TH1F* hM_Trigs_Scaler;
  TH1F* hM_numSums;
  //
  TH1F* hM_acc0_everything;    //all beam&laser states, wide range
  TH1F* hM_acc_All[NUM_ACCUM_TYPES];  //#triggers actually summed each MPS
  TH1F* hM_acc_laserOn[NUM_ACCUM_TYPES];  //histo of laser-on summed pulses
  TH1F* hM_acc_laserOff[NUM_ACCUM_TYPES];  //histo of laser-off summed pulses
  //just histo acc0 for beamOff
  TH1F* hM_acc0_beamOff_LaserOn;
  TH1F* hM_acc0_beamOff_laserOff;
  TH1F* hM_acc4_beamOff_LaserOn;
  TH1F* hM_acc4_beamOff_laserOff;
  //
  // Misc Histogram
  TH2F* hS_Asym0_raw;   //strip charge of mult raw asymmetry
};

#endif
