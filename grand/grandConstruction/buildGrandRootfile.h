#ifndef buildGrandRootfile_h
#define buildGrandRootfile_h

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "vars.h"

vector<vector<Float_t>> keys;
Int_t numCyclesAcc;
// cyc tree
// basics
Int_t snailNum, runNum, cycleNum;
Int_t firstOffStartMPS, firstOffEndMPS, onStartMPS, onEndMPS, lastOffStartMPS, lastOffEndMPS;
Int_t cycSign;
Float_t cycQW1, cycHW1, cycQW2, cycIHWP, cycVWien, cycHWien, cycSolWien;
PolVar cycLaserPol;
// mpswise
//vector<vector<TString>> cycMPSVars;
vector<DataVar> cycMPSData;
// quartetwise
//vector<vector<TString>> cycQrtVars;
vector<DataVar> cycQrtData;
//triggerwise
//Float_t meanPedF, meanErrPedF, rmsPedF, rmsErrPedF;
//Float_t meanPedL, meanErrPedL, rmsPedL, rmsErrPedL;
DataVar firstOffPedestal, lastOffPedestal;
// calculations
Float_t cycleTime;
vector<TString> cycQrtPols;
vector<PolVar> cycQrtCalc;
Float_t meanDiff0LasOn, meanErrDiff0LasOn, meanSum0LasOn, meanErrSum0LasOn;
Float_t meanDiff0LasOff, meanErrDiff0LasOff, meanDiff0LasOff1, meanErrDiff0LasOff1, meanDiff0LasOff2, meanErrDiff0LasOff2;
Float_t meanSum0LasOff, meanErrSum0LasOff, meanSum0LasOff1, meanErrSum0LasOff1, meanSum0LasOff2, meanErrSum0LasOff2;
Float_t meanDiff4LasOn, meanErrDiff4LasOn, meanSum4LasOn, meanErrSum4LasOn;
Float_t meanDiff4LasOff, meanErrDiff4LasOff, meanDiff4LasOff1, meanErrDiff4LasOff1, meanDiff4LasOff2, meanErrDiff4LasOff2;
Float_t meanSum4LasOff, meanErrSum4LasOff, meanSum4LasOff1, meanErrSum4LasOff1, meanSum4LasOff2, meanErrSum4LasOff2;
PolVar anPow;
PolVar collOffset;
Float_t alphaOn, alphaOff;
PolVar sigSubSum0, sigSubSum4;
PolVar asym0LasOnAlt, asym0LasOffAlt, asym0LasOff1Alt, asym0LasOff2Alt;
PolVar asym4LasOnAlt, asym4LasOffAlt, asym4LasOff1Alt, asym4LasOff2Alt;
PolVar asym0, asym0NGC, asym0LasOn, asym0LasOff, asym0LasOff1, asym0LasOff2, pol0;
PolVar asym4, asym4NGC, asym4LasOn, asym4LasOff, asym4LasOff1, asym4LasOff2, pol4;
Int_t cycleCut;
Float_t acc0OnMode, acc0OnOffset, acc0OffMode, acc0OffOffset;

// run tree
// basics
Int_t numRunCycles, numRunCyclesAcc;
// mpswise
//vector<vector<TString>> runMPSVars;
vector<DataVar> runMPSData;
// epicswise
//vector<vector<TString>> runEpcVars;
vector<StdVar> runEpcData;
// calculations
Float_t runTime;
Int_t runSign;
vector<Float_t> runPol0Avg, runPol0Err;
vector<Float_t> runAsym0Avg, runAsym0Err;
vector<Float_t> runAsym0NGCAvg, runAsym0NGCErr;
vector<Float_t> runAsym0OnAvg, runAsym0OnErr;
vector<Float_t> runAsym0OffAvg, runAsym0OffErr;
FitPolVar runAsym0, runAsym0NGC, runAsym0On, runAsym0Off, runPol0;
Float_t runQW1, runHW1, runQW2, runIHWP, runVWienAngle, runHWienAngle, runPhiFG;
PolVar runLaserPol;

// snl tree
// basics
Int_t numSnlCycles, numSnlCyclesAcc;
Int_t numSnlRuns;
//calculations
Float_t snailTime;
Int_t snailSign;
//vector<Float_t> snlPol0Avg, snlPol0Err, snlPol4Avg, snlPol4Err;
//vector<Float_t> snlAsym0Avg, snlAsym0Err, snlAsym4Avg, snlAsym4Err;
//FitPolVar snlAsym0, snlAsym4, snlPol0, snlPol4;
vector<Float_t> snlPol0Avg, snlPol0Err;
vector<Float_t> snlAsym0Avg, snlAsym0Err;
vector<Float_t> snlAsym0NGCAvg, snlAsym0NGCErr;
vector<Float_t> snlAsym0OnAvg, snlAsym0OnErr;
vector<Float_t> snlAsym0OffAvg, snlAsym0OffErr;
FitPolVar snlAsym0, snlAsym0NGC, snlAsym0On, snlAsym0Off, snlPol0;
vector<Float_t> qw1Lst, hw1Lst, qw2Lst, ihwpLst, VWienAngleLst, HWienAngleLst, PhiFGLst;
Float_t qw1, hw1, qw2, ihwp, VWienAngle, HWienAngle, PhiFG;
//PolVar laserPol;
Float_t laserPol;

Int_t helicityFreq(Int_t runNum){
  if(runNum > 4242 && runNum < 4622)
    return 240;
  else
    return 120;
}

/**
  This function grabs the BPM information that was calculated in cycIterSet() in
  buildGrandRootfile.C and puts it in a 2D vector for use in calculating
  analyzing power. The outside vector contains two vectors each with 4 elements 
  each. The first vector contains the BPM mean positions for the cycle being
  considered, while the second contains their statistical error for the cycle.

  The four elements represent the four positions tracked by 2A and 2B.

  They go in this order:
  Index 0 = bpmAx, Index 1 = bpmAy, Index 2 = bpmBx, Index 3 = bpmBy

  So the mean position for bpmBx (for example) in this vector would be
  cycBPMs[0][2]

  And the mean error on bpmAy (for example) would be
  cycBPMs[1][1]
**/
vector<vector<Float_t>> getCyclePositions(Int_t runNumber){
  vector<vector<Float_t>> cycBPMs;
  //vector<Float_t> bpmMeans, bpmErrs;
  ifstream infile(Form("%s/Run%i_bpms.csv", getenv("COMPMON_RUNPLOTS"), runNumber));
  string readStr;
  while(getline(infile, readStr)){
    stringstream ss(readStr);
    string token;
    vector<Float_t> bpmLine;
    while(getline(ss, token, ',')){
      bpmLine.push_back(atof(token.c_str()));
    }
    cycBPMs.push_back(bpmLine);
  }

  //Float_t axErr = TMath::Sqrt(TMath::Power(cycMPSData[10].meanErr, 2));
  //Float_t ayErr = TMath::Sqrt(TMath::Power(cycMPSData[11].meanErr, 2));
  //Float_t bxErr = TMath::Sqrt(TMath::Power(cycMPSData[12].meanErr, 2));
  //Float_t byErr = TMath::Sqrt(TMath::Power(cycMPSData[13].meanErr, 2));

  //bpmMeans.push_back(cycMPSData[10].mean); bpmErrs.push_back(axErr);
  //bpmMeans.push_back(cycMPSData[11].mean); bpmErrs.push_back(ayErr);
  //bpmMeans.push_back(cycMPSData[12].mean); bpmErrs.push_back(bxErr);
  //bpmMeans.push_back(cycMPSData[13].mean); bpmErrs.push_back(byErr);

  //cycBPMs.push_back(bpmMeans); cycBPMs.push_back(bpmErrs);

  return cycBPMs;
}

void calcCollOffset(Int_t sNum, Int_t rNum){
  Float_t collCentX = 0.0; Float_t collCentXErr = 0.0;
  Float_t collCentY = 0.0; Float_t collCentYErr = 0.0;
  if(sNum < 100 || sNum == 500){
    collCentX = 2.160; collCentXErr = 0.080;
    collCentY = 0.709; collCentYErr = 0.098;
  }
  vector<vector<Float_t>> bpms = getCyclePositions(rNum);
  Float_t collPosX = bpms[2][0] + 6*(bpms[0][0] - bpms[2][0]); 
  Float_t collPosY = bpms[3][0] + 6*(bpms[3][0] - bpms[1][0]);

  Float_t xDiffErr = 6*TMath::Sqrt(TMath::Power(bpms[0][1], 2) + TMath::Power(bpms[2][1], 2));
  Float_t yDiffErr = 6*TMath::Sqrt(TMath::Power(bpms[1][1], 2) + TMath::Power(bpms[3][1], 2));
  Float_t xProjErr = TMath::Sqrt(TMath::Power(bpms[2][1], 2) + TMath::Power(xDiffErr, 2));
  Float_t yProjErr = TMath::Sqrt(TMath::Power(bpms[3][1], 2) + TMath::Power(yDiffErr, 2));

  Float_t collOffX = collPosX - collCentX;
  Float_t collOffY = collPosY - collCentY;
  Float_t collOffXErr = TMath::Sqrt(TMath::Power(xProjErr, 2) + TMath::Power(collCentXErr, 2));
  Float_t collOffYErr = TMath::Sqrt(TMath::Power(yProjErr, 2) + TMath::Power(collCentYErr, 2));

  collOffset.mean = TMath::Sqrt(TMath::Power(collOffX, 2) + TMath::Power(collOffY, 2));
  Float_t offXsqErr = TMath::Abs(TMath::Power(collOffX, 2))*2*collOffXErr/collOffX;
  Float_t offYsqErr = TMath::Abs(TMath::Power(collOffY, 2))*2*collOffYErr/collOffY;
  Float_t offSqSumErr = TMath::Sqrt(TMath::Power(offXsqErr, 2) + TMath::Power(offYsqErr, 2));
  collOffset.meanErr = TMath::Abs(collOffset.mean)*0.5*offSqSumErr/(TMath::Power(collOffX, 2) + TMath::Power(collOffY, 2));
}

void prexAnPow(Int_t sNum){
  Float_t p0 = 0.0165612; Float_t p0Err = 9.98921E-4;
  Float_t p1 = -0.00001020; Float_t p1Err = 0.00001082;
  Float_t p2 =  0.00001591; Float_t p2Err = 0.00000267;

  Float_t collOffSq = collOffset.mean*collOffset.mean;
  Float_t collOffSqErr = TMath::Abs(collOffSq)*2*collOffset.meanErr/collOffset.mean;
  Float_t t1 = p1*collOffset.mean;
  Float_t t1Err = TMath::Abs(t1)*TMath::Sqrt(TMath::Power(p1Err/p1, 2) + TMath::Power(collOffset.meanErr/collOffset.mean, 2));
  Float_t t2 = p2*collOffSq;
  Float_t t2Err = TMath::Abs(t2)*TMath::Sqrt(TMath::Power(p2Err/p2, 2) + TMath::Power(collOffSqErr/collOffSq, 2));
  anPow.mean = p0 + t1 + t2;
  //anPow.meanErr = TMath::Sqrt(TMath::Power(p0Err, 2) + TMath::Power(t1Err, 2) + TMath::Power(t2Err, 2));
  //anPow.mean = 0.01655915;
  anPow.meanErr = 0.0;
}

/**
  The variable anPow is defined above in this file and is a PolVar.
  The definition of PolVar can be found in vars.h in this folder. In words,
  a PolVar is a struct with two Float_t's: mean and meanErr. anPow gets
  written to the cyc tree in the grand rootfile, so it is defined independently
  for each cycle.
 
  This function is called once per cycle: at the start of calcCyclePol() in
  buildGrandRootfile.C. anPow is a variable shared between these two files so
  any change in the .C file affects the .h file too.

  Right now the function sets the mean analyzing powers from simulation for
  PREX and CREX plus some legacy DVCS code I'm not sure I can get rid of.
  
  meanErr's are set to zero since there is no error from the simulations
  (or at least it is vanishingly small. I've left them in because the bpm's
  will have some uncertainty attached to them from statistics that will
  probably need to get propagated upwards. My guess is that these will too
  be vanishingly small, but it's best to be thorough.

  The BPM position function is called but never used. This is where we implement
  the position correction. All we need is a function to propagate the BPM positions
  onto the collimator surface and determine the analyzing power change from the
  offset from center.
**/
void setAnalyzingPower(Int_t sNum, Int_t run_num){
  vector<vector<Float_t>> bpmMeans = getCyclePositions(run_num);
  if(run_num < 3800){                                               //DVCS runs (added for backwards compatibility)
    anPow.mean = 0.1105;
    anPow.meanErr = 0.0;
   }                                
  else if(run_num >= 4232 && run_num <= 4929){                      //PREX runs
    //anPow.mean = 0.01655915;
    //anPow.meanErr = 0.0;
    prexAnPow(sNum);
  }
  else if(run_num >= 4930){                                         //CREX runs
    anPow.mean = 0.036017934;
    anPow.meanErr = 0.0;
  }
  else{
    printf("Run doesn't have a defined analyzing power.\n");
    anPow.mean = 0.0;
    anPow.meanErr = 0.0;
  }
}

/**
vector<vector<TString>> readVarsFile(TString outTree, TString inTree){
  ifstream infile(Form("%s/grandConfigs/%s_%s_vars.csv", getenv("COMPMON_GRAND"), outTree.Data(), inTree.Data()));
  string readStr;
  vector<vector<TString>> vars;
  while(getline(infile, readStr)){
    vector<TString> oneVar; stringstream ss(readStr);
    string token;
    while(getline(ss, token, ',')){
      oneVar.push_back(token.c_str());
    }
    vars.push_back(oneVar);
  }
  return vars;
}

vector<TString> readPolsFile(TString outTree, TString inTree){
  ifstream infile(Form("%s/grandConfigs/%s_%s_pol.csv", getenv("COMPMON_GRAND"), outTree.Data(), inTree.Data()));
  string readStr;
  vector<TString> vars;
  while(getline(infile, readStr)){
    //printf("Read from file: %s\n", readStr.c_str());
    vars.push_back(readStr.c_str());
  }
  return vars;
}
**/

vector<vector<int>> findCycles(int runNum){
  ifstream infile(Form("%s/cycles_%i.dat", getenv("COMPMON_CYCLES"), runNum));
  string readStr;
  int cyclesInThisRun = 0;
  vector<vector<int>> run;
  while(getline(infile, readStr)){
    vector<int> limits; stringstream ss(readStr);
    for(int i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    cyclesInThisRun++;
    vector<int> cycle;
    for(int i = 0; i < limits.size(); i++)
      cycle.push_back(limits[i]);
    run.push_back(cycle);
  }
  printf("Looked in %s/cycles_%i.dat\n", getenv("COMPMON_CYCLES"), runNum);
  printf("Found %i cycles\n", (Int_t)run.size());
  return run;
}

void resetChain(TChain *ch){
  if(ch){
    delete ch;
  }
}

vector<TChain *> loadChain(Int_t runnum){
  TChain *mpswise = 0; TChain *quartetwise = 0;
  TChain *pulserwise = 0; TChain *triggerwise = 0;
  TChain *epicswise = 0; TChain *runwise = 0;
  TChain *snapshots = 0;

  resetChain(mpswise); resetChain(quartetwise);
  resetChain(pulserwise); resetChain(triggerwise);
  resetChain(epicswise); resetChain(runwise);
  resetChain(snapshots);

  mpswise = new TChain("mpswise");
  quartetwise = new TChain("quartetwise");
  pulserwise = new TChain("pulserwise");
  triggerwise = new TChain("triggerwise");
  epicswise = new TChain("epicswise");
  runwise = new TChain("runwise");
  snapshots = new TChain("snapshots");
  std::vector<TChain*> chains;
  chains.push_back(mpswise);
  chains.push_back(quartetwise);
  chains.push_back(pulserwise);
  chains.push_back(triggerwise);
  chains.push_back(epicswise);
  chains.push_back(runwise);
  chains.push_back(snapshots);

  TString filesPre = Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),runnum);
  int nfiles = 0;
  for(size_t ch = 0; ch < chains.size(); ch++) {
    nfiles = chains[ch]->Add(filesPre+".root");
    nfiles += chains[ch]->Add(filesPre+"_*.root");

    if(nfiles<=0) {
	    std::cerr << "Looked for files under: " << filesPre+".root" << std::endl;
	    std::cerr << "Found no files to plot!" << std::endl;
	    return chains;
    }
  }
  return chains;
}

#endif
