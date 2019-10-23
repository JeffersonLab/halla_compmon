#ifndef GLOBAL_H
#define GLOBAL_H

#include "MyMath.h"
#include "laser.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <TChain.h>
#include <TCanvas.h>
#include <TString.h>
#include <TGraphErrors.h>

//const char ROOTFILE_PATH="../rootfiles"
//const char *ROOTFILE_PATH="/home/cornejo/scratch/compton/rootfiles";
//const char *ROOTFILE_PATH="/data/cmuwork/rootfiles/tests_2018";
//const char *ROOTFILE_PATH="/data/cmuwork/rootfiles/DVCS_Tests";
const char *ROOTFILE_PATH="$COMP_ROOTFILES";
TString gComptonOutPath;
TString gComptonROOTFILES;

// Some constants (can be changed by source files)
int kMinEntriesForPattern = 25*0+100;
TChain *gChain = 0;
int gStartEntry = 0;
int gEntries = 0;
int gRun;
bool kUseJC2Rootfiles=false;

const double kStdError = 1e6;

TChain *gChainSnapshots = 0;
int gEntriesSnapshots = 0;

bool kSplitLaserOffRegion = true;
bool kPrintPatternHistory = false;
const int kPrintLaserCycleDetails = true;
bool kPrintIntermediateCycleInfo = true;
bool kPrintBCMResults = false;
bool kSaveGraphs = true;
bool kAllowBeamOffCycles = true;
const int kMaxAccCount = 100e3;
const int kMinCycleCounts = 10*30;
const int kMaxCycleSeparation = 45*30;
const int kHelStructure = 4;
const bool kBeamOn = true;
const bool kUseRawBCM = true;
const double kMinCurrent = 90.; // When requiring beam ON only
const double kMaxCurrent = 1.25; // When requiring beam OFF only
const int kGraphColors[6] = { kRed, kGreen+1, kRed+2, kRed+1, kBlue+1, kPink+1 };
const int kGraphStyles[7] = { 20, 21, 24, 25, 22, 23};
const int kMaxDS1Bg = 250;
const double kMaxAccDS1BgFrac = 25e3*1e3;
const double kAvgMPSTime = 0.03383; // seconds
const double kAvgMPSTimeDeviation = 0.00002;  // seconds
const bool kBCMNormalizeMPSLevel = false;
const int kMinCycleBeamOff = 30;
const int kMinCycleBeamOn = 300;
const int kMinCycleMpsCount = 30*240/8.; // 30 seconds

const char *kHelTreeName[2] = {"Neg", "Pos" };
const int kNAccums = 6;


// Tree variables
double Acc[kNAccums][2];
double NAcc[kNAccums][2];
double Acc0;
double Acc4;
float posCavPowerCalibrated;
float negCavPowerCalibrated;
float posBCM;
float negBCM;
int beamState;
int laserState;
float rawBCM;
int mpsCoda;
double NAcc0;
int NAcc4;
int helicityState;
int actualHelicity;
int reportedHelicity;
int patternCounter;
int runScaler0;
int runScaler1;
int runScaler4;
// Snapshots tree variables
int numSamples;
float snap[1000];
int snapClock;
int mpsCodaSnapshots;
float bcmSnapshots;
int beamStateSnapshots;
int laserStateSnapshots;
int randomTime;

MyData gAnalyzingPower[kNAccums];

// Define some common notation
#define PATTERN_MINUS 0
#define PATTERN_PLUS 1
#define PATTERN_UNKNOWN -1
#define LASER_RIGHTON = 0
#define LASER_LEFTON = 1
#define LASER_RIGHTOFF = 2
#define LASER_LEFTOFF = 3
#define HELICITY_INFO_ACTUAL 0
#define HELICITY_INFO_REPORTED 1

#define PRINT_NAME_PADDING 8

#define NUM_STATES 5 // There are 4 states (OFF0, ON, OFF1, BKG, Sub)
const char *gStateNames[5] = {"OFF0","ON  ","OFF1","Bkg ","Sub "};
const char *gStateSafeNames[5] = {"OFF0","ON","OFF1","Bkg","Sub"};
const char *gStateSafeNamesPretty[5] = {"LEFT-OFF","ON","RIGHT-OFF","AVG BKG","SUBTRACTED"};
TString gLaserCyclesFileName;


enum BeamCheckEnum {
 BEAM_ON,
 BEAM_OFF,
 BEAM_NO_CHECK
};

enum LaserCheckEnum {
  LASER_ON_LEFT,
  LASER_OFF_LEFT,
  LASER_ON_RIGHT,
  LASER_OFF_RIGHT,
  LASER_NO_CHECK
};

enum SpinStatesLaser {
  // Laser ON
  SPIN_LEFTMINUS_ON,
  SPIN_LEFTPLUS_ON,
  SPIN_RIGHTMINUS_ON,
  SPIN_RIGHTPLUS_ON,
  // Laser OFFF
  SPIN_LEFTMINUS_OFF,
  SPIN_LEFTPLUS_OFF,
  SPIN_RIGHTMINUS_OFF,
  SPIN_RIGHTPLUS_OFF
};



/*
typedef struct {
  int leftOffStart;
  int leftOffEnd;
  int onStart;
  int onEnd;
  int rightOffStart;
  int rightOffEnd;
  int status;
} LaserPattern_t;
*/

std::vector<int> gLaserCyclesStatus;
std::vector<int> gLaserCyclesStart[3];
std::vector<int> gLaserCyclesEnd[3];
std::vector<std::vector<std::vector<double> > > gPedestalResults;
std::vector<double> gPedestals[4];


void findRange(std::vector<double> &vals, std::vector<double> &errs,
    double &min, double &max, bool reset = false)
{
  for(size_t i = 0; i < vals.size(); i++) {
    if(vals[i]-errs[i]<min) {
      min = vals[i]-errs[i];
    }
    if(vals[i]+errs[i]>max) {
      max = vals[i]+errs[i];
    }
  }
}

void paddRange(double &min, double &max, double padd = 0.1)
{
  double range = (max-min);
  min -= range*padd;
  max += range*padd;
}

void setGraphStyle(TGraphErrors * graph, int marker, int color, int width = 1)
{
  graph->SetMarkerStyle(marker);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineWidth(width);
  graph->SetFillColor(0);
}

const char * GetCycleFileName()
{
  if(gLaserCyclesFileName.IsNull())
    gLaserCyclesFileName.Form("%s/laserCycles_%d.dat",gComptonOutPath.Data(),gRun);
  return gLaserCyclesFileName.Data();
}

void loadChain(TChain *chain)
{
  if(kUseJC2Rootfiles)
    //chain->Add(TString::Format("../rootfiles/jc2_%d.root",gRun));
    chain->Add(TString::Format("%s/jc2_%d.root",ROOTFILE_PATH,gRun));
  else
    //chain->Add(TString::Format("../rootfiles/compmon_%d.root",gRun));
    chain->Add(TString::Format("%s/compmon_%d.root",gComptonROOTFILES.Data(),gRun));
}


bool readCycles()
{
  fstream in;
  int c;
  int s_off0;
  int e_off0;
  int s_on;
  int e_on;
  int s_off1;
  int e_off1;
  int status;
  in.open(gComptonOutPath+TString::Format("/laserCycles_%d.dat",
        gRun),std::ios::in);
  if(!in.is_open())
    return false;

  while( in >> c &&
      in >> s_off0 && in >> e_off0 &&
      in >> s_on && in >> e_on &&
      in >> s_off1 && in >> e_off1 &&
      in >> status && !in.eof() ) {
    gLaserCyclesStart[0].push_back(s_off0);
    gLaserCyclesEnd[0].push_back(e_off0);
    gLaserCyclesStart[1].push_back(s_on);
    gLaserCyclesEnd[1].push_back(e_on);
    gLaserCyclesStart[2].push_back(s_off1);
    gLaserCyclesEnd[2].push_back(e_off1);
    gLaserCyclesStatus.push_back(status);
  }

  if(gLaserCyclesStatus.size() > 0) {
    std::cout << "Succesfully read " << gLaserCyclesStatus.size() << 
      " laser cycles from file." << std::endl;
    return true;
  }

  return false;
}

bool checkCurrent() {
  float bcm = (posBCM + negBCM)/2.0;
  if(kUseRawBCM) {
    if(kBeamOn && rawBCM >= 1350) {
      return true;
    }

    if(!kBeamOn && rawBCM <= 505)
      return true;
  } else {
    if(kBeamOn && bcm >= kMinCurrent)
      return true;

    if(!kBeamOn && bcm <= kMaxCurrent)
      return true;
  }

  return false;
}


void storeCyclesInFile()
{
  fstream outFile;
  outFile.open(GetCycleFileName(),std::ios::out);
  std::cout << "Saving cycles to file: " << GetCycleFileName() << std::endl;
  for(size_t c = 0; c < gLaserCycles.size(); c++) {
    outFile << std::setw(3) << std::setfill('0') << c << " "
      << std::setfill(' ')
      << std::setw(1) << gLaserCycles[c].laserOn << " "
      << std::setw(7) << gLaserCycles[c].start << " "
      << std::setw(7) << gLaserCycles[c].end << " "
      << std::setw(7) << gLaserCycles[c].mpsCount << " "
      << std::setw(5) << gLaserCycles[c].beamOnCount << " "
      << std::setw(5) << gLaserCycles[c].beamOffCount << " "
      << std::setw(1) << gLaserCycles[c].beamOnStatus << " "
      << std::setw(1) << gLaserCycles[c].beamOffStatus << " "
      << std::setw(1) << gLaserCycles[c].cycleStatus << std::endl;
  }
  outFile.close();
}

void printCycle(int cycleNum, LaserCycle_t cycle)
{
  std::cout << "Cycle[" << std::setw(3) << std::setfill('0')
    << cycleNum << "] ";
  if(cycle.laserOn) {
    std::cout << "Laser ON  ";
  } else {
    std::cout << "Laser OFF ";
  }
  std::cout << std::setfill(' ');
  std::cout << "mpsRange[" << std::setw(7) << cycle.start << ", "
    << std::setw(7) << cycle.end << "] ";
  std::cout << "Good MPS Count: " << std::setw(7) << cycle.mpsCount << " ";
  std::cout << " Beam Count ON: " << std::setw(5)
    << cycle.beamOnCount << " OFF: " << std::setw(5)
    << cycle.beamOffCount << std::endl;
}

void printCyclePattern(int num, LaserPattern_t pat)
{
  std::cout << "Pattern[" << std::setw(3) << std::setfill('0') << num << "] "
    << std::setfill(' ')
    << std::setw(3) << pat.offLeftNum
    << "[" << std::setw(7) << pat.offLeft->start << ", "
    << std::setw(7) << pat.offLeft->end << "] "
    << std::setw(3) << pat.onNum
    << "[" << std::setw(7) << pat.on->start << ", "
    << std::setw(7) << pat.on->end << "] "
    << std::setw(3) << pat.offRightNum
    << "[" << std::setw(7) << pat.offRight->start << ", "
    << std::setw(7) << pat.offRight->end << "] " << std::endl;

}

bool findLaserPatterns()
{
  if(gLaserCycles.size() < 2) {
    return false;
  }

  // Now read in all the patterns
  LaserCycle_t *cyc[3];
  bool goodBeamOn;
  bool goodBeamOff;
  bool goodCycle;
  for(size_t c = 0; c < gLaserCycles.size()-2; c++) {
    goodBeamOn = goodBeamOff = goodCycle =  true;
    for(int k = 0; k < 3; k++) {
      cyc[k] = &gLaserCycles[c+k];
      if(!cyc[k]->beamOnStatus)
        goodBeamOn = false;
      if(!cyc[k]->beamOffStatus)
        goodBeamOff = false;
      if(!cyc[k]->cycleStatus) {
        goodCycle = false;
        std::cout << "Bad cycle status (skipping): " << c+k << std::endl;
      }
    }
    if( goodCycle &&
        !cyc[0]->laserOn && cyc[1]->laserOn && !cyc[2]->laserOn ) {
      if( cyc[1]->start - cyc[0]->end < kMaxCycleSeparation &&
          cyc[2]->start - cyc[1]->end < kMaxCycleSeparation) {
          LaserPattern_t pat;
          pat.offLeft = cyc[0];
          pat.on = cyc[1];
          pat.offRight = cyc[2];
          pat.offLeftNum = c;
          pat.onNum = c+1;
          pat.offRightNum = c+2;
          /*
          pat.leftOffStart = cyc[0].start;
          pat.leftOffEnd = cyc[0].end;
          pat.onStart = cyc[1].start;
          pat.onEnd = cyc[1].end;
          pat.rightOffStart = cyc[2].start;
          pat.rightOffEnd = cyc[2].end;
          */
          pat.status = goodCycle;
        if(goodBeamOn) {
          pat.patNum = gLaserPatterns.size();
          gLaserPatterns.push_back(pat);
        }
        if(goodBeamOff) {
          pat.patNum = gLaserPatternsBeamOff.size();
          gLaserPatternsBeamOff.push_back(pat);
        }
      }
    }
  }

  if(kPrintLaserCycleDetails) {
    for(size_t p = 0; p < gLaserPatterns.size(); p ++) {
      std::cout << "BEAM ON  ";
      printCyclePattern(p,gLaserPatterns[p]);
    }
    for(size_t p = 0; p < gLaserPatternsBeamOff.size(); p ++) {
      std::cout << "BEAM OFF ";
      printCyclePattern(p,gLaserPatternsBeamOff[p]);
    }
  }

  if(gLaserPatterns.size() > 0 || gLaserPatternsBeamOff.size() > 0)
    return true;

  return false;

}

bool ReadCyclesFromFile()
{
  fstream in;
  int num;
  LaserCycle_t c;
  in.open(GetCycleFileName(),std::ios::in);
  if(!in.is_open())
    return false;

  while( in >> num &&
      in >> c.laserOn &&
      in >> c.start &&
      in >> c.end &&
      in >> c.mpsCount &&
      in >> c.beamOnCount &&
      in >> c.beamOffCount &&
      in >> c.beamOnStatus &&
      in >> c.beamOffStatus &&
      in >> c.cycleStatus &&
      !in.eof() ) {
    gLaserCycles.push_back(c);
  }

  if(gLaserCycles.size() > 0) {
    std::cout << "Succesfully read " << gLaserCycles.size() << 
      " laser cycles from file. Now looking for patterns." << std::endl;
    return findLaserPatterns();
  }

  return false;
}


bool findLaserCycles()
{
  LaserCycle_t tmp;

  // Get the first entry to start it off
  int startEntry = 0;
  for(startEntry = gStartEntry; startEntry < gEntries; startEntry++) {
    gChain->GetEntry(gStartEntry);
    if(laserState < 4  && beamState < 3)
      break;
  }
  tmp.start = mpsCoda;
  tmp.end = mpsCoda;
  tmp.start_entry = startEntry;
  tmp.end_entry = startEntry;
  tmp.mpsCount = 0;
  if(laserState < 2) { // laserON
    tmp.laserOn = true;
  } else if (laserState < 4) { // laser OFF
    tmp.laserOn = false;
  }
  tmp.beamOnCount = 0;
  tmp.beamOffCount = 0;
  if(beamState==1) {
    tmp.beamOnCount = 1;
  } else if (beamState == 0) {
    tmp.beamOffCount = 1;
  }
  tmp.beamOnStatus = 1;
  tmp.beamOffStatus = 1;
  tmp.cycleStatus = 1;
  std::cout << "Start entry: " << startEntry << ". Entries: " << gEntries << std::endl;

  bool laserOn;
  // Now loop through the rest of the mps's to find the laser cycles
  for(int entry = startEntry; entry < gEntries; entry++) {
    gChain->GetEntry(entry);

    // Skip bad laser states (they should be handled by the analysis
    // scripts themselves)
    if(laserState >= 4) {
      continue;
    }

    if(laserState < 2) {
      laserOn = true;
    } else {
      laserOn = false;
    }

    // Check if this laser State is the same as the previous one
    if(laserOn == tmp.laserOn) {
      if(beamState==0) {
        tmp.beamOffCount++;
      } else if (beamState==1) {
        tmp.beamOnCount++;
      }
      tmp.end = mpsCoda;
      tmp.end_entry = entry;
      tmp.mpsCount++;
    } else { // new state, start new cycle
      if(tmp.mpsCount < kMinCycleMpsCount ) {
        tmp.cycleStatus = 0;
      }
      if(tmp.beamOnCount < kMinCycleBeamOn) {
        tmp.beamOnStatus = 0;
      }
      if(tmp.beamOffCount < kMinCycleBeamOff) {
        tmp.beamOffStatus = 0;
      }
      // Print the cycle info
      if(kPrintLaserCycleDetails) {
        printCycle(gLaserCycles.size(),tmp);
      }

      // Add cycle to vector
      gLaserCycles.push_back(tmp);

      // Reset temporary variable
      tmp.start = mpsCoda;
      tmp.end = mpsCoda;
      tmp.laserOn = laserOn;
      tmp.mpsCount = 1;
      tmp.beamOnCount = 0;
      tmp.beamOffCount = 0;
      tmp.mpsCount = 1;
      tmp.cycleStatus = 1;
      tmp.beamOnStatus = 1;
      tmp.beamOffStatus = 1;
      tmp.start_entry = entry;
      tmp.end_entry = entry;
      if(beamState==0) {
        tmp.beamOffCount++;
      } else if (beamState==1) {
        tmp.beamOnCount++;
      }
    }
  }
  // Now store the cycles to file
  storeCyclesInFile();

  return findLaserPatterns();
}


bool findLaserCycles2()
{
  std::vector<int> cycles[4];
  int lastFound = 4;
  int startEntry;

  // Find the first valid laser state (i.e., not 4)
  for(startEntry = gStartEntry; startEntry < gEntries &&
      lastFound == 4; startEntry++) {
    gChain->GetEntry(startEntry);
    lastFound = laserState;
  }

  int start = startEntry;
  int end = startEntry;
  int counts = 0;
  //int bcmCount = 0;
  int badFound = 0;
  bool endFound = false;

  // Now find all start and end of laser states
  for(int entry = startEntry; entry < gEntries; entry++) {

    gChain->GetEntry(entry);

    if(laserState == 4) {
      badFound++;
      continue;
    }

    if(lastFound != laserState || (badFound>30) ) {
      if(kPrintIntermediateCycleInfo) {
        std::cout << "Found state( " << cycles[0].size() << "): " << laserState
          << " [" << start << ", " << end << "] with " << counts << " counts."
          << std::endl;
      }
      //if(bcmCount >= counts/2) {
        cycles[0].push_back(start);
        cycles[1].push_back(end);
        cycles[2].push_back(counts);
        cycles[3].push_back(lastFound);
      //}
      start = entry;
      end = entry;
      // change to track mpsCoda
      start = mpsCoda;
      end = mpsCoda;
      counts = 1;
      //bcmCount = 1;
    } else {
      end = entry;
      end = mpsCoda;
      counts++;
      //if(checkCurrent())
        //bcmCount++;
      //if((kBeamOn&&bcm >= kMinCurrent)||(!kBeamOn&&bcm<= kMaxCurrent))
      //  bcmCount++;
    }

    lastFound = laserState;
    badFound = 0;
  }

  if(cycles[0].size() < 3)
    return false;

  size_t startCycle = 0;
  bool fewCounts;
  for(size_t i = 0; i < cycles[0].size()-2; i++) {
    fewCounts = false;
    for(size_t j = 0; j < 3; j++) {
      if(cycles[2][i+j] < kMinCycleCounts)
        fewCounts = true;
    }
    if( !fewCounts &&
        cycles[3][i] >=2 &&
        cycles[3][i] == cycles[3][i+2] &&
        cycles[3][i+1] < 2 ) {
      if( (cycles[0][i+1] - cycles[1][i]) < kMaxCycleSeparation &&
          (cycles[0][i+2] - cycles[1][i+1]) < kMaxCycleSeparation ) {
        // Found valid OFF-ON-OFF cycle!
        for(size_t k = 0; k < 3; k++) {
          gLaserCyclesStart[k].push_back(cycles[0][i+k]);
          gLaserCyclesEnd[k].push_back(cycles[1][i+k]);
        }
        gLaserCyclesStatus.push_back(1);
      } else {
        /*
        std::cout << "Cycle not good: ";
        for(int l = 0; l < 3; l++) {
          std::cout << "[" << cycles[0][l+i] << ", " << cycles[1][l+i] << "] ";
        }
        std::cout << std::endl;
        */
      }
    }
  }

  // Print out all cycles
  fstream laserOutFile;
  laserOutFile.open(gComptonOutPath+TString::Format("/laserCycles_%d.dat",
        gRun),std::ios::out);
  for(size_t m = 0; m < gLaserCyclesStart[0].size(); m++) {
    std::cout << "Laser Pattern[" << m << "]: ";
    laserOutFile << m << " ";
    for(size_t n = 0; n < 3; n++) {
      std::cout << "[" << gLaserCyclesStart[n][m] << ", "
        << gLaserCyclesEnd[n][m] << "] ";
      laserOutFile << gLaserCyclesStart[n][m] << " "
        << gLaserCyclesEnd[n][m] << " ";
    }
    laserOutFile << gLaserCyclesStatus[m];  // Default status is "good"
    laserOutFile << std::endl;
    std::cout << std::endl;
  }
  laserOutFile.close();
  std::cout << "Saved laser cycles to: " << gComptonOutPath.Data() << "/laserCycles_" << gRun << ".dat" << std::endl;
  return true;
}

class VComptonVariable {
public:
  VComptonVariable(TString name, TString units, int print,VComptonVariable *norm = 0) {
    Config(name,units,print,norm);
  };
  virtual ~VComptonVariable() {};
  const char *GetUnits() { return fUnits.Data(); }
  const char* GetName() { return fName.Data(); }
  const char* GetNameRight() { return fNamePaddedRight.Data(); }
  const char* GetNameLeft() { return fNamePaddedLeft.Data(); }
  int GetPrintType() { return fPrintType; }
  VComptonVariable* getNormVar() { return fNormVar; };
  virtual double val(int) = 0;

private:
  void Config(TString name, TString units, int print, VComptonVariable *norm) {
    SetName(name);
    fUnits = units;
    fPrintType = print;
    fNormVar = norm;
  }
  void SetName(TString name) {
    fName = name;
    fNamePaddedRight = name;
    fNamePaddedLeft = name;
    while(fNamePaddedRight.Length() < PRINT_NAME_PADDING) {
      fNamePaddedLeft += " ";
      fNamePaddedRight.Insert(0," ");
    }
  }
  TString fName;
  TString fUnits;
  TString fNamePaddedRight;
  TString fNamePaddedLeft;
  int fPrintType;
  VComptonVariable *fNormVar;
};

template<typename T>
class ComptonVariable : public VComptonVariable {
public:
  ComptonVariable(TString name, T *p, TString units,
      int print, VComptonVariable *norm = 0) :
    VComptonVariable(name,units,print,norm), fVal(p) {
  }
  virtual ~ComptonVariable() { fVal = 0; }
  virtual double val(int h) {
    if(getNormVar()) {
      return fVal[h]/getNormVar()->val(h);
    }
    return double(fVal[h]);
  }
private:
  T* fVal;
};



class ComptonVariableOld {
private:
  TString fName;
  TString fUnits;
  TString fNamePaddedRight;
  TString fNamePaddedLeft;
  int fType;
  int *fIntType;
  float *fFloatType;
  double *fDoubleType;
  bool fNormalizeToBCM;
  bool fSubtractBkg;
  bool fComputePol;
  int fPrintType;

public:
  ComptonVariableOld() : fType(-1), fIntType(0), fFloatType(0), fDoubleType(0){};
  virtual ~ComptonVariableOld(){
    fIntType = 0;
    fFloatType = 0;
    fDoubleType = 0;
  };

  ComptonVariableOld(TString name, int *p, TString units,int print, bool norm, bool sub, bool pol = false) {
    Config(name,1,units,print,norm,sub,pol);
    fIntType = p;
  }
  ComptonVariableOld(TString name, float *p, TString units,int print, bool norm, bool sub, bool pol = false) {
    Config(name,1,units,print,norm,sub,pol);
    fType = 2;
    fFloatType = p;
  }
  ComptonVariableOld(TString name, double *p, TString units,int print, bool norm, bool sub, bool pol = false) {
    Config(name,1,units,print,norm,sub,pol);
    fType = 3;
    fDoubleType = p;
  }
  void Config(TString name, int type, TString units, int print, bool norm, bool sub, bool pol) {
    SetName(name);
    fUnits = units;
    fType = type;
    fNormalizeToBCM = norm;
    fSubtractBkg = sub;
    fComputePol = pol;
    fPrintType = print;
  }
  void SetName(TString name) {
    fName = name;
    fNamePaddedRight = name;
    fNamePaddedLeft = name;
    while(fNamePaddedRight.Length() < PRINT_NAME_PADDING) {
      fNamePaddedLeft += " ";
      fNamePaddedRight.Insert(0," ");
    }
  }
  double GetValue() {
    if(fType == 1) {
      return *fIntType;
    } else if (fType == 2) {
      return *fFloatType;
    } else if (fType == 3) {
      return *fDoubleType;
    }
    return -1e6;
  }
  const char *GetUnits() { return fUnits.Data(); }
  const char* GetName() { return fName.Data(); }
  const char* GetNameRight() { return fNamePaddedRight.Data(); }
  const char* GetNameLeft() { return fNamePaddedLeft.Data(); }
  bool BCMNormalize() { return fNormalizeToBCM; }
  bool BackgroundSubtract() { return fSubtractBkg; }
  bool ComputePol() { return fComputePol; }
  int GetPrintType() { return fPrintType; }
};

void printPatternHistory(Int_t history, Int_t validEntries)
{
  std::cout << "History: ";
  bool state;
  for(int i = 31; i >= 0 ; i--) {
    state = (history>>i)&0x1;
    if(i>validEntries)
      std::cout << "X";
    else if(state)
      std::cout << "+";
    else
      std::cout << "-";
  }
  std::cout << std::endl;
}

bool findHelicityPattern(TChain *chain, int entries, int &startEntry, int &helicityState,
    bool printPattern = false)
{
  unsigned int stateHistory = 0;
  // Find the start of a pattern
  for(Int_t entry = startEntry; entry < entries; entry++) {
    chain->GetEntry(entry);
    stateHistory = (stateHistory<<1)+helicityState;
    if(printPattern)
      printPatternHistory(stateHistory,entry);
    int patCheck = stateHistory&0x1F;//&0x7;
    // (03/11/2016 jc2: modified to include more history)
    //if((patCheck == 2 || patCheck == 5) && entry >= kMinEntriesForPattern) {
    if((patCheck == 0x12 || patCheck == 0x0D) && entry >= kMinEntriesForPattern) {
      startEntry = entry;
      std::cout << "Found start of pattern at entry " << entry << std::endl;
      return true;
    }
  }

  return false;

}

int GetPatternType(int *pattern, int helStructure)
{
  if(helStructure == 4) {
    if( pattern[0] == 0 && pattern[1] == 1 && pattern[2] == 1 && pattern[3] == 0)
      return PATTERN_MINUS;

    if( pattern[0] == 1 && pattern[1] == 0 && pattern[2] == 0 && pattern[3] == 1)
      return PATTERN_PLUS;
  } else if (helStructure == 2) {
    if( pattern[0] == 0 && pattern[1] == 1)
      return PATTERN_MINUS;

    if( pattern[0] == 1 && pattern[1] == 0)
      return PATTERN_PLUS;
  }

  return PATTERN_UNKNOWN;
}

void processPattern(int type, std::vector<double> &data, std::vector<double> &result, int helStructure)
{
  double h03;
  double h12;
  double sign = type==PATTERN_MINUS ? -1. : 1.;
  double sum;
  double diff;
  double asym;

  if(helStructure == 4) {
    // To make comparisions with the pair-wise case, average
    // the H0 and H1 values
    h03 = (data[0] + data[3])/2.;
    h12 = (data[1] + data[2])/2.;
  } else if ( helStructure == 2) {
    h03 = data[0];
    h12 = data[1];
  }
  sum = h03+h12;
  diff = sign*(h03-h12);
  asym = 100.*diff/sum;

  if(type == 0) {
    result[0] = h03;
    result[1] = h12;
  } else if (type == 1) {
    result[1] = h03;
    result[0] = h12;
  }
  result[2] = sum;
  result[3] = diff;
  result[4] = asym;
  result[5] = sign;
}

std::string PatternToString(int *pat, int structure)
{
  std::string result;
  for(int i = 0; i < structure; i++) {
    result += (pat[i] == 1 ? "+" : "-");
  }
  return result;
}

void printLineF( std::fstream &out, int length = 80)
{
  for(int i = 0; i < length; i++) {
    out << "-";
  }
  out << std::endl;
}

void printLine(int length = 80)
{
  for(int i = 0; i < length; i++) {
    std::cout << "-";
  }
  std::cout << std::endl;
}

double getAnalyzingPower(int run)
{
  // 4GeV running
  if(run >= 2516 && run <= 2576)
    //return 5.00392;
    return 6.7431;
  //if(run>= 2577 )
  if(run < 4000 )
    return 11.4644;

  return 0.0166;

  //std::cerr << "For run " << run << " bad analyzing power. " << std::endl;
  return -1e6;
}

void saveHisto(TCanvas *canvas, TString name)
{
  TString dir = "results/";
  canvas->SaveAs(dir+name);
}

void printHelper(std::ostringstream &myout, std::fstream *out = 0)
{
  if(out) {
    *out << myout.str();
  } else {
    std::cout << myout.str();
  }
}

void fmtNumber(double num, std::ostringstream &myout)
{
  // The idea is to make them all so that the decimal is always at the same
  // place for all lines.
  double mult = 1e7;
  int width = 17;
  while(num < mult && mult > 1) {
    myout << " ";
    width -=1;
    mult /= 10.;
  }
  if(num>=0) {
    myout << " ";
    width -=1;
  }
  myout << std::left << std::setprecision(6) << std::fixed << std::setw(width) << num;
}

void printNumber(double num, std::fstream *out = 0)
{
  std::ostringstream myout;
  // The idea is to make them all so that the decimal is always at the same
  // place for all lines.
  double mult = 1e7;
  int width = 17;
  while(num < mult && mult > 1) {
    myout << " ";
    width -=1;
    mult /= 10.;
  }
  if(num>=0) {
    myout << " ";
    width -=1;
  }
  myout << std::left << std::setprecision(6) << std::fixed << std::setw(width) << num;
  printHelper(myout,out);
  //printf("% 12.4f",num);
  //int width = 17;
  //if(num>=0) {
  //  std::cout << " ";
  //  width -=1;
  //}
  //std::cout << std::left << std::setprecision(4) << std::fixed << std::setw(width) << num;
}

void printResult(TString string, double v, double e, TString unit, std::fstream *out = 0)
{
  std::ostringstream myout;
  myout << string.Data() << " : ";
  fmtNumber(v,myout);
  myout << " +/- ";
  fmtNumber(e,myout);
  if(!unit.IsNull())
    myout << " (" << unit.Data() << ")";
  myout << std::endl;
  printHelper(myout,out);
}

void printResult(TString string, MyHistoData dat, TString unit, std::fstream *out = 0)
{
  printResult(string,dat.val(),dat.err(),unit,out);
/*
  std::cout << string.Data() << " : ";
  printNumber(dat.val());
  std::cout << " +/- ";
  printNumber(dat.err());
  if(!unit.IsNull())
    std::cout << " (" << unit.Data() << ")";
  std::cout << std::endl;
*/
}

void printResult(TString string, MyData dat, TString unit, std::fstream *out = 0)
{
  printResult(string,dat.val(),dat.err(),unit,out);
/*
  std::cout << string.Data() << " : ";
  printNumber(dat.val());
  std::cout << " +/- ";
  printNumber(dat.err());
  if(!unit.IsNull())
    std::cout << " (" << unit.Data() << ")";
  std::cout << std::endl;
*/
}

void printResult(TString string, MyErrorWeightedData dat, TString unit, std::fstream *out = 0)
{
  printResult(string,dat.val(),dat.err(),unit,out);
/*
  std::cout << string.Data() << " : ";
  printNumber(dat.val());
  std::cout << " +/- ";
  printNumber(dat.err());
  if(!unit.IsNull())
    std::cout << " (" << unit.Data() << ")";
  std::cout << std::endl;
*/
}



void readAnalyzingPower()
{
  fstream in;
  in.open("analyzing_power.dat",std::ios::in);
  int acc,run;
  double val,err;
  for(acc = 0; acc < 6; acc++) {
    gAnalyzingPower[acc].Set(-1.0,-1.0);
  }

  std::string line;
  int ncomment = 0;
  while(std::getline(in >> std::ws,line)) {
    // Find any comments and get rid of them
    ncomment = line.find("#");
    if(ncomment != std::string::npos )
      line.erase(ncomment);
    ncomment = line.find(";");
    if(ncomment != std::string::npos )
      line.erase(ncomment);
    // If what follows has a run, acc and ath with its error, then continue
    std::istringstream iss(line);
    if( iss >> run && iss >> acc && iss >> val && iss >> err ) {
      // Is our run greater or equal to this run?
      if(run<=gRun && acc >= 0 && acc<kNAccums) {
        gAnalyzingPower[acc].Set(val,err);
      }
    }
  }

  std::cout << "Setting run " << gRun << " analyzing powers to: " << std::endl;
  for(acc = 0; acc < kNAccums; acc++ ) {
    if(gAnalyzingPower[acc].val() != -1 ) {
      printf("Acc%d : %8.5f +/- %8.5f\n",acc,gAnalyzingPower[acc].val(),
          gAnalyzingPower[acc].err());
    }
  }

  in.close();
}




bool init(int run, bool readCyclesFromFile)
{
  gRun = run;

  // Format the output and input paths
  gComptonOutPath=TString::Format("%s/runs/Run%d",getenv("COMPMON_WEB"),gRun);
  gComptonROOTFILES=TString::Format("%s",getenv("COMP_ROOTFILES"));

  // Clear gChain
  if(gChain)
    delete gChain;

  gChain = new TChain("quartetwise");
  //gChain->Add(TString::Format("../rootfiles/FadcCalo2016_%d.root",run));
  loadChain(gChain);

  // Disable all branches in the gChain
  gChain->SetBranchStatus("*",0);

  // Enable standard variables
  gChain->SetBranchStatus("helicityState",1);
  gChain->SetBranchStatus("PosHelBCM",1);
  gChain->SetBranchStatus("NegHelBCM",1);
  gChain->SetBranchStatus("beamState",1);
  gChain->SetBranchStatus("laserState",1);
  gChain->SetBranchStatus("rawBCM",1);
  gChain->SetBranchStatus("PosHelCavPowerCalibrated",1);
  gChain->SetBranchStatus("NegHelCavPowerCalibrated",1);
  gChain->SetBranchStatus("mpsCoda",1);
  for(Int_t h = 0; h < 2; h++) {
    for(Int_t acc = 0; acc < kNAccums; acc++) {
      gChain->SetBranchStatus(
          TString::Format("%sHelAcc%d",kHelTreeName[h],acc),1);
      gChain->SetBranchStatus(
          TString::Format("%sHelNSamples%d",kHelTreeName[h],acc),1);
    }
  }

  // Now setup the branch pointers
  gChain->SetBranchAddress("helicityState",&helicityState);
  gChain->SetBranchAddress("beamState",&beamState);
  gChain->SetBranchAddress("laserState",&laserState);
  gChain->SetBranchAddress("rawBCM",&rawBCM);
  gChain->SetBranchAddress("PosHelCavPowerCalibrated",&posCavPowerCalibrated);
  gChain->SetBranchAddress("NegHelCavPowerCalibrated",&negCavPowerCalibrated);
  gChain->SetBranchAddress("PosHelBCM",&posBCM);
  gChain->SetBranchAddress("NegHelBCM",&negBCM);
  gChain->SetBranchAddress("mpsCoda",&mpsCoda);
  for(Int_t h = 0; h < 2; h++) {
    for(Int_t acc = 0; acc < kNAccums; acc++) {
      gChain->SetBranchAddress(
          TString::Format("%sHelAcc%d",kHelTreeName[h],acc),&Acc[acc][h]);
      gChain->SetBranchAddress(
          TString::Format("%sHelNSamples%d",kHelTreeName[h],acc),&NAcc[acc][h]);
    }
  }

  // Get the total entries from the tree
  gEntries =  gChain->GetEntries();

  /*
   // First, find the first complete pattern in the tree
  if(!findHelicityPattern(gChain,gEntries,gStartEntry,helicityState,kPrintPatternHistory)) {
    std::cerr << "Helicity pattern not found! Exiting." << std::endl;
    return false;
  }
  */

  // Now find all the laser cycles
  bool readCycleStatus = false;
  if(readCyclesFromFile) {
    readCycleStatus = ReadCyclesFromFile();
    if(!readCycleStatus)
      std::cerr << "Could not read cycle information from data. "
        << " Will parse the tree and look for cycles again." << std::endl;
  }

  if(!readCycleStatus && !findLaserCycles()) {
    std::cerr << "No valid cycles found in tree. Exiting!" << std::endl;
    return false;
  }

  // Setup the snapshot tree
  if(gChainSnapshots)
    delete gChainSnapshots;

  gChainSnapshots = new TChain("snapshots");
  //gChainSnapshots->Add(TString::Format("../rootfiles/FadcCalo2016_%d.root",run));
  loadChain(gChainSnapshots);

  // Disable all branches
  gChainSnapshots->SetBranchStatus("*",0);

  // Enable only the two branches we want
  gChainSnapshots->SetBranchStatus("numSamples",1);
  gChainSnapshots->SetBranchStatus("snap",1);
  gChainSnapshots->SetBranchStatus("snapClock",1);
  gChainSnapshots->SetBranchStatus("mpsCoda",1);
  gChainSnapshots->SetBranchStatus("bcm",1);
  gChainSnapshots->SetBranchStatus("beamState",1);
  gChainSnapshots->SetBranchStatus("laserState",1);
  gChainSnapshots->SetBranchStatus("randomTime",1);

  gChainSnapshots->SetBranchAddress("numSamples",&numSamples);
  gChainSnapshots->SetBranchAddress("snap",&snap);
  gChainSnapshots->SetBranchAddress("snapClock",&snapClock);
  gChainSnapshots->SetBranchAddress("mpsCoda",&mpsCodaSnapshots);
  gChainSnapshots->SetBranchAddress("bcm",&bcmSnapshots);
  gChainSnapshots->SetBranchAddress("beamState",&beamStateSnapshots);
  gChainSnapshots->SetBranchAddress("laserState",&laserStateSnapshots);
  gChainSnapshots->SetBranchAddress("randomTime",&randomTime);

  gEntriesSnapshots = gChainSnapshots->GetEntries();


  // Read in the analyzing power for this run
  readAnalyzingPower();

  return true;

}

void readPedestals()
{
  fstream in;
  in.open(TString::Format("results/pedestals_%d.dat",gRun),
      std::ios::in);
  int c;
  double off0[2];
  double on[2];
  double off1[2];
  double bk[2];
  while( in >> c
      && in >> off0[0] && in >> off0[1]
      && in >> on[0] && in >> on[1]
      && in >> off1[0] && in >> off1[1]
      && in >> bk[0] && in >> bk[1]
      && !in.eof()) {
    gPedestals[0].push_back(off0[0]);
    gPedestals[0].push_back(off0[1]);
    gPedestals[1].push_back(on[0]);
    gPedestals[1].push_back(on[1]);
    gPedestals[2].push_back(off1[0]);
    gPedestals[2].push_back(off1[1]);
    gPedestals[3].push_back(bk[0]);
    gPedestals[3].push_back(bk[1]);
  }
}


#endif // GLOBAL_H
