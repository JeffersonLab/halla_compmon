#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

//#include "runs.h"
#include "laserUtils.h"
#include "../online/utils.h"

using namespace std;

vector<LaserPeriod> findLaserPeriods(TChain *quartetwise, Int_t runNum){
  int laserState, firstMPSnumber, beamState, dithering;
  vector<LaserPeriod> periods;

  quartetwise->SetBranchAddress("laserState", &laserState);
  quartetwise->SetBranchAddress("firstMPSnumber", &firstMPSnumber);
  quartetwise->SetBranchAddress("beamState", &beamState);
  quartetwise->SetBranchAddress("dithering", &dithering);

  int currentLaserState = -1; int stateNum = 0;
  LaserPeriod myLaserPeriod; myLaserPeriod.reset();
  for(int i = 0; i < quartetwise->GetEntries(); i++){
    quartetwise->GetEntry(i);
    myLaserPeriod.runNum = runNum;
    if(i == 0){
      currentLaserState = laserState;
      myLaserPeriod.entryStart = i;
      myLaserPeriod.mpsStart = firstMPSnumber;
      myLaserPeriod.laserState = currentLaserState;
    }
    if(currentLaserState != laserState){
      myLaserPeriod.entryEnd = i - 1;
      myLaserPeriod.mpsEnd = firstMPSnumber - 1;
      myLaserPeriod.periodNum = stateNum++;
      periods.push_back(myLaserPeriod);
      currentLaserState = laserState;
      myLaserPeriod.reset();
      myLaserPeriod.mpsStart = firstMPSnumber;
      myLaserPeriod.entryStart = i;
      myLaserPeriod.laserState = currentLaserState;
    }
    if(i == quartetwise->GetEntries() - 1){
      myLaserPeriod.entryEnd = i;
      myLaserPeriod.mpsEnd = firstMPSnumber + 7;
      myLaserPeriod.periodNum = stateNum;
      periods.push_back(myLaserPeriod);
    }
    myLaserPeriod.incBeam(beamState, dithering);
  }
  printf("  Found %i distinct laser periods\n", ((int)periods.size()));
  return periods;
}

/**
  This method cuts out any brief laser periods, and joins
  laser periods on other sides

  EXAMPLE: if the laser lock flashes on for less than three
  seconds and then flashes off, the resultant laser cycle
  will not have sufficient statistics to make a measurement
  of asymmetry. So we combine the flash with the two
  surrounding laser off states to make one big laser off
  period. (The laser on average will be cut by the
  analysis later.)
**/
vector<LaserPeriod> trimLaserPeriods(vector<LaserPeriod> periods){
  vector<LaserPeriod> trimPeriods, trimPeriods2;
  
  for(LaserPeriod lp : periods){
    if((trimPeriods.size() == 0 && (lp.laserState == 0 || lp.laserState == 1 || lp.laserState == 4)) ||
      lp.duration() < 3.0 || lp.laserState == 4)
      continue;
    trimPeriods.push_back(lp);
  }
  
  int periodNum = 0;
  //Loop through vector count periods where adjacent periods have the same laserState
  vector<vector<int>> coalescePeriods;
  int sameConsecPeriods = 0;
  int currentLaserState = -1;
  vector<int> currentPeriods;
  for(int i = 0; i < trimPeriods.size(); i++){
    if(trimPeriods[i].laserState != currentLaserState){
      if(i == 0){
        currentLaserState = trimPeriods[i].laserState;
      }
      else{
        currentLaserState = trimPeriods[i].laserState;
        coalescePeriods.push_back(currentPeriods);
        currentPeriods.clear();
      }
    }
    currentPeriods.push_back(i);
  }
  if(currentPeriods.size() > 0){
    coalescePeriods.push_back(currentPeriods);
  }

  int periodCount = 0;
  for(vector<int> periodList : coalescePeriods){
    LaserPeriod combinePeriod;
    combinePeriod.mpsStart = trimPeriods[periodList.front()].mpsStart;
    combinePeriod.entryStart = trimPeriods[periodList.front()].entryStart;
    combinePeriod.mpsEnd = trimPeriods[periodList.back()].mpsEnd;
    combinePeriod.entryEnd = trimPeriods[periodList.back()].entryEnd;
    combinePeriod.laserState = trimPeriods[periodList.front()].laserState;
    combinePeriod.periodNum = periodCount++;
    combinePeriod.runNum = trimPeriods[periodList.front()].runNum;
    
    int beamOn = 0, beamOff = 0, beamUnk = 0;
    for(int index : periodList){
      beamOn  += trimPeriods[index].beamOnEntry;
      beamOff += trimPeriods[index].beamOffEntry;
      beamUnk += trimPeriods[index].beamUnkEntry;
    }
    combinePeriod.beamOnEntry = beamOn;
    combinePeriod.beamOffEntry = beamOff;
    combinePeriod.beamUnkEntry = beamUnk;

    if(combinePeriod.laserState == 4){
      printf("  Unknown laser state lasting %.4f seconds\n", combinePeriod.duration());
    }

    trimPeriods2.push_back(combinePeriod);
  }
  printf("  Trimmed it down to %i valid periods\n", ((int)trimPeriods2.size()));
  return trimPeriods2;
}

vector<LaserCycle> findLaserCycles(vector<LaserPeriod> trimPeriods, TChain *quartetwise, vector<vector<int>> eventCuts, Int_t runNum){
  vector<LaserCycle> cycles;
  ofstream errFile; errFile.open(Form("%s/errorCodes/Run%i_errCodes.txt", getenv("COMPMON_LASERCYCLES"), runNum));
  int cycleNum = 0;
  if(trimPeriods.size() < 3) return cycles;

  for(int i = 1; i < trimPeriods.size() - 1; i++){
    if( trimPeriods[i].isOn() && trimPeriods[i - 1].isOff() && trimPeriods[i + 1].isOff() ){
      printf("  Found a candidate cycle\n");
      LaserCycle cycle;
      cycle.firstOff = trimPeriods[i - 1];
      cycle.on = trimPeriods[i];
      cycle.lastOff = trimPeriods[i + 1];
      //vector<Int_t> cycParams = balanceLaserOffs(quartetwise, cycle);
      //cycle.entryStart = cycParams[0]; cycle.mpsStart = cycParams[1]; cycle.entryEnd = cycParams[2]; cycle.mpsEnd = cycParams[3];
      //cycle.firstOffBeamOnEntry = cycParams[4]; cycle.firstOffBeamOffEntry = cycParams[5]; cycle.firstOffBeamUnkEntry = cycParams[6];
      //cycle.lastOffBeamOnEntry = cycParams[7]; cycle.lastOffBeamOffEntry = cycParams[8]; cycle.lastOffBeamUnkEntry = cycParams[9];
      cycle.entryStart = trimPeriods[i - 1].entryStart; cycle.mpsStart = trimPeriods[i - 1].mpsStart;
      cycle.entryEnd = trimPeriods[i + 1].entryEnd; cycle.mpsEnd = trimPeriods[i + 1].mpsEnd; 
      cycle.firstOffBeamOnEntry = trimPeriods[i - 1].beamOnEntry; cycle.firstOffBeamOffEntry = trimPeriods[i - 1].beamOffEntry; cycle.firstOffBeamUnkEntry = trimPeriods[i - 1].beamUnkEntry;
      cycle.lastOffBeamOnEntry = trimPeriods[i + 1].beamOnEntry; cycle.lastOffBeamOffEntry = trimPeriods[i + 1].beamOffEntry; cycle.lastOffBeamUnkEntry = trimPeriods[i + 1].beamUnkEntry;
      //cycle.firstOffBeamOnEntry  = quartetwise->GetEntries(Form("Entry$>=%i && Entry$<=%i && beamState==1", cycle.startEntry(), cycle.firstOff.entryEnd));
      //cycle.firstOffBeamOffEntry = quartetwise->GetEntries(Form("Entry$>=%i && Entry$<=%i && beamState==0", cycle.startEntry(), cycle.firstOff.entryEnd));
      //cycle.firstOffBeamUnkEntry = quartetwise->GetEntries(Form("Entry$>=%i && Entry$<=%i && beamState==2", cycle.startEntry(), cycle.firstOff.entryEnd));
      //cycle.lastOffBeamOnEntry  = quartetwise->GetEntries(Form("Entry$>=%i && Entry$<=%i && beamState==1", cycle.lastOff.entryStart, cycle.endEntry()));
      //cycle.lastOffBeamOffEntry = quartetwise->GetEntries(Form("Entry$>=%i && Entry$<=%i && beamState==0", cycle.lastOff.entryStart, cycle.endEntry()));
      //cycle.lastOffBeamUnkEntry = quartetwise->GetEntries(Form("Entry$>=%i && Entry$<=%i && beamState==2", cycle.lastOff.entryStart, cycle.endEntry()));
      //TH1F *hFirstOff = new TH1F("hFirstOff", "", 500, -50, 50);
      //TH1F *hOn = new TH1F("hOn", "", 500, -50, 50);
      //TH1F *hLastOff = new TH1F("hLastOff", "", 500, -50, 50);
      //mpswise->Project("hFirstOff", "Acc0/NAcc0", Form("mpsCoda>=%i && mpsCoda<=%i", cycle.startMPS(), cycle.firstOffEndMPS()));
      //mpswise->Project("hOn", "Acc0/NAcc0", Form("mpsCoda>=%i && mpsCoda<=%i", cycle.onStartMPS(), cycle.onEndMPS()));
      //mpswise->Project("hLastOff", "Acc0/NAcc0", Form("mpsCoda>=%i && mpsCoda<=%i", cycle.lastOffStartMPS(), cycle.endMPS()));
      //cycle.adjPedFirstOff = hFirstOff->GetMean(); cycle.adjPedOn = hOn->GetMean(); cycle.adjPedLastOff = hLastOff->GetMean();
      //cycle.pedRMSFirstOff = hFirstOff->GetRMS(); cycle.pedRMSOn = hOn->GetRMS(); cycle.pedRMSLastOff = hLastOff->GetRMS();
      //delete hFirstOff; delete hOn; delete hLastOff;
      Int_t valid = cycle.isValidCycle(eventCuts);
      if(valid == 0){
        printf("    Found valid laser cycle %i\n", cycleNum + 1);
        cycle.cycleNum = cycleNum++;
        cycles.push_back(cycle);
      }
      else if((valid & 0x20) != 0){printf("    Cycle manually cut by events\n");}
      else{printf("    Cycle was invalid with error code %04X...\n", valid);}
      errFile<<Form("%i\n", valid);
    }
  }
  errFile.close();
  return cycles;
}

vector<LaserCycle> runCycles(int runNum){
  vector<TChain *> runChains = loadChain(runNum);
  Int_t quartetwise = 1;
  //TFile *fin = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  //TTree *quartetwise = (TTree *)fin->Get("quartetwise");
  //TTree *mpswise = (TTree *)fin->Get("mpswise");
  vector<vector<int>> eventCuts = readEventCuts(runNum);
  return findLaserCycles(trimLaserPeriods(findLaserPeriods(runChains[quartetwise], runNum)), runChains[quartetwise], eventCuts, runNum);
}

void printCycles(vector<LaserCycle> cycles, int runNum){
  for(LaserCycle cycle : cycles){
    printf("Run %i, cycle %i: PreOff: %i-%i, On: %i-%i, PostOff: %i-%i\n", runNum, cycle.cycleNum,
            cycle.mpsStart, cycle.firstOff.mpsEnd, cycle.on.mpsStart, cycle.on.mpsEnd, cycle.lastOff.mpsStart, cycle.mpsEnd);
    printf("  First off: On: %i, Off: %i, Unk: %i\n", cycle.firstOff.beamOnEntry, cycle.firstOff.beamOffEntry, cycle.firstOff.beamUnkEntry);
    printf("  On: On: %i, Off: %i, Unk: %i\n", cycle.on.beamOnEntry, cycle.on.beamOffEntry, cycle.on.beamUnkEntry);
    printf("  Last off: On: %i, Off: %i, Unk: %i\n", cycle.lastOff.beamOnEntry, cycle.lastOff.beamOffEntry, cycle.lastOff.beamUnkEntry);
  }
}

void writeCycles(vector<LaserCycle> cycles, int runNum){
  ofstream output; output.open(Form("%s/cycles_%i.dat", getenv("COMPMON_CYCLES"), runNum));
  for(LaserCycle cycle : cycles){
    output<<cycle.mpsStart<<","<<cycle.firstOffEndMPS()<<","<<cycle.onStartMPS()<<","<<cycle.onEndMPS()<<","<<cycle.lastOffStartMPS()<<","<<cycle.mpsEnd<<"\n";
  }
  output.close();
  printf("Created cycles file %s/cycles_%i.dat\n", getenv("COMPMON_CYCLES"), runNum);
}


void laserCycles(int runNum){
  printf("Finding laser cycles for run %i...\n", runNum);
  vector<LaserCycle> allCycles = runCycles(runNum);
  writeCycles(allCycles, runNum);
}

/**
void laserCycles(Int_t runNum){
  TFile *fin = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *quartetwise = (TTree *)fin->Get("quartetwise");
  vector<LaserPeriod> allPeriods = trimLaserPeriods(findLaserPeriods(quartetwise, runNum));
  for(LaserPeriod p : allPeriods){
    printf("Period %i, Laser %i: MPS Limits: %i-%i, Entry Counts: (%i, %i, %i)\n", p.periodNum, p.laserState, p.mpsStart, p.mpsEnd,
           p.beamOnEntry, p.beamOffEntry, p.beamUnkEntry);
  }
}
**/
/**
void laserCyclesAJZ(int runNum){
  TFile *fin = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *quartetwise = (TTree *)fin->Get("quartetwise");

  vector<LaserPeriod> periods = findLaserPeriods(quartetwise);
  vector<LaserPeriod> periods2 = trimLaserPeriods(periods);
  vector<LaserCycle> cycles = findLaserCycles(periods2);
  for(LaserPeriod period : periods){
    printf("Period number %2i has state %s from MPS %6i to MPS %6i with duration %3.2f and with %.2f%% beam on\n", 
            period.periodNum, period.laserStatus().Data(), period.mpsStart, period.mpsEnd, period.duration(), period.fracBeamOn()*100);
  }

  printf("\n*****************************\n\n");

  for(LaserPeriod period : periods2){
    printf("Period2 number %2i has state %s from MPS %6i to MPS %6i with duration %3.2f and with %.2f%% beam on\n", 
            period.periodNum, period.laserStatus().Data(), period.mpsStart, period.mpsEnd, period.duration(), period.fracBeamOn()*100);
  }

  printf("\n*****************************\n\n");

  for(LaserCycle cycle : cycles){
    printf("Created cycle %2i from %6i to %6i lasting %3.2f seconds.\n",
           cycle.cycleNum, cycle.startMPS(), cycle.endMPS(), cycle.duration());
  }

  printf("\n*****************************\n\n");

  vector<vector<int>> runList = productionRunList(1);
  for(int snail = 0; snail < runList.size(); snail++){
    printf("Snail %i:\n", snail + 1);
    for(int run = 0; run < runList[snail].size(); run++){
      printf("    Run %i\n", runList[snail][run]);
    }
  }
  
}
**/
