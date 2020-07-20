#ifndef LASERUTILS_H
#define LASERUTILS_H

#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

typedef struct{
  int mpsStart;
  int mpsEnd;
  int entryStart;
  int entryEnd;
  int laserState;
  int periodNum;
  int beamOnEntry;
  int beamOffEntry;
  int beamUnkEntry;
  int runNum;

  Float_t helicityFreq(){
    if(runNum > 4242 && runNum < 4621){return 240.0;}
    else{return 120.0;}
  }
  bool mpsIsInPeriod(int mpsCoda){return (mpsCoda >= mpsStart && mpsCoda <= mpsEnd);}
  bool entryIsInState(int entry) {return (entry >= entryStart && entry <= entryEnd);}
  float duration(){return (mpsEnd - mpsStart)*1.0/helicityFreq();}
  float fracBeamOn(){return beamOnEntry*1.0/(beamOnEntry + beamOffEntry + beamUnkEntry);}
  bool isOn(){return laserState == 0 || laserState == 1;}
  bool isOff(){return laserState == 2 || laserState == 3;}

  void reset(){
    mpsStart = -1; mpsEnd = -1; entryStart = -1; entryEnd = -1; laserState = -1;
    periodNum = -1; beamOnEntry = 0; beamOffEntry = 0; beamUnkEntry = 0; runNum = -1;
  }

  void incBeam(int beamState){
    if(beamState == 0) beamOffEntry++;
    else if(beamState == 1) beamOnEntry++;
    else beamUnkEntry++;
  }

  TString laserStatus(){
    if(laserState == 0 || laserState == 1) return "ON ";
    else if(laserState == 2 || laserState == 3) return "OFF";
    else return "UNK";
  }
} LaserPeriod;

typedef struct{
  LaserPeriod firstOff;
  LaserPeriod on;
  LaserPeriod lastOff;
  int cycleNum;
  //Float_t adjPedFirstOff, adjPedOn, adjPedLastOff;
  //Float_t pedRMSFirstOff, pedRMSOn, pedRMSLastOff;
  Int_t firstOffBeamOnEntry, firstOffBeamOffEntry, firstOffBeamUnkEntry;
  Int_t lastOffBeamOnEntry, lastOffBeamOffEntry, lastOffBeamUnkEntry;
  Int_t entryStart, entryEnd, mpsStart, mpsEnd;

  //int startMPS(){return (firstOff.mpsStart + firstOff.mpsEnd)/2;}
  //int endMPS()  {return (lastOff.mpsStart + lastOff.mpsEnd)/2 - 1;}
  //int startEntry(){return (firstOff.entryStart + firstOff.entryEnd)/2;}
  //int endEntry()  {return (lastOff.entryStart + lastOff.entryEnd)/2 - 1;}
  //int startMPS(){return firstOff.mpsStart;}
  //int endMPS()  {return lastOff.mpsEnd;}
  //int startEntry(){return firstOff.entryStart;}
  //int endEntry()  {return lastOff.entryEnd;}
  int firstOffEndMPS(){return firstOff.mpsEnd;}
  int onStartMPS(){return on.mpsStart;}
  int onEndMPS(){return on.mpsEnd;}
  int lastOffStartMPS(){return lastOff.mpsStart;}
  float duration(){return (mpsEnd - mpsStart)*1.0/on.helicityFreq();}
  float firstOffFracBeamOn(){return firstOffBeamOnEntry*1.0/(firstOffBeamOnEntry + firstOffBeamOffEntry + firstOffBeamUnkEntry);}
  float onFracBeamOn(){return on.fracBeamOn();}
  float lastOffFracBeamOn(){return lastOffBeamOnEntry*1.0/(lastOffBeamOnEntry + lastOffBeamOffEntry + lastOffBeamUnkEntry);}
  bool evaluateEventCuts(vector<vector<int>> eventCuts){
    bool noEventCuts = true;
    for(vector<int> cut : eventCuts){
      noEventCuts = noEventCuts && !((cut[0] <= mpsStart && cut[1] >= mpsStart) || (cut[0] <= mpsEnd && cut[1] >= mpsEnd) || 
                                     (cut[0] <= mpsStart && cut[1] >= mpsEnd) || (cut[0] >= mpsStart && cut[1] <= mpsEnd));
    }
    return noEventCuts;
  }
  Int_t isValidCycle(vector<vector<int>> eventCuts){
    bool correctPattern = firstOff.isOff() && on.isOn() && lastOff.isOff();
    bool laserOnBeamOn = on.duration()*on.fracBeamOn() > 3.0;
    //bool laserOffBeamOn = firstOffFracBeamOn()*100 > 10.0 && lastOffFracBeamOn()*100 > 10.0;
    bool laserOffBeamOn = firstOffBeamOnEntry*1.0/30.0 > 3.0 && lastOffBeamOnEntry > 3.0;
    bool periodSeparation12 = (on.mpsStart - firstOff.mpsEnd)*1.0/on.helicityFreq() < 10.0;
    bool periodSeparation23 = (lastOff.mpsStart - on.mpsEnd)*1.0/on.helicityFreq() < 10.0;
    //bool narrowPedestal = pedRMSFirstOff < 0.1 && pedRMSOn < 0.1 && pedRMSLastOff < 0.1;
    bool noEventCuts = evaluateEventCuts(eventCuts);
    Int_t cutFlag = 0x1*(Int_t)(!correctPattern) + 0x2*(Int_t)(!laserOnBeamOn) + 0x4*(Int_t)(!laserOffBeamOn);
    cutFlag += 0x8*(Int_t)(!periodSeparation12) + 0x10*(Int_t)(!periodSeparation23) + 0x20*(Int_t)(!noEventCuts);
    //printf("Cut Summary: \n");
    //printf("  Correct Pattern: %i\n  Laser On, Beam On: %i\n", correctPattern, laserOnBeamOn);
    //printf("  Laser Off, Beam On: %i\n  1-2 Separation: %i\n", laserOffBeamOn, periodSeparation12);
    //printf("  2-3 Separation: %i\n  No Event Cuts: %i\n", periodSeparation23, noEventCuts);
    //printf("Period Separations: \n");
    //printf("    MPS limits: %i-%i, %i-%i, %i-%i\n", mpsStart, firstOff.mpsEnd, on.mpsStart, on.mpsEnd, lastOff.mpsStart, mpsEnd);
    //printf("    Entry limits: %i-%i, %i-%i, %i-%i\n", entryStart, firstOff.entryEnd, on.entryStart, on.entryEnd, lastOff.entryStart, entryEnd);
    //printf("    First Off Counts: %i, %i, %i\n", firstOffBeamOnEntry, firstOffBeamOffEntry, firstOffBeamUnkEntry);
    //printf("    Last Off Counts: %i, %i, %i\n", lastOffBeamOnEntry, lastOffBeamOffEntry, lastOffBeamUnkEntry);
    return cutFlag;
  }

  void reset(){
    cycleNum = -1; //adjPedFirstOff = 0; adjPedOn = 0; adjPedLastOff = 0;
    //pedRMSFirstOff = 0; pedRMSOn = 0; pedRMSLastOff = 0;
    firstOffBeamOnEntry = 0; firstOffBeamOffEntry = 0; firstOffBeamUnkEntry = 0;
    lastOffBeamOnEntry = 0; lastOffBeamOffEntry = 0; lastOffBeamUnkEntry = 0;
    entryStart = 0; entryEnd = 0; mpsStart = 0; mpsEnd = 0;
  }
} LaserCycle;

vector<vector<int>> readEventCuts(int runNum){
  vector<vector<int>> eventCuts;
  ifstream mapfile(Form("%s/compmon_event_cuts_%i.map", getenv("COMPMON_MAPS"), runNum));
  if(!mapfile.good()){return eventCuts;}
  string readStr;
  while(getline(mapfile, readStr)){
    vector<int> cutRegion; stringstream ss(readStr);
    for(int i; ss >> i;){
      cutRegion.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    eventCuts.push_back(cutRegion);
  }
  mapfile.close();
  printf("Applying eventcuts file for run %i\n", runNum);
  return eventCuts;
}

vector<Int_t> balanceLaserOffs(TTree *quartetwise, LaserCycle cyc){
  Int_t firstCur = cyc.firstOff.entryEnd; Int_t lastCur = cyc.lastOff.entryStart;

  Int_t beamState, laserState, firstMPS;

  quartetwise->SetBranchAddress("beamState", &beamState);
  quartetwise->SetBranchAddress("laserState", &laserState);
  quartetwise->SetBranchAddress("firstMPSnumber", &firstMPS);
  Int_t firstValidEnts = 0; Int_t lastValidEnts = 0;
  Int_t selEntry = 0;
  Int_t fBeamOn = 0; Int_t fBeamOff = 0; Int_t fBeamUnk = 0;
  Int_t lBeamOn = 0; Int_t lBeamOff = 0; Int_t lBeamUnk = 0;
  
  while(firstCur > cyc.firstOff.entryStart && lastCur < cyc.lastOff.entryEnd){
    if(firstValidEnts > lastValidEnts){quartetwise->GetEntry(lastCur); selEntry = lastCur++;}
    else{quartetwise->GetEntry(firstCur); selEntry = firstCur--;}
    if(selEntry <= cyc.firstOff.entryEnd && (laserState==2 || laserState==3)){
      if(beamState==1){firstValidEnts++; fBeamOn++;}
      else if(beamState==0){fBeamOff++;}
      else{fBeamUnk++;}
    }
    else if(selEntry >= cyc.lastOff.entryStart && (laserState==2 || laserState==3)){
      if(beamState==1){lastValidEnts++; lBeamOn++;}
      else if(beamState==0){lBeamOff++;}
      else{lBeamUnk++;}
    }
    else if(selEntry > cyc.firstOff.entryEnd && selEntry < cyc.lastOff.entryStart){
      printf("  Wrong event range specified on cycle balancing!\n");
      exit(1);
    }
    //printf("      Loop Summary: Sel: %i (MPS %i), First: %i, Last: %i; First Count: %i, Last Count: %i\n", selEntry, firstMPS, firstCur, lastCur, firstValidEnts, lastValidEnts);
  }
  
  vector<Int_t> cycParams;
  quartetwise->GetEntry(firstCur); cycParams.push_back(firstCur); cycParams.push_back(firstMPS); 
  quartetwise->GetEntry(lastCur); cycParams.push_back(lastCur); cycParams.push_back(firstMPS);
  cycParams.push_back(fBeamOn); cycParams.push_back(fBeamOff); cycParams.push_back(fBeamUnk);
  cycParams.push_back(lBeamOn); cycParams.push_back(lBeamOff); cycParams.push_back(lBeamUnk);
  //cyc.firstOffBeamOnEntry = fBeamOn; cyc.firstOffBeamOffEntry = fBeamOff; cyc.firstOffBeamUnkEntry = fBeamUnk;
  //cyc.lastOffBeamOnEntry = lBeamOn; cyc.lastOffBeamOffEntry = lBeamOff; cyc.lastOffBeamUnkEntry = lBeamUnk;
  //cyc.entryStart = firstCur; cyc.mpsStart = fMPS; cyc.entryEnd = lastCur; cyc.mpsEnd = lMPS;
  //printf("  Balance Summary:\n");
  //printf("    Found %i first entries, and %i last entries\n", firstValidEnts, lastValidEnts);
  //printf("    MPS Limits: %i-%i, %i-%i, %i-%i\n", cyc.mpsStart, cyc.firstOff.mpsEnd, cyc.on.mpsStart, cyc.on.mpsEnd, cyc.lastOff.mpsStart, cyc.mpsEnd);
  //printf("    Entry Limits: %i-%i, %i-%i, %i-%i\n", cyc.entryStart, cyc.firstOff.entryStart, cyc.on.entryStart, cyc.on.entryEnd, cyc.lastOff.entryStart, cyc.entryEnd);
  return cycParams;
}

#endif
