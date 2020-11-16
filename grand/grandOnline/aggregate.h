#ifndef AGGREGATE_H
#define AGGREGATE_H

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "../grand/vars.h"

Int_t snailNum, runNum, cycleNum;
Int_t fOffMPS1, fOffMPS2, onMPS1, onMPS2, lOffMPS1, lOffMPS2;
Int_t cycCut;
Float_t cycleTime, anPow;

const Int_t nPols = 8;
TString polNames0[nPols] = {"Pol0", "Asym0", "Asym0NGC", "Asym0LasOn", "Asym0LasOff1", "Asym0LasOff2", "Asym0LasOnAlt", "SigSubSum0"};
Bool_t pol0Pcts[nPols] = {true, true, true, true, true, true, false, false};
Float_t pol0Mults[nPols] = {100.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1.0};
vector<PolVar> pols0; vector<TF1 *> fits0;
vector<TH1F *> pol0Hists; vector<TH1F *> pol0Pulls;
vector<TH1F *> pol0Pull2Ds; vector<TPaveText *> pol0Texts;
TString polNames4[nPols] = {"Pol4", "Asym4", "Asym4NGC", "Asym4LasOn", "Asym4LasOff1", "Asym4LasOff2", "Asym4LasOnAlt", "SigSubSum4"};
Bool_t pol4Pcts[nPols] = {true, true, true, true, true, true, false, false};
Float_t pol4Mults[nPols] = {100.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1.0};
vector<PolVar> pols4; vector<TF1 *> fits4;
vector<TH1F *> pol4Hists; vector<TH1F *> pol4Pulls;
vector<TH1F *> pol4Pull2Ds; vector<TPaveText *> pol4Texts;

const Int_t nAccVars = 20;
TString varNames0[nAccVars] = {"Acc0LasOn", "Acc0LasOff1", "Acc0LasOff2", "Acc0BeamOff",
                               "PosAcc0LasOn",   "NegAcc0LasOn",   "DiffAcc0LasOn",   "SumAcc0LasOn",
                               "PosAcc0LasOff1", "NegAcc0LasOff1", "DiffAcc0LasOff1", "SumAcc0LasOff1",
                               "PosAcc0LasOff2", "NegAcc0LasOff2", "DiffAcc0LasOff2", "SumAcc0LasOff2",
                               "PosAcc0BeamOff", "NegAcc0BeamOff", "DiffAcc0BeamOff", "SumAcc0BeamOff"};
vector<DataVar> vars0;
vector<TH1F *> var0Hists;
TString varNames4[nAccVars] = {"Acc4LasOn", "Acc4LasOff1", "Acc4LasOff2", "Acc4BeamOff",
                               "PosAcc4LasOn",   "NegAcc4LasOn",   "DiffAcc4LasOn",   "SumAcc4LasOn",
                               "PosAcc4LasOff1", "NegAcc4LasOff1", "DiffAcc4LasOff1", "SumAcc4LasOff1",
                               "PosAcc4LasOff2", "NegAcc4LasOff2", "DiffAcc4LasOff2", "SumAcc4LasOff2",
                               "PosAcc4BeamOff", "NegAcc4BeamOff", "DiffAcc4BeamOff", "SumAcc4BeamOff"};
vector<DataVar> vars4;
vector<TH1F *> var4Hists;

const Int_t nVars = 20;
TString varNames[nVars] = {"LaserPower", "BeamCurrent", "bpmAx", "bpmAy", "bpmBx", "bpmBy", 
                           "CentralRateLasOn", "CentralRateLasOff", "USbg1", "USbg2","DSbg1", "DSbg2", 
                           "HFingerLasOn", "HFingerLasOff", "VFingerLasOn", "VFingerLasOff",
                           "diff_bpmAx", "diff_bpmAy", "diff_bpmBx", "diff_bpmBy"};
vector<DataVar> vars;
vector<TH1F *> varHists;

vector<Bool_t> isPol;

Int_t nCycles;
Int_t nCyclesCut;

void initHists(Int_t snlNum){
  for(Int_t i = 0; i < 2*nPols; i++){
    TString cutAdd("");
    if(i % 2 == 1) cutAdd = "NoCut";
    TH1F *h; TH1F *hPull2D;
    if(i % 2 == 0){h = new TH1F(Form("h%s_%s", cutAdd.Data(), polNames0[(Int_t)(i/2)].Data()), Form("Snail %i: %s vs Cycle", snlNum, polNames0[(Int_t)(i/2)].Data()), nCyclesCut, 0, nCyclesCut);}
    else{h = new TH1F(Form("h%s_%s", cutAdd.Data(), polNames0[(Int_t)(i/2)].Data()), Form("Snail %i: %s vs Cycle (No CycleCut)", snlNum, polNames0[(Int_t)(i/2)].Data()), nCycles, 0, nCycles);}
    TH1F *hPull = new TH1F(Form("hPull%s_%s", cutAdd.Data(), polNames0[(Int_t)(i/2)].Data()), Form("%s Pull Plot", polNames0[(Int_t)(i/2)].Data()), 40, -8, 8);
    if(i % 2 == 0){hPull2D = new TH1F(Form("hPull2D%s_%s", cutAdd.Data(), polNames0[(Int_t)(i/2)].Data()), "", nCyclesCut, 0, nCyclesCut);}
    else{hPull2D = new TH1F(Form("hPull2D%s_%s", cutAdd.Data(), polNames0[(Int_t)(i/2)].Data()), "", nCycles, 0, nCycles);}
    hPull2D->SetLineColor(kGreen); hPull2D->SetFillColor(kGreen);
    pol0Hists.push_back(h); pol0Pulls.push_back(hPull); pol0Pull2Ds.push_back(hPull2D);
  }
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, varNames0[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail%i: %s RMS vs Cycle", snlNum, varNames0[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames0[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    var0Hists.push_back(h);
  }
  for(Int_t i = 0; i < 2*nPols; i++){
    TString cutAdd("");
    if(i % 2 == 1) cutAdd = "NoCut";
    TH1F *h; TH1F *hPull2D;
    if(i % 2 == 0){h = new TH1F(Form("h%s_%s", cutAdd.Data(), polNames4[(Int_t)(i/2)].Data()), Form("Snail %i: %s vs Cycle", snlNum, polNames4[(Int_t)(i/2)].Data()), nCyclesCut, 0, nCyclesCut);}
    else{h = new TH1F(Form("h%s_%s", cutAdd.Data(), polNames4[(Int_t)(i/2)].Data()), Form("Snail %i: %s vs Cycle (No CycleCut)", snlNum, polNames4[(Int_t)(i/2)].Data()), nCycles, 0, nCycles);}
    TH1F *hPull = new TH1F(Form("hPull%s_%s", cutAdd.Data(), polNames4[(Int_t)(i/2)].Data()), Form("%s Pull Plot", polNames4[(Int_t)(i/2)].Data()), 40, -8, 8);
    if(i % 2 == 0){hPull2D = new TH1F(Form("hPull2D%s_%s", cutAdd.Data(), polNames4[(Int_t)(i/2)].Data()), "", nCyclesCut, 0, nCyclesCut);}
    else{hPull2D = new TH1F(Form("hPull2D%s_%s", cutAdd.Data(), polNames4[(Int_t)(i/2)].Data()), "", nCycles, 0, nCycles);}
    hPull2D->SetLineColor(kGreen); hPull2D->SetFillColor(kGreen);
    pol4Hists.push_back(h); pol4Pulls.push_back(hPull); pol4Pull2Ds.push_back(hPull2D);
  }
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, varNames4[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s RMS vs Cycle", snlNum, varNames4[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames4[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    var4Hists.push_back(h);
  }
  for(Int_t i = 0; i < 2*nVars; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, varNames[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s RMS vs Cycle", snlNum, varNames[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    varHists.push_back(h);
  }
}

void initFits(){
  for(Int_t i = 0; i < 2*nPols; i++){
    TString fName0 = Form("f_%s", polNames0[(Int_t)i/2].Data());
    TString fName4 = Form("f_%s", polNames4[(Int_t)i/2].Data());
    TF1 *f0 = new TF1(fName0.Data(), "pol0");
    TF1 *f4 = new TF1(fName4.Data(), "pol0");
    fits0.push_back(f0);
    fits4.push_back(f4);

    TPaveText *pt0 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
    TPaveText *pt4 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
    pol0Texts.push_back(pt0);
    pol4Texts.push_back(pt4);
  }
}

#endif
