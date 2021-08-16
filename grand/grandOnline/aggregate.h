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

#include "../grandConstruction/vars.h"

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
vector<vector<TLine *>> pol0Lines;
vector<Float_t> pol0Min, pol0Max;
TString polNames4[nPols] = {"Pol4", "Asym4", "Asym4NGC", "Asym4LasOn", "Asym4LasOff1", "Asym4LasOff2", "Asym4LasOnAlt", "SigSubSum4"};
Bool_t pol4Pcts[nPols] = {true, true, true, true, true, true, false, false};
Float_t pol4Mults[nPols] = {100.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1.0};
vector<PolVar> pols4; vector<TF1 *> fits4;
vector<TH1F *> pol4Hists; vector<TH1F *> pol4Pulls;
vector<TH1F *> pol4Pull2Ds; vector<TPaveText *> pol4Texts;
vector<vector<TLine *>> pol4Lines;
vector<Float_t> pol4Min, pol4Max;

const Int_t nAccVars = 20;
TString varNames0[nAccVars] = {"Acc0LasOn", "Acc0LasOff1", "Acc0LasOff2", "Acc0BeamOff",
                               "PosAcc0LasOn",   "NegAcc0LasOn",   "DiffAcc0LasOn",   "SumAcc0LasOn",
                               "PosAcc0LasOff1", "NegAcc0LasOff1", "DiffAcc0LasOff1", "SumAcc0LasOff1",
                               "PosAcc0LasOff2", "NegAcc0LasOff2", "DiffAcc0LasOff2", "SumAcc0LasOff2",
                               "PosAcc0BeamOff", "NegAcc0BeamOff", "DiffAcc0BeamOff", "SumAcc0BeamOff"};
vector<DataVar> vars0;
vector<TH1F *> var0Hists;
vector<vector<TLine *>> var0Lines;
vector<Float_t> var0Min, var0Max;
TString varNames4[nAccVars] = {"Acc4LasOn", "Acc4LasOff1", "Acc4LasOff2", "Acc4BeamOff",
                               "PosAcc4LasOn",   "NegAcc4LasOn",   "DiffAcc4LasOn",   "SumAcc4LasOn",
                               "PosAcc4LasOff1", "NegAcc4LasOff1", "DiffAcc4LasOff1", "SumAcc4LasOff1",
                               "PosAcc4LasOff2", "NegAcc4LasOff2", "DiffAcc4LasOff2", "SumAcc4LasOff2",
                               "PosAcc4BeamOff", "NegAcc4BeamOff", "DiffAcc4BeamOff", "SumAcc4BeamOff"};
vector<DataVar> vars4;
vector<TH1F *> var4Hists;
vector<vector<TLine *>> var4Lines;
vector<Float_t> var4Min, var4Max;

const Int_t nVars = 23;
TString varNames[nVars] = {"LaserPower", "BeamCurrent", "bpmAx", "bpmAy", "bpmBx", "bpmBy", 
                           "CentralRateLasOn", "CentralRateLasOff", "USbg1", "USbg2","DSbg1", "DSbg2", 
                           "HFingerLasOn", "HFingerLasOff", "VFingerLasOn", "VFingerLasOff",
                           "diff_bpmAx", "diff_bpmAy", "diff_bpmBx", "diff_bpmBy",
                           "AsymBCMLasOn", "AsymBCMLasOff1", "AsymBCMLasOff2"};
vector<DataVar> vars;
vector<TH1F *> varHists;
vector<vector<TLine *>> varLines;
vector<Float_t> varMin, varMax;

const Int_t nBursts = 20;
TString burstNames[nBursts] = {"BurstPosAcc0LasOn",   "BurstNegAcc0LasOn",   "BurstDiffAcc0LasOn",   "BurstSummAcc0LasOn",   "BurstAsym0LasOn",
                               "BurstPosAcc0LasOff1", "BurstNegAcc0LasOff1", "BurstDiffAcc0LasOff1", "BurstSummAcc0LasOff1", "BurstAsym0LasOff1",
                               "BurstPosAcc0LasOff2", "BurstNegAcc0LasOff2", "BurstDiffAcc0LasOff2", "BurstSummAcc0LasOff2", "BurstAsym0LasOff2",
                               "BurstPosAcc0LasOff",  "BurstNegAcc0LasOff",  "BurstDiffAcc0LasOff",  "BurstSummAcc0LasOff",  "BurstAsym0LasOff"};
vector<FitPolVar> burstVars;
vector<TH1F *> burstHists;
vector<vector<TLine *>> burstLines;
vector<Float_t> burstMin, burstMax;

const Int_t nCombos = 3;
TString comboNames[nCombos] = {"BurstAsym0NGC", "BurstAsym0", "BurstPol0"};
Bool_t comboPcts[nCombos] = {true, true, true};
Float_t comboMults[nCombos] = {1000.0, 1000.0, 100.0};
vector<PolVar> comboVars; vector<TF1 *> comboFits;
vector<TH1F *> comboHists; vector<TH1F *> comboPulls;
vector<TH1F *> comboPull2Ds; vector<TPaveText *> comboTexts;
vector<vector<TLine *>> comboLines;
vector<Float_t> comboMin, comboMax;

vector<Bool_t> isPol;

Int_t nCycles;
Int_t nCyclesCut;

TString exptName(Int_t snailNum){
  if(snailNum < 100 || snailNum == 500) return "prex";
  else return "crex";
}

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
    vector<TLine *> lines;
    pol0Lines.push_back(lines);
    pol0Min.push_back(1e16); pol0Max.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, varNames0[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail%i: %s RMS vs Cycle", snlNum, varNames0[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames0[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    var0Hists.push_back(h);
    vector<TLine *> lines;
    var0Lines.push_back(lines);
    var0Min.push_back(1e16); var0Max.push_back(-1e16);
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
    vector<TLine *> lines;
    pol4Lines.push_back(lines);
    pol4Min.push_back(1e16); pol4Max.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, varNames4[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s RMS vs Cycle", snlNum, varNames4[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames4[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    var4Hists.push_back(h);
    vector<TLine *> lines;
    var4Lines.push_back(lines);
    var4Min.push_back(1e16); var4Max.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nVars; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, varNames[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s RMS vs Cycle", snlNum, varNames[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    varHists.push_back(h);
    vector<TLine *> lines;
    varLines.push_back(lines);
    varMin.push_back(1e16); varMax.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nBursts; i++){
    TString hName = Form("Snail %i: %s vs Cycle", snlNum, burstNames[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s Chi2 / NDF vs Cycle", snlNum, burstNames[(Int_t)i/2].Data()); hNameAdd = "Chi2";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), burstNames[(Int_t)i/2].Data()), hName.Data(), nCycles, 0, nCycles);
    burstHists.push_back(h);
    vector<TLine *> lines;
    burstLines.push_back(lines);
    burstMin.push_back(1e16); burstMax.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nCombos; i++){
    TString cutAdd("");
    if(i % 2 == 1) cutAdd = "NoCut";
    TH1F *h; TH1F *hPull2D;
    if(i % 2 == 0){h = new TH1F(Form("h%s_%s", cutAdd.Data(), comboNames[(Int_t)(i/2)].Data()), Form("Snail %i: %s vs Cycle", snlNum, comboNames[(Int_t)(i/2)].Data()), nCyclesCut, 0, nCyclesCut);}
    else{h = new TH1F(Form("h%s_%s", cutAdd.Data(), comboNames[(Int_t)(i/2)].Data()), Form("Snail %i: %s vs Cycle (No CycleCut)", snlNum, comboNames[(Int_t)(i/2)].Data()), nCycles, 0, nCycles);}
    TH1F *hPull = new TH1F(Form("hPull%s_%s", cutAdd.Data(), comboNames[(Int_t)(i/2)].Data()), Form("%s Pull Plot", comboNames[(Int_t)(i/2)].Data()), 40, -8, 8);
    if(i % 2 == 0){hPull2D = new TH1F(Form("hPull2D%s_%s", cutAdd.Data(), comboNames[(Int_t)(i/2)].Data()), "", nCyclesCut, 0, nCyclesCut);}
    else{hPull2D = new TH1F(Form("hPull2D%s_%s", cutAdd.Data(), comboNames[(Int_t)(i/2)].Data()), "", nCycles, 0, nCycles);}
    hPull2D->SetLineColor(kGreen); hPull2D->SetFillColor(kGreen);
    comboHists.push_back(h); comboPulls.push_back(hPull); comboPull2Ds.push_back(hPull2D);
    vector<TLine *> lines;
    comboLines.push_back(lines);
    comboMin.push_back(1e16); comboMax.push_back(-1e16);
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
  for(Int_t i = 0; i < 2*nCombos; i++){
    TString fName0 = Form("f_%s", comboNames[(Int_t)i/2].Data());
    TF1 *f0 = new TF1(fName0.Data(), "pol0");
    comboFits.push_back(f0);

    TPaveText *pt0 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
    comboTexts.push_back(pt0);
  }
}

void initBranches(TTree* cyc){
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycleNum);
  cyc->SetBranchAddress("firstOffStartMPS", &fOffMPS1);
  cyc->SetBranchAddress("firstOffEndMPS", &fOffMPS2);
  cyc->SetBranchAddress("onStartMPS", &onMPS1);
  cyc->SetBranchAddress("onEndMPS", &onMPS2);
  cyc->SetBranchAddress("lastOffStartMPS", &lOffMPS1);
  cyc->SetBranchAddress("lastOffEndMPS", &lOffMPS2);
  cyc->SetBranchAddress("cycleTime", &cycleTime);
  cyc->SetBranchAddress("AnalyzingPower", &anPow);
  cyc->SetBranchAddress("CycleCut", &cycCut);

  for(Int_t i = 0; i < nPols; i++){PolVar pol; pols0.push_back(pol);}
  for(Int_t i = 0; i < nAccVars; i++){DataVar data; vars0.push_back(data);}
  for(Int_t i = 0; i < nPols; i++){PolVar pol; pols4.push_back(pol);}
  for(Int_t i = 0; i < nAccVars; i++){DataVar data; vars4.push_back(data);}
  for(Int_t i = 0; i < nVars; i++){DataVar data; vars.push_back(data);}
  for(Int_t i = 0; i < nBursts; i++){FitPolVar data; burstVars.push_back(data);}
  for(Int_t i = 0; i < nCombos; i++){PolVar data; comboVars.push_back(data);}

  for(Int_t i = 0; i < nPols; i++){cyc->SetBranchAddress(polNames0[i].Data(), &pols0[i]);}
  for(Int_t i = 0; i < nAccVars; i++){cyc->SetBranchAddress(varNames0[i].Data(), &vars0[i]);}
  for(Int_t i = 0; i < nPols; i++){cyc->SetBranchAddress(polNames4[i].Data(), &pols4[i]);}
  for(Int_t i = 0; i < nAccVars; i++){cyc->SetBranchAddress(varNames4[i].Data(), &vars4[i]);}
  for(Int_t i = 0; i < nVars; i++){cyc->SetBranchAddress(varNames[i].Data(), &vars[i]);}
  for(Int_t i = 0; i < nBursts; i++){cyc->SetBranchAddress(burstNames[i].Data(), &burstVars[i]);}
  for(Int_t i = 0; i < nCombos; i++){cyc->SetBranchAddress(comboNames[i].Data(), &comboVars[i]);}
}

void setPol0Strs(Int_t i){
  Float_t par = fits0[i]->GetParameter(0);
  Float_t parErr = fits0[i]->GetParError(0);
  Float_t chi2 = fits0[i]->GetChisquare();
  Int_t ndf = fits0[i]->GetNDF();

  pol0Texts[i]->AddText(Form("--------Fit Results--------"));
  pol0Texts[i]->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  if(pol0Pcts[(Int_t)i/2])
    pol0Texts[i]->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  pol0Texts[i]->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  pol0Texts[i]->SetBorderSize(1); pol0Texts[i]->SetFillColor(0);
}

void setPol4Strs(Int_t i){
  Float_t par = fits4[i]->GetParameter(0);
  Float_t parErr = fits4[i]->GetParError(0);
  Float_t chi2 = fits4[i]->GetChisquare();
  Int_t ndf = fits4[i]->GetNDF();

  pol4Texts[i]->AddText(Form("--------Fit Results--------"));
  pol4Texts[i]->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  if(pol4Pcts[(Int_t)i/2])
    pol4Texts[i]->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  pol4Texts[i]->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  pol4Texts[i]->SetBorderSize(1); pol4Texts[i]->SetFillColor(0);
}

void setComboStrs(Int_t i){
  Float_t par = comboFits[i]->GetParameter(0);
  Float_t parErr = comboFits[i]->GetParError(0);
  Float_t chi2 = comboFits[i]->GetChisquare();
  Int_t ndf = comboFits[i]->GetNDF();

  comboTexts[i]->AddText(Form("--------Fit Results--------"));
  comboTexts[i]->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  if(comboPcts[(Int_t)i/2])
    comboTexts[i]->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  comboTexts[i]->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  comboTexts[i]->SetBorderSize(1); comboTexts[i]->SetFillColor(0);
}

void makePol0Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nPols; i++){
    TString cutAdd("");
    if(i % 2 == 1){cutAdd = "NoCut";}
    TCanvas *c = new TCanvas(Form("cPol%s_%s", cutAdd.Data(), polNames0[(Int_t)i/2].Data()), "Pol Canvas", 1200, 600);
    TPad *pPol1 = new TPad(Form("pPol1%s_%s", cutAdd.Data(), polNames0[(Int_t)i/2].Data()), "Pol Avg", 0.0, 0.3, 0.7, 1.0);
    pPol1->SetGridx(1); pPol1->SetGridy(1);
    TPad *pPol2 = new TPad(Form("pPol2%s_%s", cutAdd.Data(), polNames0[(Int_t)i/2].Data()), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
    TPad *pPol3 = new TPad(Form("pPol3%s_%s", cutAdd.Data(), polNames0[(Int_t)i/2].Data()), "Pol Pull Graph", 0.0, 0.0, 0.7, 0.3);
    pPol3->SetGridx(1); pPol3->SetGridy(1);
    pPol1->Draw(); pPol2->Draw(); pPol3->Draw();

    pPol1->cd();
    pol0Hists[i]->Fit(fits0[i], "Q", "", 0, nCycles);
    Float_t par = fits0[i]->GetParameter(0);
    for(Int_t j = 0; j < nCycles; j++){
      Float_t nDiffs = (pol0Hists[i]->GetBinContent(j + 1) - par)/pol0Hists[i]->GetBinError(j + 1);
      pol0Pull2Ds[i]->SetBinContent(j + 1, nDiffs); pol0Pulls[i]->Fill(nDiffs);
    }
    
    setPol0Strs(i);
    pol0Hists[i]->SetStats(0);
    pol0Hists[i]->Draw("P");
    for(Int_t j = 0; j < pol0Lines[i].size(); j++){
      pol0Lines[i][j]->Draw("same");
    }
    pol0Texts[i]->Draw("same");

    pPol2->cd();
    pol0Pulls[i]->SetStats(220);
    pol0Pulls[i]->Draw();
    
    pPol3->cd();
    pol0Pull2Ds[i]->SetStats(0);
    pol0Pull2Ds[i]->Draw();

    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
  //return msmt;
}

void makePol0Lines(Int_t snlNum, vector<Float_t> bounds, vector<Float_t> cutBounds){
  for(Int_t i = 0; i < pol0Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      Float_t xpos = (i % 2 == 1) ? bounds[j] : cutBounds[j];
      TLine *line = new TLine(xpos, pol0Min[i], xpos, pol0Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      pol0Lines[i].push_back(line);
    }
  }
}

void makePol4Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nPols; i++){
    TString cutAdd("");
    if(i % 2 == 0){cutAdd = "NoCut";}
    TCanvas *c = new TCanvas(Form("cPol%s_%s", cutAdd.Data(), polNames4[(Int_t)i/2].Data()), "Pol Canvas", 1200, 600);
    TPad *pPol1 = new TPad(Form("pPol1%s_%s", cutAdd.Data(), polNames4[(Int_t)i/2].Data()), "Pol Avg", 0.0, 0.3, 0.7, 1.0);
    pPol1->SetGridx(1); pPol1->SetGridy(1);
    TPad *pPol2 = new TPad(Form("pPol2%s_%s", cutAdd.Data(), polNames4[(Int_t)i/2].Data()), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
    TPad *pPol3 = new TPad(Form("pPol3%s_%s", cutAdd.Data(), polNames4[(Int_t)i/2].Data()), "Pol Pull Graph", 0.0, 0.0, 0.7, 0.3);
    pPol3->SetGridx(1); pPol3->SetGridy(1);
    pPol1->Draw(); pPol2->Draw(); pPol3->Draw();

    pPol1->cd();
    pol4Hists[i]->Fit(fits4[i], "Q", "", 0, nCycles);
    Float_t par = fits4[i]->GetParameter(0);
    for(Int_t j = 0; j < nCycles; j++){
      Float_t nDiffs = (pol4Hists[i]->GetBinContent(j + 1) - par)/pol4Hists[i]->GetBinError(j + 1);
      pol4Pull2Ds[i]->SetBinContent(j + 1, nDiffs); pol4Pulls[i]->Fill(nDiffs);
    }
    
    setComboStrs(i);
    pol4Hists[i]->SetStats(0);
    pol4Hists[i]->Draw("P");
    for(Int_t j = 0; j < pol4Lines[i].size(); j++){
      pol4Lines[i][j]->Draw("same");
    }
    pol4Texts[i]->Draw("same");

    pPol2->cd();
    pol4Pulls[i]->SetStats(220);
    pol4Pulls[i]->Draw();
    
    pPol3->cd();
    pol4Pull2Ds[i]->SetStats(0);
    pol4Pull2Ds[i]->Draw();

    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
  //return msmt;
}

void makePol4Lines(Int_t snlNum, vector<Float_t> bounds, vector<Float_t> cutBounds){
  for(Int_t i = 0; i < pol4Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      Float_t xpos = (i % 2 == 1) ? bounds[j] : cutBounds[j];
      TLine *line = new TLine(xpos, pol4Min[i], xpos, pol4Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      pol4Lines[i].push_back(line);
    }
  }
}

void makeAnPowPlots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  TString canName("cAnPow_Cut");
}

void makeVar0Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString canName = Form("cAccVar_%s", varNames0[(Int_t)i/2].Data());
    if(i % 2 == 1)
      canName = Form("cAccVar_%s_RMS", varNames0[(Int_t)i/2].Data());
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Acc Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    var0Hists[i]->SetStats(0);
    var0Hists[i]->Draw("P");
    for(Int_t j = 0; j < var0Lines[i].size(); j++){
      var0Lines[i][j]->Draw("same");
    }
    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeVar0Lines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < var0Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], var0Min[i], bounds[j], var0Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      var0Lines[i].push_back(line);
    }
  }
}

void makeVar4Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString canName = Form("cAccVar_%s", varNames4[(Int_t)i/2].Data());
    if(i % 2 == 1)
      canName = Form("cAccVar_%s_RMS", varNames4[(Int_t)i/2].Data());
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Acc Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    var4Hists[i]->SetStats(0);
    var4Hists[i]->Draw("P");
    for(Int_t j = 0; j < var4Lines[i].size(); j++){
      var4Lines[i][j]->Draw("same");
    }
    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeVar4Lines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < var4Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], var4Min[i], bounds[j], var4Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      var4Lines[i].push_back(line);
    }
  }
}

void makeVarPlots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nVars; i++){
    TString canName = Form("cAccVar_%s", varNames[(Int_t)i/2].Data());
    if(i % 2 == 1)
      canName = Form("cAccVar_%s_RMS", varNames[(Int_t)i/2].Data());
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    varHists[i]->SetStats(0);
    varHists[i]->Draw("P");
    for(Int_t j = 0; j < varLines[i].size(); j++){
      varLines[i][j]->Draw("same");
    }
    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeVarLines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < varHists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], varMin[i], bounds[j], varMax[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      varLines[i].push_back(line);
    }
  }
}

void makeBurstPlots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString canName = Form("cBurstVar_%s", burstNames[(Int_t)i/2].Data());
    if(i % 2 == 1){
      canName = Form("cBurstVar_%s_Chi2", burstNames[(Int_t)i/2].Data());
      burstHists[i]->SetMarkerStyle(3);
    }
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Burst Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    burstHists[i]->SetStats(0);
    burstHists[i]->Draw("P");
    for(Int_t j = 0; j < burstLines[i].size(); j++){
      burstLines[i][j]->Draw("same");
    }
    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeBurstLines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < burstHists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], burstMin[i], bounds[j], burstMax[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      burstLines[i].push_back(line);
    }
  }
}

void makeComboPlots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nCombos; i++){
    TString cutAdd("");
    if(i % 2 == 0){cutAdd = "NoCut";}
    TCanvas *c = new TCanvas(Form("cPol%s_%s", cutAdd.Data(), comboNames[(Int_t)i/2].Data()), "Pol Canvas", 1200, 600);
    TPad *pPol1 = new TPad(Form("pPol1%s_%s", cutAdd.Data(), comboNames[(Int_t)i/2].Data()), "Pol Avg", 0.0, 0.3, 0.7, 1.0);
    pPol1->SetGridx(1); pPol1->SetGridy(1);
    TPad *pPol2 = new TPad(Form("pPol2%s_%s", cutAdd.Data(), comboNames[(Int_t)i/2].Data()), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
    TPad *pPol3 = new TPad(Form("pPol3%s_%s", cutAdd.Data(), comboNames[(Int_t)i/2].Data()), "Pol Pull Graph", 0.0, 0.0, 0.7, 0.3);
    pPol3->SetGridx(1); pPol3->SetGridy(1);
    pPol1->Draw(); pPol2->Draw(); pPol3->Draw();

    pPol1->cd();
    comboHists[i]->Fit(comboFits[i], "Q", "", 0, nCycles);
    Float_t par = comboFits[i]->GetParameter(0);
    for(Int_t j = 0; j < nCycles; j++){
      Float_t nDiffs = (comboHists[i]->GetBinContent(j + 1) - par)/comboHists[i]->GetBinError(j + 1);
      comboPull2Ds[i]->SetBinContent(j + 1, nDiffs); comboPulls[i]->Fill(nDiffs);
    }
    
    setComboStrs(i);
    comboHists[i]->SetStats(0);
    comboHists[i]->Draw("P");
    for(Int_t j = 0; j < comboLines[i].size(); j++){
      comboLines[i][j]->Draw("same");
    }
    comboTexts[i]->Draw("same");

    pPol2->cd();
    comboPulls[i]->SetStats(220);
    comboPulls[i]->Draw();
    
    pPol3->cd();
    comboPull2Ds[i]->SetStats(0);
    comboPull2Ds[i]->Draw();

    c->Print(Form("%s/agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
  //return msmt;
}

void makeComboLines(Int_t snlNum, vector<Float_t> bounds, vector<Float_t> cutBounds){
  for(Int_t i = 0; i < comboHists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      Float_t xpos = (i % 2 == 1) ? bounds[j] : cutBounds[j];
      TLine *line = new TLine(xpos, comboMin[i], xpos, comboMax[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      comboLines[i].push_back(line);
    }
  }
}

#endif
