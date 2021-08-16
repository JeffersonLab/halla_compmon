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

#include "aggregate.h"

using namespace std;


void aggregateGrand(Int_t snlNum){
  TFile *grand = new TFile(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), exptName(snlNum).Data()), "READ");
  TTree *cyc = (TTree *)grand->Get("cyc");

  nCycles = (Int_t)cyc->GetEntries(Form("snailNum==%i", snlNum));
  nCyclesCut = (Int_t)cyc->GetEntries(Form("snailNum==%i && CycleCut==0", snlNum));
  printf("Number of cycles: %i\n", nCycles);
  initHists(snlNum);
  initFits();
  initBranches(cyc);
  
  Int_t nBin = 1; Int_t nBinCut = 1;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(snailNum != snlNum) continue;
    for(Int_t j = 0; j < 2*nPols; j+=2){
      pol0Hists[j+1]->SetBinContent(nBin, pols0[(Int_t)(j/2)].mean*pol0Mults[(Int_t)(j/2)]);
      pol0Hists[j+1]->SetBinError(nBin, pols0[(Int_t)(j/2)].meanErr*pol0Mults[(Int_t)(j/2)]);
      pol4Hists[j+1]->SetBinContent(nBin, pols4[(Int_t)(j/2)].mean*pol4Mults[(Int_t)(j/2)]);
      pol4Hists[j+1]->SetBinError(nBin, pols4[(Int_t)(j/2)].meanErr*pol4Mults[(Int_t)(j/2)]);
      if(nBin % 3 == 1){
        pol0Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        pol4Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
      if(cycCut == 0){
        pol0Hists[j]->SetBinContent(nBinCut, pols0[(Int_t)(j/2)].mean*pol0Mults[(Int_t)(j/2)]);
        pol0Hists[j]->SetBinError(nBinCut, pols0[(Int_t)(j/2)].meanErr*pol0Mults[(Int_t)(j/2)]);
        pol4Hists[j]->SetBinContent(nBinCut, pols4[(Int_t)(j/2)].mean*pol4Mults[(Int_t)(j/2)]);
        pol4Hists[j]->SetBinError(nBinCut, pols4[(Int_t)(j/2)].meanErr*pol4Mults[(Int_t)(j/2)]);
      }
      if(nBinCut % 3 == 1 && cycCut == 0){
        pol0Hists[j]->GetXaxis()->SetBinLabel(nBinCut, Form("%i.%i", runNum, cycleNum));
        pol4Hists[j]->GetXaxis()->SetBinLabel(nBinCut, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nAccVars; j+=2){
      var0Hists[j]->SetBinContent(nBin, vars0[(Int_t)j/2].mean);
      var0Hists[j]->SetBinError(nBin, vars0[(Int_t)j/2].meanErr);
      var4Hists[j]->SetBinContent(nBin, vars4[(Int_t)j/2].mean);
      var4Hists[j]->SetBinError(nBin, vars4[(Int_t)j/2].meanErr);
      var0Hists[j+1]->SetBinContent(nBin, vars0[(Int_t)j/2].rms);
      var0Hists[j+1]->SetBinError(nBin, vars0[(Int_t)j/2].rmsErr);
      var4Hists[j+1]->SetBinContent(nBin, vars4[(Int_t)j/2].rms);
      var4Hists[j+1]->SetBinError(nBin, vars4[(Int_t)j/2].rmsErr);
      if(nBin % 3 == 1){
        var0Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        var4Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        var0Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        var4Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nVars; j+=2){
      varHists[j]->SetBinContent(nBin, vars[(Int_t)j/2].mean);
      varHists[j]->SetBinError(nBin, vars[(Int_t)j/2].meanErr);
      varHists[j+1]->SetBinContent(nBin, vars[(Int_t)j/2].rms);
      varHists[j+1]->SetBinError(nBin, vars[(Int_t)j/2].rmsErr);
      if(nBin % 3 == 1){
        varHists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        varHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nBursts; j+=2){
      burstHists[j]->SetBinContent(nBin, burstVars[(Int_t)j/2].mean);
      burstHists[j]->SetBinError(nBin, burstVars[(Int_t)j/2].meanErr);
      if(burstVars[(Int_t)j/2].NDF == 0){burstHists[j+1]->SetBinContent(nBin, 0.0);}
      else{burstHists[j+1]->SetBinContent(nBin, burstVars[(Int_t)j/2].Chi2*1.0/burstVars[(Int_t)j/2].NDF);}
      if(nBin % 3 == 1){
        burstHists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        burstHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nCombos; j+=2){
      comboHists[j+1]->SetBinContent(nBin, comboVars[(Int_t)(j/2)].mean*comboMults[(Int_t)(j/2)]);
      comboHists[j+1]->SetBinError(nBin, comboVars[(Int_t)(j/2)].meanErr*comboMults[(Int_t)(j/2)]);
      if(nBin % 3 == 1){
        comboHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
      if(cycCut == 0){
        comboHists[j]->SetBinContent(nBinCut, comboVars[(Int_t)(j/2)].mean*comboMults[(Int_t)(j/2)]);
        comboHists[j]->SetBinError(nBinCut, comboVars[(Int_t)(j/2)].meanErr*comboMults[(Int_t)(j/2)]);
      }
      if(nBinCut % 3 == 1 && cycCut == 0){
        comboHists[j]->GetXaxis()->SetBinLabel(nBinCut, Form("%i.%i", runNum, cycleNum));
      }
    }
    nBin++;
    if(cycCut == 0){nBinCut++;}
  }

  Int_t msmt = 0;
  //gStyle->SetOptStat(0);
  makePol0Plots(snlNum, &msmt);
  makeVar0Plots(snlNum, &msmt);
  //makePol4Plots(snlNum, &msmt);
  //makeVar4Plots(snlNum, &msmt);
  makeVarPlots(snlNum, &msmt);
  makeBurstPlots(snlNum, &msmt);
  makeComboPlots(snlNum, &msmt);

  TString outputDir = Form("%s/snails/snail%i", getenv("COMPMON_WEB"), snlNum);
  gSystem->Exec(Form("pdfunite %s/agg_plots_*.pdf %s/snail%i_agg_plots.pdf", outputDir.Data(), outputDir.Data(), snlNum));
  gSystem->Exec(Form("rm -f %s/agg_plots_*.pdf", outputDir.Data()));
}
