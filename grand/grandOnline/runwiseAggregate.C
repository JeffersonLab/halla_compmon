#include "runwiseAggregate.h"

using namespace std;


void runwiseAggregate(Int_t snlNum){
  TFile *grand = new TFile(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), exptName(snlNum).Data()), "READ");
  TTree *run = (TTree *)grand->Get("run");

  nRuns = (Int_t)run->GetEntries(Form("snailNum==%i", snlNum));
  printf("Number of runs: %i\n", nRuns);
  initHists(snlNum);
  initFits();
  initBranches(run);
  
  Int_t nBin = 1;
  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    if(snailNum != snlNum) continue;
    for(Int_t j = 0; j < nPols; j++){
      pol0Hists[j]->SetBinContent(nBin, pols0[j].mean*pol0Mults[j]);
      pol0Hists[j]->SetBinError(nBin, pols0[j].meanErr*pol0Mults[j]);
      // pol4Hists[j]->SetBinContent(nBin, pols4[j].mean*pol4Mults[j]);
      // pol4Hists[j]->SetBinError(nBin, pols4[j].meanErr*pol4Mults[j]);
      pol0Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
      // pol4Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
    }
    for(Int_t j = 0; j < 2*nAccVars; j+=2){
      var0Hists[j]->SetBinContent(nBin, vars0[(Int_t)j/2].mean);
      var0Hists[j]->SetBinError(nBin, vars0[(Int_t)j/2].meanErr);
      // var4Hists[j]->SetBinContent(nBin, vars4[(Int_t)j/2].mean);
      // var4Hists[j]->SetBinError(nBin, vars4[(Int_t)j/2].meanErr);
      var0Hists[j+1]->SetBinContent(nBin, vars0[(Int_t)j/2].rms);
      var0Hists[j+1]->SetBinError(nBin, vars0[(Int_t)j/2].rmsErr);
      // var4Hists[j+1]->SetBinContent(nBin, vars4[(Int_t)j/2].rms);
      // var4Hists[j+1]->SetBinError(nBin, vars4[(Int_t)j/2].rmsErr);
      var0Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
      // var4Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
      var0Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
      // var4Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
    }
    for(Int_t j = 0; j < 2*nVars; j+=2){
      varHists[j]->SetBinContent(nBin, vars[(Int_t)j/2].mean);
      varHists[j]->SetBinError(nBin, vars[(Int_t)j/2].meanErr);
      varHists[j+1]->SetBinContent(nBin, vars[(Int_t)j/2].rms);
      varHists[j+1]->SetBinError(nBin, vars[(Int_t)j/2].rmsErr);
      varHists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
      varHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i", runNum));
    }
    nBin++;
  }

  Int_t msmt = 0;
  //gStyle->SetOptStat(0);
  makePol0Plots(snlNum, &msmt);
  makeVar0Plots(snlNum, &msmt);
  //makePol4Plots(snlNum, &msmt);
  //makeVar4Plots(snlNum, &msmt);
  makeVarPlots(snlNum, &msmt);
  makeComboPlots(snlNum, &msmt);

  TString outputDir = Form("%s/snails/snail%i", getenv("COMPMON_WEB"), snlNum);
  gSystem->Exec(Form("pdfunite %s/run_agg_plots_*.pdf %s/snail%i_run_agg_plots.pdf", outputDir.Data(), outputDir.Data(), snlNum));
  gSystem->Exec(Form("rm -f %s/run_agg_plots_*.pdf", outputDir.Data()));
}
