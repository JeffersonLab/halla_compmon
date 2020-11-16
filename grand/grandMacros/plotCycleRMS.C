#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "../vars.h"

using namespace std;

void plotCycleRMS(Int_t prexOrCrex){
  TString exptStr("");
  TString exptName("");
  if(prexOrCrex == 1){
    exptStr = "prex";
    exptName = "PREX-II";
  }
  else if(prexOrCrex == 2){
    exptStr = "crex";
    exptName = "CREX";
  }
  else{
    printf("Invalid Experiment Code\n");
    exit(1);
  }
  
  TFile *f = TFile::Open(Form("%s/%sGrandCompton.root", getenv("COMPMON_GRAND"), exptStr.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t runNum, cycNum;
  DataVar on, fOff, lOff;
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum); 
  cyc->SetBranchAddress("Acc0LasOn", &on);
  cyc->SetBranchAddress("Acc0LasOff1", &fOff);
  cyc->SetBranchAddress("Acc0LasOff2", &lOff);
  
  TH1F *hF = new TH1F("hF", Form("%s Acc0 (Las Off) RMS by Cycle", exptName.Data()),
                      (Int_t)cyc->GetEntries(), 0, (Int_t)cyc->GetEntries());
  TH1F *hL = new TH1F("hL", Form("%s Acc0 (Las Off) RMS by Cycle", exptName.Data()),
                      (Int_t)cyc->GetEntries(), 0, (Int_t)cyc->GetEntries());
  TH1F *hO = new TH1F("hO", Form("%s Acc0 (Las On) RMS by Cycle", exptName.Data()),
                      (Int_t)cyc->GetEntries(), 0, (Int_t)cyc->GetEntries());

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    hF->SetBinContent(i+1, fOff.rms);
    hF->SetBinError(i+1, fOff.rmsErr);
    hL->SetBinContent(i+1, lOff.rms);
    hL->SetBinError(i+1, lOff.rmsErr);
    hO->SetBinContent(i+1, on.rms);
    hO->SetBinError(i+1, on.rmsErr);
    if(i % 4 == 0){
      hF->GetXaxis()->SetBinLabel(i+1, Form("%i.%i", runNum, cycNum));
      hL->GetXaxis()->SetBinLabel(i+1, Form("%i.%i", runNum, cycNum));
      hO->GetXaxis()->SetBinLabel(i+1, Form("%i.%i", runNum, cycNum));
    }
  }

  gStyle->SetOptStat(0);
  TCanvas *cRMS = new TCanvas("cRMS", "RMS Graph", 1200, 800);
  cRMS->SetGridx(); cRMS->SetGridy();
  TCanvas *cOnRMS = new TCanvas("cOnRMS", "RMS Graph (Las On)", 1200, 800);
  cOnRMS->SetGridx(); cOnRMS->SetGridy();
  cRMS->cd();
  hF->GetXaxis()->SetTitle("cycle num"); hL->GetXaxis()->SetTitle("cycle num");
  hF->GetYaxis()->SetTitle("Acc0 RMS"); hL->GetYaxis()->SetTitle("Acc0 RMS");
  hF->SetMarkerColor(kRed); hL->SetMarkerColor(kRed + 2);
  hF->SetLineColor(kRed); hL->SetLineColor(kRed + 2);
  hF->Draw(); hL->Draw("same");
  cOnRMS->cd();
  hO->GetXaxis()->SetTitle("cycle num"); hO->GetYaxis()->SetTitle("Acc0 RMS");
  hO->SetMarkerColor(kGreen + 2); hO->SetLineColor(kGreen + 2);
  hO->Draw();
}
