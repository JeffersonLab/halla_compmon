#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "../buildGrandRootfile.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

using namespace std;

void runwisePolPlot(){
  TFile *grand = new TFile(Form("%s/prexGrandCompton.root", getenv("COMPMON_GRAND")), "READ");
  TTree *run = (TTree *)grand->Get("run");

  TCanvas *c = new TCanvas("c", "Pol Canvas", 1500, 800);
  c->Divide(1, 2);
  
  FitPolVar pol, asym0LasOff;
  Int_t runNum;
  run->SetBranchAddress("Pol0", &pol);
  run->SetBranchAddress("Asym0LasOff", &asym0LasOff);
  run->SetBranchAddress("runNum", &runNum);
  
  gStyle->SetOptStat(0);
  TH1F *hPol = new TH1F("hPol", "Polarization by Run (New Calc)", run->GetEntries(), 0, run->GetEntries());
  TH1F *hAsymOff = new TH1F("hAsymOff", "Laser Off Asymmetry", run->GetEntries(), 0, run->GetEntries());
  //TH1F *hAltPol = new TH1F("hAltPol", "Polarization by Snail (Old Calc), snl->GetEntries(), 0, snl->GetEntries());"
  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    hPol->SetBinContent(i + 1, 100*TMath::Abs(pol.mean)); hPol->SetBinError(i + 1, 100*TMath::Abs(pol.meanErr));
    hAsymOff->SetBinContent(i + 1, 1000*asym0LasOff.mean); hAsymOff->SetBinError(i + 1, 1000*TMath::Abs(asym0LasOff.meanErr));
    hPol->GetXaxis()->SetBinLabel(i + 1, Form("%i", runNum));
    hAsymOff->GetXaxis()->SetBinLabel(i + 1, Form("%i", runNum));
    printf("%i\t%.3f +/- %.3f\t%.3f +/- %.3f\n", runNum, 100*pol.mean, 100*pol.meanErr,
            1000*asym0LasOff.mean, 1000*asym0LasOff.meanErr);
  }

  hPol->GetXaxis()->SetTitle("runNum");
  hAsymOff->GetXaxis()->SetTitle("runNum");
  hPol->GetYaxis()->SetTitle("|pol| (pct)");
  hAsymOff->GetYaxis()->SetTitle("Las Off Asym (ppt)");
  c->cd(1); hPol->Draw();
  c->cd(2); hAsymOff->Draw();
}
