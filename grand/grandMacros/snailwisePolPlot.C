#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "../buildGrandRootfile.h"
#include "makePlots.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

using namespace std;

void snailwisePolPlot(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *grand = new TFile(Form("%s/%sGrandCompton.root", getenv("COMPMON_GRAND"), expt.Data()), "READ");
  TTree *snl = (TTree *)grand->Get("snl");

  //TCanvas *c = new TCanvas("c", "Pol Canvas", 1500, 800);
  //c->Divide(1, 2);
  
  FitPolVar pol;
  Int_t snailNum;
  snl->SetBranchAddress("Pol0", &pol);
  snl->SetBranchAddress("snailNum", &snailNum);
  
  gStyle->SetOptStat(0);
  TH1F *hPol = new TH1F("hPol", "Polarization by Snail (New Calc)", snl->GetEntries(), 0, snl->GetEntries());
  //TH1F *hAltPol = new TH1F("hAltPol", "Polarization by Snail (Old Calc), snl->GetEntries(), 0, snl->GetEntries());"
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    hPol->SetBinContent(i + 1, 100*TMath::Abs(pol.mean)); hPol->SetBinError(i + 1, 100*TMath::Abs(pol.meanErr));
    hPol->GetXaxis()->SetBinLabel(i + 1, Form("%i", snailNum));
    printf("%i\t%.3f +/- %.3f\n", snailNum, 100*pol.mean, 100*pol.meanErr);
  }

  hPol->GetXaxis()->SetTitle("snailNum");
  hPol->GetYaxis()->SetTitle("|pol| (pct)");
  hPol->Draw();
}
