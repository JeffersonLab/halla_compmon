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

void cyclewisePolPlot(Int_t minRun, Int_t maxRun){
  TFile *grand = new TFile(Form("%s/prexGrandCompton.root", getenv("COMPMON_GRAND")), "READ");
  TTree *cyc = (TTree *)grand->Get("cyc");

  TCanvas *c = new TCanvas("c", "Pol Canvas", 1500, 800);
  //c->Divide(1, 2);
  
  FitPolVar pol, asym0LasOff1, asym0LasOff2;
  Int_t runNum; Int_t cycNum;
  cyc->SetBranchAddress("Pol0", &pol);
  cyc->SetBranchAddress("Asym0LasOff1", &asym0LasOff1);
  cyc->SetBranchAddress("Asym0LasOff2", &asym0LasOff2);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  
  gStyle->SetOptStat(0);
  Int_t nEntries = (Int_t)cyc->GetEntries(Form("runNum>=%i && runNum<=%i", minRun, maxRun));
  Int_t nBin = 1;
  TH1F *hPol = new TH1F("hPol", "Polarization by Cycle", nEntries, 0, nEntries);
  TH1F *hAsymOff1 = new TH1F("hAsymOff1", "Laser Off Asymmetry", nEntries, 0, nEntries);
  TH1F *hAsymOff2 = new TH1F("hAsymOff2", "Laser Off Asymmetry", nEntries, 0, nEntries);
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(runNum < minRun || runNum > maxRun) continue;
    hPol->SetBinContent(nBin + 1, 100*TMath::Abs(pol.mean)); hPol->SetBinError(nBin + 1, 100*TMath::Abs(pol.meanErr));
    hAsymOff1->SetBinContent(nBin + 1, 1000*asym0LasOff1.mean); hAsymOff1->SetBinError(nBin + 1, 1000*asym0LasOff1.meanErr);
    hAsymOff2->SetBinContent(nBin + 1, 1000*asym0LasOff2.mean); hAsymOff2->SetBinError(nBin + 1, 1000*asym0LasOff2.meanErr);
    if(cycNum==1){
      hPol->GetXaxis()->SetBinLabel(nBin + 1, Form("%i", runNum));
      hAsymOff1->GetXaxis()->SetBinLabel(nBin + 1, Form("%i", runNum));
      hAsymOff2->GetXaxis()->SetBinLabel(nBin + 1, Form("%i", runNum));
    }
    printf("%i.%i\t%.3f +/- %.3f\t%.3f +/- %.3f\t%.3f +/- %.3f\n", runNum, cycNum, 100*pol.mean, 100*pol.meanErr,
            1000*asym0LasOff1.mean, 1000*asym0LasOff1.meanErr, 1000*asym0LasOff1.mean, 1000*asym0LasOff2.meanErr);
    nBin++;
  }

  hPol->GetXaxis()->SetTitle("runNum");
  hAsymOff1->GetXaxis()->SetTitle("runNum");
  hAsymOff2->GetXaxis()->SetTitle("runNum");
  hPol->GetYaxis()->SetTitle("|pol| (pct)");
  hAsymOff1->GetYaxis()->SetTitle("Las Off Asym (ppt)");
  hAsymOff2->GetYaxis()->SetTitle("Las Off Asym (ppt)");
  hAsymOff2->SetLineColor(kRed);
  //c->cd(1); hPol->Draw();
  hAsymOff1->Draw(); hAsymOff2->Draw("same");
}
