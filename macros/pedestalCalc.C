#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

using namespace std;

void pedestalCalc(Int_t runNum){
  TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *triggerwise = (TTree *)f->Get("triggerwise");

  Float_t sum, preSum, postSum;
  Int_t mpsCoda, beamState, laserState; Int_t nPts = 0;
  Float_t ymin = 0; Float_t ymax = 0;

  TGraph *gPed = new TGraph();
  TH2F *hPed = new TH2F("hPed", "Calculated Pedestal from Snapshots", 500, 0, triggerwise->GetEntries(), 400, 3700, 3850);

  triggerwise->SetBranchAddress("sum", &sum);
  triggerwise->SetBranchAddress("sumPre", &preSum);
  triggerwise->SetBranchAddress("sumPost", &postSum);
  triggerwise->SetBranchAddress("mpsCoda", &mpsCoda);
  triggerwise->SetBranchAddress("beamState", &beamState);
  triggerwise->SetBranchAddress("laserState", &laserState);

  for(Int_t i = 0; i < triggerwise->GetEntries(); i++){
    triggerwise->GetEntry(i);
    Float_t pedestal = (preSum + postSum)/2.0;
    bool correctPedestal = TMath::Abs(preSum - postSum)/pedestal < 0.03 && pedestal < 3900;
    bool narrowPedestal = TMath::Abs(preSum - postSum) < 50;
    bool beamOn = beamState == 1;
    bool laserOn = laserState == 0 || laserState == 1;

    if(correctPedestal && narrowPedestal && beamOn){
      gPed->SetPoint(nPts++, i, pedestal);
      hPed->Fill(i, pedestal);
    }
  }

  gStyle->SetOptStat(0);
  gPed->SetTitle("Calculated Pedestal from Snapshots");
  gPed->GetXaxis()->SetTitle("Entry"); gPed->GetYaxis()->SetTitle("Pedestal (RAU)");
  hPed->GetXaxis()->SetTitle("Entry"); hPed->GetYaxis()->SetTitle("Pedestal (RAU)");

  TCanvas *c = new TCanvas("c", "Pedestal Display", 1200, 800);
  c->Divide(1, 2);
  c->cd(1); gPed->Draw("ap");
  c->cd(2); hPed->Draw("colz");


}
