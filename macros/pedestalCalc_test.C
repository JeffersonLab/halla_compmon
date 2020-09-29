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

void pedestalCalc_test(Int_t runNum){
  TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *triggerwise = (TTree *)f->Get("triggerwise");

  Float_t sum, preSum, postSum;
  Int_t mpsCoda, beamState, laserState; Int_t nPts = 0;
  Float_t ymin = 0; Float_t ymax = 0;
  Bool_t sumIsRandom;

  // TGraph *gPed = new TGraph();
  // TGraph *hPed = new TGraph();
  //TGraph *tPed = new TGraph();
  //TGraph *nPed = new TGraph();
  // TH1F *gPed = new TH1F("gPed","preSum laser on",400,3660,3820);
  // TH1F *hPed = new TH1F("hPed","preSum laser on",400,3660,3820);
  // TH1F *tPed = new TH1F("tPed","preSum laser on",400,3660,3820);
  // TH1F *nPed = new TH1F("nPed","preSum laser on",400,3660,3820);

  TH1F *gPed = new TH1F("gPed","postSum laser on",100,3770,3800);
  TH1F *hPed = new TH1F("hPed","preSum laser on",100,3770,3800);
  TH1F *tPed = new TH1F("tPed","preSum laser on",100,3770,3800);
  TH1F *nPed = new TH1F("nPed","postSum laser on",100,3770,3800);
  // TH1F *lPed = new TH1F("lPed","preSum laser on",120,-60,60);
  //TH1F *mPed = new TH1F("mPed","preSum laser on",120,-60,60);
  TH2F *lPed = new TH2F("lPed","postVspre", 1000,-500,1500,20,-10,10);
  TH2F *mPed = new TH2F("mPed","postVspre", 1000,-500,1500,20,-10,10);
  //  TH2F *hPed = new TH2F("hPed", "Calculated Pedestal from Snapshots", 500, 0, triggerwise->GetEntries(), 400, 3700, 3850);

  triggerwise->SetBranchAddress("sum", &sum);
  triggerwise->SetBranchAddress("sumPre", &preSum);
  triggerwise->SetBranchAddress("sumPost", &postSum);
  triggerwise->SetBranchAddress("mpsCoda", &mpsCoda);
  triggerwise->SetBranchAddress("beamState", &beamState);
  triggerwise->SetBranchAddress("laserState", &laserState);
  triggerwise->SetBranchAddress("sumIsRandom",&sumIsRandom);

  for(Int_t i = 0; i < triggerwise->GetEntries(); i++){
    triggerwise->GetEntry(i);
    Float_t pedestal = (preSum + postSum)/2.0;
    bool correctPedestal = TMath::Abs(preSum - postSum)/pedestal < 0.03 && pedestal < 3900;
    bool narrowPedestal = TMath::Abs(preSum - postSum) < 50;
    bool beamOn = beamState == 1;
    bool laserOn = laserState == 0 || laserState == 1;
    bool laserOff = laserState ==2 || laserState == 3;
    bool Randoms = sumIsRandom == 1;		 
    bool notRandoms = sumIsRandom == 0;

    if(laserOn && beamOn && Randoms){
       // gPed->SetPoint(nPts++, i, preSum);
       hPed->Fill(postSum);
       gPed->Fill(preSum);
       // lPed->Fill(postSum,preSum);
       lPed->Fill(sum,(preSum-postSum));
       
    }
    if(laserOff && beamOn && Randoms){
      tPed->Fill(preSum);
      nPed->Fill(postSum);
      // mPed->Fill(postSum,preSum);
      mPed->Fill(sum,(preSum-postSum));
    }
   }

  gStyle->SetOptStat();
  gStyle->SetStatFormat("6.6g");
  gPed->SetTitle("preSum laser on");
  hPed->SetTitle("postSum laser on");
  tPed->SetTitle("preSum laser off");
  nPed->SetTitle("postSum laser off");
  lPed->SetTitle("laser on");
  mPed->SetTitle("laser off");
  
  gPed->GetXaxis()->SetTitle("RAU"); gPed->GetYaxis()->SetTitle("");
  hPed->GetXaxis()->SetTitle("RAU"); hPed->GetYaxis()->SetTitle("");
  tPed->GetXaxis()->SetTitle("RAU"); tPed->GetYaxis()->SetTitle("");
  nPed->GetXaxis()->SetTitle("RAU"); nPed->GetYaxis()->SetTitle("");
  lPed->GetXaxis()->SetTitle("sum"); lPed->GetYaxis()->SetTitle("preSum-PostSum");
  mPed->GetXaxis()->SetTitle("sum"); mPed->GetYaxis()->SetTitle("preSum-PostSum");
    
  TCanvas *c = new TCanvas("c", "pre & post sum", 1200, 800);
  c->Divide(2, 2);
  c->cd(1); gPed->Draw("");
  c->cd(2); hPed->Draw("");
  c->cd(3); tPed->Draw("");
  c->cd(4); nPed->Draw("");

TCanvas *c1 = new TCanvas("c1", "pre & post sum", 1200, 800);
   c1->Divide(1,2);
   c1->cd(1); lPed->Draw("colz");
   c1->cd(2); mPed->Draw("colz");







}
