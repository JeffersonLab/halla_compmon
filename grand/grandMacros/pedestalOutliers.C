#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "../vars.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

using namespace std;

void pedestalOutliers(Int_t prexOrCrex){
  TFile *f;
  if(prexOrCrex == 1){
    f = TFile::Open(Form("%s/prexGrandCompton.root", getenv("COMPMON_GRAND")));
  }
  else if(prexOrCrex == 2){
    f = TFile::Open(Form("%s/crexGrandCompton.root", getenv("COMPMON_GRAND")));
  }
  else{
    printf("ERROR: invalid experiment code\n");
    exit(1);
  }

  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t runNum, cycNum;
  //Float_t fPedMean, fPedErr, lPedMean, lPedErr;
  DataVar fPed, lPed;

  TString fOffErr("PedestalMeanFirstOff.meanErr");
  TString lOffErr("PedestalMeanLastOff.meanErr");
  TString errStr = Form("sqrt(%s*%s + %s*%s)", fOffErr.Data(), fOffErr.Data(), lOffErr.Data(), lOffErr.Data());
  TString diffStr("(PedestalMeanFirstOff.mean - PedestalMeanLastOff.mean)");
  TString relStr = Form("%s/%s", diffStr.Data(), errStr.Data());

  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  //cyc->SetBranchAddress("PedestalMeanFirstOff/mean", &fPedMean);
  //cyc->SetBranchAddress("PedestalMeanFirstOff/meanErr", &fPedErr);
  //cyc->SetBranchAddress("PedestalMeanLastOff/mean", &lPedMean);
  //cyc->SetBranchAddress("PedestalMeanLastOff/meanErr", &lPedErr);
  cyc->SetBranchAddress("PedestalMeanFirstOff", &fPed);
  cyc->SetBranchAddress("PedestalMeanLastOff", &lPed);

  TGraphErrors *g = new TGraphErrors(); g->SetName("Err_RMS");
  g->SetTitle("Mean Err vs RMS for Pedestal Fits");

  Int_t nCut = 0; Int_t nPts = 0;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    Float_t combErr = TMath::Sqrt(fPed.meanErr*fPed.meanErr + lPed.meanErr*lPed.meanErr);
    Float_t pedDiff = fPed.mean - lPed.mean;
    Float_t invRelErr = pedDiff*1.0/combErr;
    //if(TMath::Abs(relErr) < 0.02){
    if(fPed.rms > 2.9 || lPed.rms > 2.9){
      nCut++;
      printf("%04i: Outlier Found - Cycle %04i.%02i - Ped Diff - %0.2f - F Off MeanErr: %0.4f, L Off MeanErr: %0.4f\n", 
              nCut, runNum, cycNum, pedDiff, fPed.rms, lPed.rms);
    }

    g->SetPoint(i, fPed.rms, lPed.rms);
  }

  TCanvas *c = new TCanvas("c", "Err Correlation Graph");
  c->cd();
  cyc->Draw("PedestalMeanFirstOff.mean - PedestalMeanLastOff.mean");
  //g->Draw("ap");
}
