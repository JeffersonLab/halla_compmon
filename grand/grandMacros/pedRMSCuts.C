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

void pedRMSCuts(){
  //TFile *f = TFile::Open(Form("%s/prexGrandCompton.root", getenv("COMPMON_GRAND")));
  TFile *f = TFile::Open("~/ajzec/grand/prexGrandCompton.root");
  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t runNum, cycNum;
  DataVar acc0Off1, acc0Off2, acc0On, pedOff1, pedOff2;
  
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("Acc0LasOff1", &acc0Off1);
  cyc->SetBranchAddress("Acc0LasOff2", &acc0Off2);
  cyc->SetBranchAddress("Acc0LasOn", &acc0On);
  cyc->SetBranchAddress("PedestalMeanFirstOff", &pedOff1);
  cyc->SetBranchAddress("PedestalMeanLastOff", &pedOff2);

  Int_t nCut = 0;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    Int_t rmsCut = (Int_t)(acc0Off1.rms > 0.66 || acc0Off2.rms > 0.66);
    Float_t rmsDiff = TMath::Abs(acc0Off1.rms - acc0Off2.rms);
    Float_t rmsErr = TMath::Sqrt(acc0Off1.rmsErr*acc0Off1.rmsErr + acc0Off2.rmsErr*acc0Off2.rmsErr);
    //Int_t rmsDiffCut = (Int_t)(rmsDiff*1.0/rmsErr > 10.0);
    Int_t rmsDiffCut = (Int_t)(rmsDiff > 0.2);
    Int_t pedDiffCut = (Int_t)(TMath::Abs(pedOff1.mean - pedOff2.mean) > 1.5);

    printf("Raw Ped Diff: %0.4f\n", TMath::Abs(pedOff1.mean));
    Int_t cutFlag = rmsCut + 0x2*rmsDiffCut + 0x4*pedDiffCut;
    if(cutFlag != 0){
      printf("%04i | Cycles %i.%i cut with flag %X\n", ++nCut, runNum, cycNum, cutFlag);
    }
  }

  TString rmsDiffStr("abs(Acc0LasOff1.rms - Acc0LasOff2.rms)");
  TString rmsErrStr("sqrt(Acc0LasOff1.rmsErr*Acc0LasOff1.rmsErr + Acc0LasOff2.rmsErr*Acc0LasOff2.rmsErr)");
  TString rmsRelStr = Form("%s/%s", rmsDiffStr.Data(), rmsErrStr.Data());
  TCanvas *cRMS = new TCanvas("cRMS", "RMS Info", 1200, 800);
  cRMS->cd();
  cyc->Draw(rmsDiffStr.Data());
}
