#include "../grandOnline/makePlots.h"

using namespace std;

void chargeAsymConsistency(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  Int_t runNum, cycNum, cycCut;
  DataVar aqOn, aqOff1, aqOff2;

  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("AsymBCMLasOn", &aqOn);
  cyc->SetBranchAddress("AsymBCMLasOff1", &aqOff1);
  cyc->SetBranchAddress("AsymBCMLasOff2", &aqOff2);

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(cycCut != 0) continue;
    Float_t err12 = TMath::Sqrt(TMath::Power(aqOff1.meanErr, 2) + TMath::Power(aqOn.meanErr, 2));
    Float_t err13 = TMath::Sqrt(TMath::Power(aqOff1.meanErr, 2) + TMath::Power(aqOff2.meanErr, 2));
    Float_t err23 = TMath::Sqrt(TMath::Power(aqOn.meanErr, 2) + TMath::Power(aqOff2.meanErr, 2));

    Float_t diff12 = TMath::Abs(aqOff1.mean - aqOn.mean);
    Float_t diff13 = TMath::Abs(aqOff1.mean - aqOff2.mean);
    Float_t diff23 = TMath::Abs(aqOn.mean - aqOff2.mean);

    if(diff12*1.0/err12 > 3.0 || diff13*1.0/err13 > 3.0 || diff23*1.0/err23 > 3.0){
      printf("Outlier cycle found! It's Run %i Cycle %i\n", runNum, cycNum);
    }
  }
} 
