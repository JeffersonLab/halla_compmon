#include "../grandOnline/makePlots.h"

using namespace std;


void gainShiftEstimate(Int_t prexOrCrex){
  TFile *f = TFile::Open(Form("%s/aggregates/crexGrandCompton.root", getenv("COMPMON_WEB")));
  TTree *cyc = (TTree *)f->Get("cyc");

  Float_t alpha = 0.01246;
  TH1F *hCorr = new TH1F("hCorr", Form("Gain Shift Asym Correction (alpha=%.5f)", alpha), 200, -0.01, 0.01);
  FitPolVar asym0, pol0;
  DataVar diffOn, diffOff1, diffOff2, sumOn, sumOff1, sumOff2;
  Int_t sign, cycCut, snailNum;

  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("DiffAcc0LasOn", &diffOn);
  cyc->SetBranchAddress("DiffAcc0LasOff1", &diffOff1);
  cyc->SetBranchAddress("DiffAcc0LasOff2", &diffOff2);
  cyc->SetBranchAddress("SumAcc0LasOn", &sumOn);
  cyc->SetBranchAddress("SumAcc0LasOff1", &sumOff1);
  cyc->SetBranchAddress("SumAcc0LasOff2", &sumOff2);
  cyc->SetBranchAddress("Asym0", &asym0);
  cyc->SetBranchAddress("Pol0", &pol0);

  Int_t minSnail = 116;
  Int_t maxSnail = 121;
  Float_t errorSums = 0.0;

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(cycCut != 0 || snailNum<minSnail || snailNum >maxSnail) continue;

    Float_t diffOff = (diffOff1.mean + diffOff2.mean)/2.0;
    Float_t sumOff = (sumOff1.mean + sumOff2.mean)/2.0;
    Float_t f = 1.0/(sumOn.mean - sumOff);

    Float_t asymCorr = (asym0.mean + f*alpha*diffOff)/(1 + f*alpha*sumOff);
    Float_t pctDiff = (asymCorr - asym0.mean)*1.0/asym0.mean;

    hCorr->Fill(pctDiff);

    errorSums += 1.0/TMath::Power(pol0.meanErr, 2);
  }

  TCanvas *c = new TCanvas("c", "Gain Shift Correction Canvas", 1200, 800);
  c->cd();
  hCorr->GetXaxis()->SetTitle("Corr Asym Pct Diff [%]");
  hCorr->Draw();

  printf("Pol0 Uncert: %.4f\n", 100.0*TMath::Sqrt(1.0/errorSums));
}
