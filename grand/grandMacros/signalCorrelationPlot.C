#include <TCanvas.h>
#include <TPad.h>
#include <TObject.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TString.h>
#include <THStack.h>
#include <vector>
#include <string>

#include "../grandOnline/makePlots.h"
#include "../grandConstruction/vars.h"

void signalCorrelationPlot(Int_t prexOrCrex, Int_t runLo, Int_t runHi){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t snailNum, runNum, cycNum, cycCut, sign;
  Float_t ihwp;
  DataVar accLasOn, accLasOff1, accLasOff2;
  PolVar pol0;
  DataVar diffAx, diffAy, diffBx, diffBy;
  const Int_t nBPMs = 4;
  TString bpmNames[nBPMs] = {"diff_bpmAx", "diff_bpmAy", "diff_bpmBx", "diff_bpmBy"};
  
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("Acc0LasOn", &accLasOn);
  cyc->SetBranchAddress("Acc0LasOff1", &accLasOff1);
  cyc->SetBranchAddress("Acc0LasOff2", &accLasOff2);
  cyc->SetBranchAddress("diff_bpmAx", &diffAx);
  cyc->SetBranchAddress("diff_bpmAy", &diffAy);
  cyc->SetBranchAddress("diff_bpmBx", &diffBx);
  cyc->SetBranchAddress("diff_bpmBy", &diffBy);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("ihwp", &ihwp);
  cyc->SetBranchAddress("Asym0NGC", &pol0);
  cyc->SetBranchAddress("CycleCut", &cycCut);

  TGraphErrors *gX = new TGraphErrors();
  gX->SetMarkerStyle(8);
  TGraphErrors *gY = new TGraphErrors();
  gY->SetMarkerStyle(8);

  Int_t nPtsX = 0;
  Int_t nPtsY = 0;
  Float_t xminX = 1e10; Float_t xmaxX = -1e10;
  Float_t xminY = 1e10; Float_t xmaxY = -1e10;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    // Float_t lasOffMeanErr = TMath::Sqrt(TMath::Power(accLasOff1.meanErr, 2) + TMath::Power(accLasOff2.meanErr, 2))/2.0;
    // Float_t valMeanErr = TMath::Sqrt(TMath::Power(accLasOn.meanErr, 2) + TMath::Power(lasOffMeanErr, 2));
    // Float_t val = accLasOn.mean - (accLasOff1.mean + accLasOff2.mean)/2.0;

    Float_t vals[nBPMs] = {diffAx.mean, diffAy.mean, diffBx.mean, diffBy.mean};
    Float_t valErrs[nBPMs] = {diffAx.meanErr, diffAy.meanErr, diffBx.meanErr, diffBy.meanErr}; 
    Float_t valX = (vals[0] + vals[2])/2.0;
    Float_t valMeanErrX = TMath::Sqrt(TMath::Power(valErrs[0], 2) + TMath::Power(valErrs[2], 2))/2.0;
    Float_t valY = (vals[1] + vals[3])/2.0;
    Float_t valMeanErrY = TMath::Sqrt(TMath::Power(valErrs[1], 2) + TMath::Power(valErrs[3], 2))/2.0;
    
    if(cycCut == 0 && sign != 0 && runNum>=runLo && runNum<=runHi){
      gX->SetPoint(nPtsX, sign*valX, sign*pol0.mean);
      gX->SetPointError(nPtsX++, 0.0, pol0.meanErr);
      gY->SetPoint(nPtsY, sign*valY, sign*pol0.mean);
      gY->SetPointError(nPtsY++, 0.0, pol0.meanErr);
      if(sign*valX < xminX){xminX = sign*valX;}
      if(sign*valX > xmaxX){xmaxX = sign*valX;}
      if(sign*valY < xminY){xminY = sign*valY;}
      if(sign*valY > xmaxY){xmaxY = sign*valY;}
    }
  }

  printf("Range X: %.6f - %.6f\n", xminX, xmaxX);
  printf("Range Y: %.6f - %.6f\n", xminY, xmaxY);

  TCanvas *c = new TCanvas("cCorr", "Correlation", 1800, 600);
  // c->Divide(2, 1);
  // c->cd(1);
  // gX->SetTitle(Form("Asym0NGC vs BPM diff x avg Correlation (Runs %i-%i)", runLo, runHi));
  // gX->GetXaxis()->SetTitle("BPM diff x avg [mm]");
  // gX->GetYaxis()->SetTitle("Asym0NGC");
  // gX->Draw("ap");

  // TF1 *fitX = new TF1("fitX", "pol1");
  // gX->Fit("fitX", "Q", "", xminX, xmaxX);

  // TPaveText *ptX = new TPaveText(0.7, 0.75, 0.9, 0.9, "blNDC");
  // ptX->AddText("--------Fit Parameters--------");
  // ptX->AddText(Form("p0: %.4f +/- %.4f", fitX->GetParameter(0), fitX->GetParError(0)));
  // ptX->AddText(Form("p1: %.4f +/- %.4f", fitX->GetParameter(1), fitX->GetParError(1)));
  // ptX->AddText(Form("Chi2 / NDF: %.3f / %i", fitX->GetChisquare(), fitX->GetNDF()));
  // ptX->AddText(Form("Correlation: %.6f", gX->GetCorrelationFactor()));
  // ptX->SetFillColor(0); ptX->SetBorderSize(1);
  // ptX->Draw("same");

  // c->cd(2);
  c->cd(2);
  gY->SetTitle(Form("Asym0NGC vs BPM diff y avg Correlation (Runs %i-%i)", runLo, runHi));
  gY->GetXaxis()->SetTitle("BPM diff y avg [mm]");
  gY->GetYaxis()->SetTitle("Asym0NGC");
  gY->Draw("ap");

  TF1 *fitY = new TF1("fitY", "pol1");
  gY->Fit("fitY", "Q", "", xminY, xmaxY);

  TPaveText *ptY = new TPaveText(0.7, 0.75, 0.9, 0.9, "blNDC");
  ptY->AddText("--------Fit Parameters--------");
  ptY->AddText(Form("p0: %.4f +/- %.4f", fitY->GetParameter(0), fitY->GetParError(0)));
  ptY->AddText(Form("p1: %.4f +/- %.4f", fitY->GetParameter(1), fitY->GetParError(1)));
  ptY->AddText(Form("Chi2 / NDF: %.3f / %i", fitY->GetChisquare(), fitY->GetNDF()));
  ptY->AddText(Form("Correlation: %.6f", gY->GetCorrelationFactor()));
  ptY->SetFillColor(0); ptY->SetBorderSize(1);
  ptY->Draw("same");
}
