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

#include "makePlots.h"
#include "../vars.h"

void signalCorrelationPlot(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t snailNum, runNum, cycNum, cycCut, sign;
  DataVar accLasOn, accLasOff1, accLasOff2;
  PolVar pol0;
  
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("Acc0LasOn", &accLasOn);
  cyc->SetBranchAddress("Acc0LasOff1", &accLasOff1);
  cyc->SetBranchAddress("Acc0LasOff2", &accLasOff2);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("Pol0", &pol0);
  cyc->SetBranchAddress("CycleCut", &cycCut);

  TGraphErrors *g = new TGraphErrors();
  g->SetMarkerStyle(8);

  Int_t nPts = 0;
  Float_t xmin = 1e10; Float_t xmax = -1e10;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    Float_t lasOffMeanErr = TMath::Sqrt(TMath::Power(accLasOff1.meanErr, 2) + TMath::Power(accLasOff2.meanErr, 2))/2.0;
    Float_t sigSubMeanErr = TMath::Sqrt(TMath::Power(accLasOn.meanErr, 2) + TMath::Power(lasOffMeanErr, 2));
    
    if(cycCut == 0){
      Float_t lasOffMean = accLasOn.mean - (accLasOff1.mean + accLasOff2.mean)/2.0;
      g->SetPoint(nPts, lasOffMean, sign*pol0.mean);
      g->SetPointError(nPts++, sigSubMeanErr, pol0.meanErr);
      if(lasOffMean < xmin){xmin = lasOffMean;}
      if(lasOffMean > xmax){xmax = lasOffMean;}
    }
  }

  TCanvas *c = new TCanvas("cCorr", "Correlation", 1200, 800);
  c->cd();
  g->SetTitle("Signal Size Pol0 Correlation");
  g->GetXaxis()->SetTitle("Acc0LasOn - Acc0LasOff [RAU]");
  g->GetYaxis()->SetTitle("Pol0");
  g->Draw("ap");

  TF1 *fit1 = new TF1("fit1", "pol1");
  g->Fit("fit1", "Q", "", xmin, xmax);

  TPaveText *pt = new TPaveText(0.7, 0.75, 0.9, 0.9, "blNDC");
  pt->AddText("--------Fit Parameters--------");
  pt->AddText(Form("p0: %.4f +/- %.4f", fit1->GetParameter(0), fit1->GetParError(0)));
  pt->AddText(Form("p1: %.4f +/- %.4f", fit1->GetParameter(1), fit1->GetParError(1)));
  pt->AddText(Form("Chi2 / NDF: %.3f / %i", fit1->GetChisquare(), fit1->GetNDF()));
  pt->SetFillColor(0); pt->SetBorderSize(1);
  pt->Draw("same");
}
