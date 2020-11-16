#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TChain.h>
#include <TColor.h>
#include <TString.h>
#include <THStack.h>

#include <vector>
#include <fstream>

#include "../vars.h"
#include "makePlots.h"

void backAsymCorr(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/%sGrandCompton.root", getenv("COMPMON_GRAND"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  PolVar asym1, asym2;
  Int_t cycCut;

  cyc->SetBranchAddress("Asym0LasOff1", &asym1);
  cyc->SetBranchAddress("Asym0LasOff2", &asym2);
  cyc->SetBranchAddress("CycleCut", &cycCut);

  TGraphErrors *g = new TGraphErrors();
  TGraphErrors *dd = new TGraphErrors();
  TGraphErrors *ddNorm = new TGraphErrors();
  TH1F *hOff1 = new TH1F("hOff1", "Laser Off Asyms Raw", cyc->GetEntries("CycleCut==0"), 0, cyc->GetEntries("CycleCut==0"));
  TH1F *hOff2 = new TH1F("hOff2", "Laser Off Asyms Raw", cyc->GetEntries("CycleCut==0"), 0, cyc->GetEntries("CycleCut==0"));
  TH1F *hDD = new TH1F("hDD", "Laser Off Asyms Double Difference", 200, -0.02, 0.02);
  TH1F *hDDNorm = new TH1F("hDDNorm", "Laser Off Asyms Double Diff (Error Normalized)", 200, -4.0, 4.0);
  Int_t graphCount = 0; Int_t ddCount = 0;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    Bool_t nonzero1 = TMath::Abs(asym1.mean) > 2.0*asym1.meanErr;
    Bool_t nonzero2 = TMath::Abs(asym2.mean) > 2.0*asym2.meanErr;
    Float_t asymDiff = asym1.mean - asym2.mean;
    Float_t asymDiffErr = TMath::Sqrt(TMath::Power(asym1.meanErr, 2) + TMath::Power(asym2.meanErr, 2));
    Bool_t outsideErr = TMath::Abs(asymDiff) > 1.0*asymDiffErr;
    if(cycCut == 0){
      dd->SetPoint(ddCount, ddCount, asymDiff);
      //dd->SetPointError(ddCount, 0.0, asymDiffErr);
      ddNorm->SetPoint(ddCount, ddCount, asymDiff/asymDiffErr);
      hOff1->SetBinContent(ddCount+1, asym1.mean);
      hOff1->SetBinError(ddCount+1, asym1.meanErr);
      hOff2->SetBinContent(ddCount+1, asym2.mean);
      hOff2->SetBinError(ddCount+1, asym2.meanErr);
      hDD->Fill(asymDiff); hDDNorm->Fill(asymDiff/asymDiffErr);
      ddCount++;
      if(outsideErr){
        g->SetPoint(graphCount, asym1.mean, asym2.mean);
        g->SetPointError(graphCount++, asym1.meanErr, asym2.meanErr);
      }
    }
  }

  TCanvas *c = new TCanvas("c", "Correlation Canvas", 1200, 800);
  c->cd(); c->SetGridx(1); c->SetGridy(1);
  g->GetXaxis()->SetTitle("Laser Off 1 Asym");
  g->GetYaxis()->SetTitle("Laser Off 2 Asym");
  g->SetTitle("Laser Off Asymmetry Correlations");
  g->Draw("ap");

  TF1 *line = new TF1("f", "pol1", -1.0, 1.0);
  line->SetParameter(0, 0.0);
  line->SetParameter(1, 1.0);
  line->SetLineColor(kBlue);
  line->Draw("same");

  printf("%i cycles left after cuts\n", graphCount);

  TCanvas *c2 = new TCanvas("c2", "Double Difference Canvas", 1600, 800);
  c2->Divide(1, 2); c2->cd(1);
  dd->SetTitle("Laser Off Asyms Double Difference");
  dd->GetXaxis()->SetTitle("CycleNum");
  dd->GetYaxis()->SetTitle("AsymOff1 - AsymOff2");
  dd->SetMarkerStyle(5);
  dd->Draw("ap");

  c2->cd(2);
  ddNorm->SetTitle("Laser Off Asyms Double Diff (Error Normalized)");
  ddNorm->GetXaxis()->SetTitle("CycleNum");
  ddNorm->GetYaxis()->SetTitle("(AsymOff1 - AsymOff2)/(Asym Off Uncert)");
  ddNorm->SetMarkerStyle(3);
  ddNorm->Draw("ap");

  TCanvas *c3 = new TCanvas("c3", "Raw Asyms Canvas", 1600, 800);
  c3->cd();
  hOff1->SetMarkerStyle(22); hOff1->SetMarkerColor(kRed);
  hOff2->SetMarkerStyle(23); hOff2->SetMarkerColor(kRed+2);
  hOff1->GetXaxis()->SetTitle("CycleNum"); hOff1->GetYaxis()->SetTitle("Raw Asym");
  hOff2->GetXaxis()->SetTitle("CycleNum"); hOff2->GetYaxis()->SetTitle("Raw Asym");
  hOff1->SetStats(0); hOff2->SetStats(0);
  TLegend *leg = new TLegend(0.82, 0.8, 0.9, 0.9, "", "blNDC");
  leg->AddEntry(hOff1, "Laser Off 1");
  leg->AddEntry(hOff2, "Laser Off 2");
  hOff1->Draw("P"); hOff2->Draw("P && same");
  leg->Draw("same");

  TCanvas *c4 = new TCanvas("c4", "Double Diff 1D Canvas", 1200, 800);
  c4->cd();
  hDD->GetXaxis()->SetTitle("AsymOff1 - AsymOff2");
  hDD->Draw();

  TCanvas *c5 = new TCanvas("c5", "Double Diff Err 1D Canvas", 1200, 800);
  c5->cd();
  hDDNorm->GetXaxis()->SetTitle("Asym Double Diff Normalized");
  hDDNorm->Draw();
}
