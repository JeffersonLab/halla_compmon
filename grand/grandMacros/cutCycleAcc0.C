#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
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

using namespace std;

void cutCycleAcc0(Int_t prexOrCrex, Int_t startRun=0, Int_t startCyc=0){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/%sGrandCompton.root", getenv("COMPMON_GRAND"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  Int_t runNum, cycNum, cycCut;
  Int_t fOffStart, fOffEnd, lOffStart, lOffEnd, onStart, onEnd;
  Float_t onMode, onOffset, offMode, offOffset;
  DataVar accOn, accOff1, accOff2;

  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("firstOffStartMPS", &fOffStart);
  cyc->SetBranchAddress("firstOffEndMPS", &fOffEnd);
  cyc->SetBranchAddress("onStartMPS", &onStart);
  cyc->SetBranchAddress("onEndMPS", &onEnd);
  cyc->SetBranchAddress("lastOffStartMPS", &lOffStart);
  cyc->SetBranchAddress("lastOffEndMPS", &lOffEnd);
  cyc->SetBranchAddress("Acc0LasOn", &accOn);
  cyc->SetBranchAddress("Acc0LasOff1", &accOff1);
  cyc->SetBranchAddress("Acc0LasOff2", &accOff2);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("Acc0OnMode", &onMode);
  cyc->SetBranchAddress("Acc0OnOffset", &onOffset);
  cyc->SetBranchAddress("Acc0OffMode", &offMode);
  cyc->SetBranchAddress("Acc0OffOffset", &offOffset);

  TCanvas *c = new TCanvas("c", "Cut Cycle Acc0 Data", 1200, 800);
  c->Divide(2, 2);

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(runNum < startRun || (runNum == startRun && cycNum < startCyc) || cycCut == 0) continue;
    TFile *runfile = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
    TTree *mpswise = (TTree *)runfile->Get("mpswise");

    TString per1 = Form("mpsCoda >= %i && mpsCoda <= %i && (laserState==2 || laserState==3) && beamState==1", fOffStart, fOffEnd);
    TString per2 = Form("mpsCoda >= %i && mpsCoda <= %i && (laserState==0 || laserState==1) && beamState==1", onStart, onEnd);
    TString per3 = Form("mpsCoda >= %i && mpsCoda <= %i && (laserState==2 || laserState==3) && beamState==1", lOffStart, lOffEnd);
    TString cyc = Form("mpsCoda >= %i && mpsCoda <= %i && laserState < 4 && beamState==1", fOffStart, lOffEnd);

    //TString hOnName = Form("hON_run%i_cyc%i", runNum, cycNum);
    //TString hOff1Name = Form("hOFF1_run%i_cyc%i", runNum, cycNum);
    //TString hOff2Name = Form("hOFF2_run%i_cyc%i", runNum, cycNum);
    TString hOnName("hON");
    TString hOff1Name("hOFF1");
    TString hOff2Name("hOFF2");

    mpswise->Draw(Form("Acc0/NAcc0>>%s", hOff1Name.Data()), per1.Data(), "goff");
    mpswise->Draw(Form("Acc0/NAcc0>>%s", hOnName.Data()), per2.Data(), "goff");
    mpswise->Draw(Form("Acc0/NAcc0>>%s", hOff2Name.Data()), per3.Data(), "goff");

    TH1F *hOff1 = (TH1F *)gDirectory->Get(hOff1Name.Data());
    TH1F *hOn = (TH1F *)gDirectory->Get(hOnName.Data());
    TH1F *hOff2 = (TH1F *)gDirectory->Get(hOff2Name.Data());

    hOff1->SetLineColor(kRed);
    hOn->SetLineColor(kGreen + 2);
    hOff2->SetLineColor(kRed + 2);

    TString off1Cut = hOff1->GetRMS() > offMode + offOffset ? "Yes" : "No";
    TString onCut = hOn->GetRMS() > onMode + onOffset ? "Yes" : "No";
    TString off2Cut = hOff2->GetRMS() > offMode + offOffset ? "Yes" : "No";

    //THStack *hs = new THStack(Form("hs_run%i_cyc%i", runNum, cycNum), Form("Run %i Cycle %i: Acc0/NAcc0", runNum, cycNum));
    THStack *hs = new THStack("hs", Form("Run %i Cycle %i: Acc0/NAcc0", runNum, cycNum));
    hs->Add(hOn); hs->Add(hOff1); hs->Add(hOff2);
    c->cd(1); hs->Draw("nostack");

    c->cd(2);
    mpswise->Draw("Acc0/NAcc0:mpsCoda", cyc.Data());

    c->cd(3);
    TPaveText *pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
    pt->SetFillColor(0); pt->SetBorderSize(0);
    pt->AddText(Form("Run %i, Cycle %i Has Been Cut", runNum, cycNum));
    pt->AddText(Form("Off 1 Mean: %.3f +/- %.3f", hOff1->GetMean(), hOff1->GetMeanError()));
    pt->AddText(Form("Off 1 RMS: %.3f +/- %.3f", hOff1->GetRMS(), hOff1->GetRMSError()));
    pt->AddText(Form("On Mean: %.3f +/- %.3f", hOn->GetMean(), hOn->GetMeanError()));
    pt->AddText(Form("On RMS: %.3f +/- %.3f", hOn->GetRMS(), hOn->GetRMSError()));
    pt->AddText(Form("Off 2 Mean: %.3f +/- %.3f", hOff2->GetMean(), hOff2->GetMeanError()));
    pt->AddText(Form("Off 2 RMS: %.3f +/- %.3f", hOff2->GetRMS(), hOff2->GetRMSError()));
    pt->Draw();

    c->cd(4);
    TPaveText *pt2 = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
    pt2->SetFillColor(0); pt2->SetBorderSize(0);
    pt2->AddText(Form("Run %i Cycle %i Cut Thresholds", runNum, cycNum));
    pt2->AddText(Form("Laser on mode: %.3f; Laser on offset: %.3f", onMode, onOffset));
    pt2->AddText(Form("Laser off mode: %.3f; Laser off offset: %.3f", offMode, offOffset));
    pt2->AddText(Form("Laser on cut thresh: %.3f", onMode + onOffset));
    pt2->AddText(Form("Laser off cut thresh: %.3f", offMode + offOffset));
    pt2->AddText(Form("Off 1 Cut? %s; On Cut? %s; Off 2 Cut? %s", off1Cut.Data(), onCut.Data(), off2Cut.Data()));
    pt2->Draw();

    c->Print(Form("%s/grandMacros/run%i_cycle%02i.png", getenv("COMPMON_GRAND"), runNum, cycNum), "png");

    runfile->Close();
  }

  f->Close();

  gSystem->Exec(Form("convert %s/grandMacros/run*_cycle*.png %s/grandMacros/all_cut_cycles.pdf", getenv("COMPMON_GRAND"), getenv("COMPMON_GRAND")));
  gSystem->Exec(Form("rm -f %s/grandMacros/run*_cycle*.png", getenv("COMPMON_GRAND")));

}
