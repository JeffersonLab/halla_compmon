#include "../grandOnline/makePlots.h"

using namespace std;

void bpmDiffs(Int_t prexOrCrex, Int_t snailMode=1){
  TString expt = experimentCode(prexOrCrex);
  TString id = (snailMode == 0) ? "Run" : "Snail";

  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *tree = (TTree *)f->Get("snl");
  if(snailMode == 0){
    tree = (TTree *)f->Get("run");
  }
  TTree *cyc = (TTree *)f->Get("cyc");

  Int_t snailNum, runNum, cycNum, cycCut, sign;
  FitPolVar pol0;
  Int_t treeNum, numCycs;
  
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("Pol0", &pol0);
  TString varName = (snailMode == 0) ? "runNum" : "snailNum";
  tree->SetBranchAddress(varName.Data(), &treeNum);
  tree->SetBranchAddress("numCyclesAcc", &numCycs);

  const Int_t nPlots = 4;
  TGraphErrors *gDiffs[nPlots];
  TGraph *gChiSqs[nPlots];
  TString dataName[nPlots] = {"diff_bpmAx", "diff_bpmAy", "diff_bpmBx", "diff_bpmBy"};
  DataVar diffs[nPlots];
  Float_t xmins[nPlots]; Float_t xmaxs[nPlots];
  Float_t ymins[nPlots]; Float_t ymaxs[nPlots];
  Float_t yminsChi[nPlots]; Float_t ymaxsChi[nPlots];

  for(Int_t i = 0; i < nPlots; i++){
    cyc->SetBranchAddress(dataName[i].Data(), &(diffs[i]));
    gDiffs[i] = new TGraphErrors(); gDiffs[i]->SetName(Form("g_%s", dataName[i].Data()));
    gDiffs[i]->SetTitle(Form("%s by %s", dataName[i].Data(), id.Data()));
    gDiffs[i]->GetXaxis()->SetTitle(Form("%sNum", id.Data()));
    gDiffs[i]->GetYaxis()->SetTitle(Form("%s [nm]", dataName[i].Data()));
    gDiffs[i]->SetMarkerStyle(31);
    gChiSqs[i] = new TGraph(); gChiSqs[i]->SetName(Form("gChi_%s", dataName[i].Data()));
    gChiSqs[i]->SetTitle(Form("%s (Chi2 / NDF) by %s", dataName[i].Data(), id.Data()));
    gChiSqs[i]->GetXaxis()->SetTitle(Form("%sNum", id.Data()));
    gChiSqs[i]->GetYaxis()->SetTitle(Form("%s [Chi2 / NDF]", dataName[i].Data()));
    gChiSqs[i]->SetMarkerStyle(31);
    xmins[i] = 1e16; xmaxs[i] = -1e16;
    ymins[i] = 1e16; ymaxs[i] = -1e16;
    yminsChi[i] = 1e16; ymaxsChi[i] = -1e16;
  }

  const Int_t nEnts = (Int_t)tree->GetEntries();
  TH1F *hists[nEnts][nPlots];
  TF1 *fits[nEnts][nPlots];

  Int_t entryCount = 0;
  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(numCycs < 2){
      printf("Skipping Entry %i\n", treeNum);
      continue;
    }
    for(Int_t j = 0; j < nPlots; j++){
      hists[i][j] = new TH1F(Form("hEntry%i_%s", treeNum, dataName[j].Data()), 
                             Form("%s %i, %s", id.Data(), treeNum, dataName[j].Data()), 
                             numCycs, 0, numCycs);
      fits[i][j] = new TF1(Form("fEntry%i_%s", treeNum, dataName[j].Data()), "pol0");
    }
    
    Int_t binNum = 1;
    for(Int_t j = 0; j < cyc->GetEntries(); j++){
      cyc->GetEntry(j);
      
      Bool_t entryCheck = (snailMode == 0) ? treeNum != runNum : treeNum != snailNum;
      if(entryCheck || cycCut != 0) continue;
      for(Int_t k = 0; k < nPlots; k++){
        hists[i][k]->SetBinContent(binNum, diffs[k].mean);
        hists[i][k]->SetBinError(binNum, diffs[k].meanErr);
        
      }
      binNum++;
    }

    for(Int_t j = 0; j < nPlots; j++){
      hists[i][j]->Fit(fits[i][j], "Q", "", 0, numCycs);
      Float_t yVal = 1e6/2.*fits[i][j]->GetParameter(0);
      Float_t yErr = 1e6/2.*fits[i][j]->GetParError(0);
      Float_t chiSq = fits[i][j]->GetChisquare()*1.0/fits[i][j]->GetNDF();
      //Float_t yVal = fits[i][j]->GetChisquare()*1.0/fits[i][j]->GetNDF();
      //Float_t yErr = 0.0;
      gDiffs[j]->SetPoint(entryCount, treeNum, yVal);
      gDiffs[j]->SetPointError(entryCount, 0.0, yErr);
      gChiSqs[j]->SetPoint(entryCount, treeNum, chiSq);
      xmins[j] = (treeNum - 1 < xmins[j]) ? treeNum - 1 : xmins[j];
      xmaxs[j] = (treeNum + 1 > xmaxs[j]) ? treeNum + 1 : xmaxs[j];
      ymins[j] = (yVal - yErr < ymins[j]) ? yVal - yErr : ymins[j];
      ymaxs[j] = (yVal + yErr > ymaxs[j]) ? yVal + yErr : ymaxs[j];
      yminsChi[j] = (chiSq < yminsChi[j]) ? chiSq : yminsChi[j];
      ymaxsChi[j] = (chiSq > ymaxsChi[j]) ? chiSq : ymaxsChi[j];
    }
    entryCount++;
  }

  TCanvas *cDiffs = new TCanvas("cDiffs", "BPM Diffs Canvas", 1200, 800);
  TCanvas *cChiSq = new TCanvas("cChiSq", "BPM Chi2 Canvas", 1200, 800);
  cDiffs->Divide(2, 2);
  cChiSq->Divide(2, 2);

  for(Int_t i = 0; i < nPlots; i++){
    cDiffs->cd(i+1)->SetGridx();
    cDiffs->cd(i+1)->SetGridy();
    gDiffs[i]->Draw("ap");
    cChiSq->cd(i+1)->SetGridx();
    cChiSq->cd(i+1)->SetGridy();
    gChiSqs[i]->Draw("ap");
  }
}

void bpmDiffs_alt2(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);

  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  Int_t snailNum, runNum, cycNum, cycCut, sign;
  FitPolVar pol0;
  
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("Pol0", &pol0);

  const Int_t nPlots = 4;
  TGraphErrors *gDiffs[nPlots];
  TString dataName[nPlots] = {"diff_bpmAx", "diff_bpmAy", "diff_bpmBx", "diff_bpmBy"};
  DataVar diffs[nPlots];
  Float_t xmins[nPlots]; Float_t xmaxs[nPlots];
  Float_t ymins[nPlots]; Float_t ymaxs[nPlots];

  for(Int_t i = 0; i < nPlots; i++){
    cyc->SetBranchAddress(dataName[i].Data(), &(diffs[i]));
    gDiffs[i] = new TGraphErrors(); gDiffs[i]->SetName(Form("g_%s", dataName[i].Data()));
    gDiffs[i]->SetTitle(Form("%s vs Cycle Pol0", dataName[i].Data()));
    gDiffs[i]->GetXaxis()->SetTitle(Form("%s [nm]", dataName[i].Data()));
    gDiffs[i]->GetYaxis()->SetTitle("Pol0 [pct]");
    gDiffs[i]->SetMarkerStyle(31);
    xmins[i] = 1e16; xmaxs[i] = -1e16;
    ymins[i] = 1e16; ymaxs[i] = -1e16;
  }

  Int_t nPts = 0;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(cycCut != 0 || snailNum < 5) continue;
    for(Int_t j = 0; j < nPlots; j++){
      gDiffs[j]->SetPoint(nPts, 1e6*sign*diffs[j].mean, 100*sign*pol0.mean);
      gDiffs[j]->SetPointError(nPts, 1e6*diffs[j].meanErr, 100*pol0.meanErr);
    }
    nPts++;
  }

  TCanvas *cDiffs = new TCanvas("cDiffs", "BPM Diffs Canvas", 1200, 800);
  cDiffs->Divide(2, 2);

  for(Int_t i = 0; i < nPlots; i++){
    cDiffs->cd(i+1)->SetGridx();
    cDiffs->cd(i+1)->SetGridy();
    gDiffs[i]->Draw("ap");
  }
}
