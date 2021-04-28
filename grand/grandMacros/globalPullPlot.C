#include "../grandOnline/makePlots.h"

using namespace std;

void globalPullPlot(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *snl = (TTree *)f->Get("snl");
  TTree *cyc = (TTree *)f->Get("cyc");

  const Int_t nFitRanges = 4;
  const Int_t nSnails = (Int_t)snl->GetEntries();
  Int_t rangeStarts[nFitRanges] = {101, 122, 138, 177};
  Int_t rangeEnds[nFitRanges] = {121, 137, 174, 219};
  Float_t outPar0[nFitRanges] = {86.6696, 93.5844, 91.7549, 87.7156};
  Float_t outPar1[nFitRanges] = {0.0, -0.0550, -0.0295, 0.0};
  Float_t inPar0[nFitRanges] = {86.5754, 96.4337, 90.3010, 87.3286};
  Float_t inPar1[nFitRanges] = {0.0, -0.0773, -0.0210, 0.0};
  Float_t outPar0Err[nFitRanges] = {0.0806, 1.1633, 0.4368, 0.0379};
  Float_t outPar1Err[nFitRanges] = {0.0, 0.0091, 0.0028, 0.0};
  Float_t inPar0Err[nFitRanges] = {0.0703, 1.1604, 0.4578, 0.0389};
  Float_t inPar1Err[nFitRanges] = {0.0, 0.0089, 0.0029, 0.0};
  Int_t rangeCycs[nFitRanges];
  TH1F *rangeHists[nFitRanges];
  TH1F *rangePulls[nFitRanges];
  Int_t snailCycs[nSnails];

  Int_t snlTreeNum, numCycs;
  Int_t snlNum, runNum, cycNum, sign, cycCut;
  Float_t ihwp;
  FitPolVar pol0;

  snl->SetBranchAddress("snailNum", &snlTreeNum);
  snl->SetBranchAddress("numCyclesAcc", &numCycs);
  cyc->SetBranchAddress("snailNum", &snlNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("ihwp", &ihwp);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("Pol0", &pol0);

  Int_t curRangeCycles = 0; Int_t rangeInd = 0;
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(sign == 0 || snlTreeNum == 150 || snlTreeNum == 151 || snlTreeNum == 159 || snlTreeNum == 160 || snlTreeNum == 220 || snlTreeNum == 221 || snlTreeNum == 175 || snlTreeNum == 176){
      snailCycs[i] = 0;
      continue;
    }
    if(snlTreeNum == rangeStarts[rangeInd]){
      curRangeCycles = 0;
    }
    snailCycs[i] = numCycs;
    curRangeCycles += numCycs;
    if(snlTreeNum == rangeEnds[rangeInd]){
      rangeCycs[rangeInd] = curRangeCycles;
      rangeInd++;
    }
  }
  
  for(Int_t i = 0; i < nFitRanges; i++){
    TString histName = Form("hRange_%i-%i", rangeStarts[i], rangeEnds[i]);
    TString histTitle = Form("Snails %i-%i: Cycle Pol0 Pulls", rangeStarts[i], rangeEnds[i]);
    TString pullName = Form("hPull_%i-%i", rangeStarts[i], rangeEnds[i]);
    TString pullTitle = Form("Snails %i-%i: Cycle Pol0 Pulls by Time", rangeStarts[i], rangeEnds[i]);
    rangeHists[i] = new TH1F(histName.Data(), histTitle.Data(), 100, -8.0, 8.0);
    rangeHists[i]->GetXaxis()->SetTitle("Cycle Pol0 Pull");
    rangePulls[i] = new TH1F(pullName.Data(), pullTitle.Data(), rangeCycs[i], 0, rangeCycs[i]);
    rangePulls[i]->GetYaxis()->SetTitle("Cycle Pol0 Pull");
    rangePulls[i]->SetLineColor(kGreen); 
    rangePulls[i]->SetFillColor(kGreen);
  }

  Int_t curInd = 0; Int_t histInd = 0; Int_t lastInd = -1;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(cycCut != 0 || snlNum == 150 || snlNum == 151 || snlNum == 159 || snlNum == 160 || snlNum == 220 || snlNum == 221 || snlNum == 175 || snlNum == 176) continue;
    for(Int_t j = 0; j < nFitRanges; j++){
      if(snlNum >= rangeStarts[j] && snlNum <= rangeEnds[j]){
        curInd = j;
        break;
      }
      if(curInd != lastInd){
        lastInd = curInd;
        histInd = 0;
      }
    }
    if(curInd >= nFitRanges) break;
    
    Float_t outFitVal = outPar1[curInd]*snlNum + outPar0[curInd];
    Float_t inFitVal = inPar1[curInd]*snlNum + inPar0[curInd];
    Float_t signCorrPol = 100.0*sign*pol0.mean;
    Float_t fitVal = 0.0;
    if(ihwp < 0.5) fitVal = outFitVal;
    else fitVal = inFitVal;

    Float_t pull = (signCorrPol - fitVal)/(100.0*pol0.meanErr);
    rangeHists[curInd]->Fill(pull);
    rangePulls[curInd]->SetBinContent(histInd + 1, pull);

    Int_t labelFreq = (Int_t)1.0/50.0*rangeCycs[curInd];
    if(histInd % labelFreq == 0){
      rangePulls[curInd]->GetXaxis()->SetBinLabel(histInd + 1, Form("%i.%i", runNum, cycNum));
    }

    histInd++;
  }

  for(Int_t i = 0; i < nFitRanges; i++){
    TCanvas *c = new TCanvas(Form("cRange_%i", i), Form("Canvas for Range %i-%i", rangeStarts[i], rangeEnds[i]), 1400, 800);
    c->Divide(1, 2);
    c->cd(1);
    rangeHists[i]->SetStats(220);
    rangeHists[i]->Draw();
    c->cd(2)->SetGridx();
    c->cd(2)->SetGridy();
    rangePulls[i]->SetStats(0);
    rangePulls[i]->Draw();
  }
}
