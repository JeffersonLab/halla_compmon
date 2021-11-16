#include "../grandOnline/makePlots.h"

using namespace std;


const Int_t nFitRanges = 4;
Int_t rangeStarts[nFitRanges] = {101, 122, 138, 177};
Int_t rangeEnds[nFitRanges] = {121, 137, 174, 221};
Float_t outPar0[nFitRanges] = {86.6820, 93.2302, 91.2236, 87.7844};
Float_t outPar1[nFitRanges] = {0.0, -0.0534, -0.0268, 0.0};
Float_t inPar0[nFitRanges] = {86.5754, 96.5867, 89.8823, 87.4129};
Float_t inPar1[nFitRanges] = {0.0, -0.0798, -0.0191, 0.0};
Float_t outPar0Err[nFitRanges] = {0.0820, 1.1912, 0.4411, 0.0383};
Float_t outPar1Err[nFitRanges] = {0.0, 0.0093, 0.0029, 0.0};
Float_t inPar0Err[nFitRanges] = {0.0715, 1.1919, 0.4661, 0.0395};
Float_t inPar1Err[nFitRanges] = {0.0, 0.0092, 0.0030, 0.0};
Int_t rangeCycs[nFitRanges];
Int_t rangeRuns[nFitRanges];
TH1F *rangeHists[nFitRanges];
TH1F *rangePulls[nFitRanges];
TH1F *rangeRunHists[nFitRanges];
TH1F *rangeRunPulls[nFitRanges];

//Correlation histogram values
Int_t nBinsXCyc = 100;
Int_t nBinsXRun = 25;
Int_t nBinsY = 100;
Float_t xLoCyc = 0.0; Float_t xHiCyc = 0.06;
Float_t xLoRun = 0.0; Float_t xHiRun = 0.06;
Float_t yLoCyc = -8.0; Float_t yHiCyc = 8.0;
Float_t yLoRun = -5.0; Float_t yHiRun = 5.0;
// TH2F *hCycOut = new TH2F("hCycOut", "", nBinsX, xLoCyc, xHiCyc, nBinsY, yLoCyc, yHiCyc);
// TH2F *hRunOut = new TH2F("hRunOut", "", nBinsX, xLoCyc, xHiCyc, nBinsY, yLoCyc, yHiCyc);
// TH2F *hCycIn  = new TH2F("hCycIn",  "", nBinsX, xLoRun, xHiRun, nBinsY, yLoRun, yHiRun);
// TH2F *hRunIn  = new TH2F("hRunIn",  "", nBinsX, xLoRun, xHiRun, nBinsY, yLoRun, yHiRun);
TProfile *hCycPos = new TProfile("hCycPos", "", nBinsXCyc, xLoCyc, xHiCyc, yLoCyc, yHiCyc);
TProfile *hRunPos = new TProfile("hRunPos", "", nBinsXRun, xLoRun, xHiRun, yLoRun, yHiRun);
TProfile *hCycNeg = new TProfile("hCycNeg", "", nBinsXCyc, xLoCyc, xHiCyc, yLoCyc, yHiCyc);
TProfile *hRunNeg = new TProfile("hRunNeg", "", nBinsXRun, xLoRun, xHiRun, yLoRun, yHiRun);
TF1 *fCycPos = new TF1("fCycPos", "pol1");
TF1 *fRunPos = new TF1("fRunPos", "pol1");
TF1 *fCycNeg = new TF1("fCycNeg", "pol1");
TF1 *fRunNeg = new TF1("fRunNeg", "pol1");
TPaveText *ptCycPos = new TPaveText(0.10, 0.70, 0.35, 0.90, "blNDC");
TPaveText *ptRunPos = new TPaveText(0.10, 0.70, 0.35, 0.90, "blNDC");
TPaveText *ptCycNeg = new TPaveText(0.10, 0.70, 0.35, 0.90, "blNDC");
TPaveText *ptRunNeg = new TPaveText(0.10, 0.70, 0.35, 0.90, "blNDC");

Int_t snlTreeNum, numCycs, snlSign, numRuns;
Int_t cycSnlNum, cycRunNum, cycNum, cycSign, cycCut;
Int_t runSnlNum, runTreeNum, runSign;
Float_t cycIHWP;
StdVar runIHWP;
FitPolVar runPol0;
PolVar cycPol0;

// Correlation Variables
FitPolVar runAsymOff;
PolVar cycAsymOff;
DataVar cycAccOn, cycAccOff1, cycAccOff2, runAccOn, runAccOff;
DataVar cycDSbg1, cycDSbg2, runDSbg1, runDSbg2;
DataVar cycRateOn, cycRateOff, runRateOn, runRateOff;
DataVar cycBPMAy, cycBPMBy, runBPMAy, runBPMBy;
DataVar cycLasPow, runLasPow;
PolVar cycCollOff;

TString corrVarName = "BPM Y Pos RMS";

void initTrees(TTree *snl, TTree *run, TTree *cyc){
  snl->SetBranchAddress("snailNum", &snlTreeNum);
  snl->SetBranchAddress("numCyclesAcc", &numCycs);
  snl->SetBranchAddress("numRuns", &numRuns);
  snl->SetBranchAddress("sign", &snlSign);

  run->SetBranchAddress("snailNum", &runSnlNum);
  run->SetBranchAddress("runNum", &runTreeNum);
  run->SetBranchAddress("sign", &runSign);
  run->SetBranchAddress("ihwp", &runIHWP);
  run->SetBranchAddress("Acc0LasOn", &runAccOn);
  run->SetBranchAddress("Acc0LasOff", &runAccOff);
  run->SetBranchAddress("bpmAy", &runBPMAy);
  run->SetBranchAddress("bpmBy", &runBPMBy);
  run->SetBranchAddress("CentralRateLasOn", &runRateOn);
  run->SetBranchAddress("CentralRateLasOff", &runRateOff);
  run->SetBranchAddress("LaserPower", &runLasPow);
  run->SetBranchAddress("DSbg1", &runDSbg1);
  run->SetBranchAddress("DSbg2", &runDSbg2);
  run->SetBranchAddress("Asym0LasOff", &runAsymOff);
  run->SetBranchAddress("Pol0", &runPol0);

  cyc->SetBranchAddress("snailNum", &cycSnlNum);
  cyc->SetBranchAddress("runNum", &cycRunNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("sign", &cycSign);
  cyc->SetBranchAddress("ihwp", &cycIHWP);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("Acc0LasOn", &cycAccOn);
  cyc->SetBranchAddress("Acc0LasOff1", &cycAccOff1);
  cyc->SetBranchAddress("Acc0LasOff2", &cycAccOff2);
  cyc->SetBranchAddress("bpmAy", &cycBPMAy);
  cyc->SetBranchAddress("bpmBy", &cycBPMBy);
  cyc->SetBranchAddress("CentralRateLasOn", &cycRateOn);
  cyc->SetBranchAddress("CentralRateLasOff", &cycRateOff);
  cyc->SetBranchAddress("LaserPower", &cycLasPow);
  cyc->SetBranchAddress("DSbg1", &cycDSbg1);
  cyc->SetBranchAddress("DSbg2", &cycDSbg2);
  cyc->SetBranchAddress("CollimatorOffset", &cycCollOff);
  cyc->SetBranchAddress("Asym0LasOff", &cycAsymOff);
  cyc->SetBranchAddress("Pol0", &cycPol0);
}


void fillCorrPlot(Float_t pull, Int_t runMode){
  Float_t cycXVal = (cycBPMAy.rms + cycBPMBy.rms)/2.0;
  Float_t runXVal = (runBPMAy.rms + runBPMBy.rms)/2.0;
  if(runMode == 0){
    if(cycSign == 1){
      hCycPos->Fill(cycXVal, pull);
    }
    else{
      hCycNeg->Fill(cycXVal, pull);
    }
  }
  else if(runMode == 1){
    if(runSign == 1){
      hRunPos->Fill(runXVal, pull);
    }
    else{
      hRunNeg->Fill(runXVal, pull);
    }
  }
}


void fillCycPullPlot(TTree *cyc, Int_t snlLo, Int_t snlHi){
  Int_t curInd = 0; Int_t histInd = 0; Int_t lastInd = -1;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(cycCut != 0 || cycSnlNum == 175 || cycSnlNum == 176) continue;
    for(Int_t j = 0; j < nFitRanges; j++){
      if(cycSnlNum >= rangeStarts[j] && cycSnlNum <= rangeEnds[j]){
        curInd = j;
        break;
      }
      if(curInd != lastInd){
        lastInd = curInd;
        histInd = 0;
      }
    }
    if(curInd >= nFitRanges) break;
    
    Float_t outFitVal = outPar1[curInd]*cycSnlNum + outPar0[curInd];
    Float_t inFitVal = inPar1[curInd]*cycSnlNum + inPar0[curInd];
    Float_t signCorrPol = 100.0*cycSign*cycPol0.mean;
    Float_t fitVal = 0.0;
    if(cycIHWP < 0.5) fitVal = outFitVal;
    else fitVal = inFitVal;

    Float_t pull = (signCorrPol - fitVal)/(100.0*cycPol0.meanErr);
    rangeHists[curInd]->Fill(pull);
    rangePulls[curInd]->SetBinContent(histInd + 1, pull);

    Int_t labelFreq = (Int_t)1.0/50.0*rangeCycs[curInd];
    if(histInd % labelFreq == 0){
      rangePulls[curInd]->GetXaxis()->SetBinLabel(histInd + 1, Form("%i.%i", cycRunNum, cycNum));
    }

    if(cycSnlNum >= snlLo && cycSnlNum <= snlHi){
      fillCorrPlot(pull, 0);
    }

    histInd++;
  }
}


void fillRunPullPlot(TTree *run, Int_t snlLo, Int_t snlHi){
  Int_t curInd = 0; Int_t histInd = 0; Int_t lastInd = -1;
  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    if(runSnlNum == 175 || runSnlNum == 176) continue;
    for(Int_t j = 0; j < nFitRanges; j++){
      if(runSnlNum >= rangeStarts[j] && runSnlNum <= rangeEnds[j]){
        curInd = j;
        if(curInd != lastInd){
          lastInd = curInd;
          histInd = 0;
        }
        break;
      }
    }
    if(curInd >= nFitRanges) break;
    
    Float_t outFitVal = outPar1[curInd]*runSnlNum+ outPar0[curInd];
    Float_t inFitVal = inPar1[curInd]*runSnlNum + inPar0[curInd];
    Float_t signCorrPol = 100.0*runSign*runPol0.mean;
    Float_t fitVal = (runIHWP.mean < 0.5) ? outFitVal : inFitVal;

    Float_t pull = (signCorrPol - fitVal)/(100.0*runPol0.meanErr);
    rangeRunHists[curInd]->Fill(pull);
    rangeRunPulls[curInd]->SetBinContent(histInd + 1, pull);
    if(TMath::Abs(pull) > 3){
      printf("Run %i has a pull of %.4f\n", runTreeNum, pull);
    }

    Int_t labelFreq = (Int_t)1.0/50.0*rangeRuns[curInd];
    if(histInd % labelFreq == 0){
      rangeRunPulls[curInd]->GetXaxis()->SetBinLabel(histInd + 1, Form("%i", runTreeNum));
    }

    if(runSnlNum >= snlLo && runSnlNum <= snlHi){
      fillCorrPlot(pull, 1);
    }

    histInd++;
  }
}


void graphCanvas(){
  Float_t loFitFrac = 1.0/4.0;
  Float_t hiFitFrac = 3.0/4.0;

  hCycPos->GetXaxis()->SetTitle(corrVarName.Data());
  hCycPos->GetYaxis()->SetTitle("Pol0 pull");
  hCycPos->SetTitle(Form("Cyclewise: %s vs Pol0 Pull (sign = +1)", corrVarName.Data()));
  hCycPos->Fit(fCycPos, "Q", "", xLoCyc + loFitFrac*(xHiCyc - xLoCyc), xLoCyc + hiFitFrac*(xHiCyc - xLoCyc));
  ptCycPos->SetFillColor(0); ptCycPos->SetBorderSize(1);
  ptCycPos->AddText("-------- Fit Params --------");
  ptCycPos->AddText(Form("p0: %.4f +/- %.4f", fCycPos->GetParameter(0), fCycPos->GetParError(0)));
  ptCycPos->AddText(Form("p1: %.4f +/- %.4f", fCycPos->GetParameter(1), fCycPos->GetParError(1)));
  ptCycPos->AddText(Form("Chi2 / NDF: %.3f / %i", fCycPos->GetChisquare(), fCycPos->GetNDF()));

  hCycNeg->GetXaxis()->SetTitle(corrVarName.Data());
  hCycNeg->GetYaxis()->SetTitle("Pol0 pull");
  hCycNeg->SetTitle(Form("Cyclewise: %s vs Pol0 Pull (sign = -1)", corrVarName.Data()));
  hCycNeg->Fit(fCycNeg, "Q", "", xLoCyc + loFitFrac*(xHiCyc - xLoCyc), xLoCyc + hiFitFrac*(xHiCyc - xLoCyc));
  ptCycNeg->SetFillColor(0); ptCycNeg->SetBorderSize(1);
  ptCycNeg->AddText("-------- Fit Params --------");
  ptCycNeg->AddText(Form("p0: %.4f +/- %.4f", fCycNeg->GetParameter(0), fCycNeg->GetParError(0)));
  ptCycNeg->AddText(Form("p1: %.4f +/- %.4f", fCycNeg->GetParameter(1), fCycNeg->GetParError(1)));
  ptCycNeg->AddText(Form("Chi2 / NDF: %.3f / %i", fCycNeg->GetChisquare(), fCycNeg->GetNDF()));
  
  hRunPos->GetXaxis()->SetTitle(corrVarName.Data());
  hRunPos->GetYaxis()->SetTitle("Pol0 pull");
  hRunPos->SetTitle(Form("Runwise: %s vs Pol0 Pull (sign = +1)", corrVarName.Data()));
  hRunPos->Fit(fRunPos, "Q", "", xLoRun + loFitFrac*(xHiRun - xLoRun), xLoRun + hiFitFrac*(xHiRun - xLoRun));
  ptRunPos->SetFillColor(0); ptRunPos->SetBorderSize(1);
  ptRunPos->AddText("-------- Fit Params --------");
  ptRunPos->AddText(Form("p0: %.4f +/- %.4f", fRunPos->GetParameter(0), fRunPos->GetParError(0)));
  ptRunPos->AddText(Form("p1: %.4f +/- %.4f", fRunPos->GetParameter(1), fRunPos->GetParError(1)));
  ptRunPos->AddText(Form("Chi2 / NDF: %.3f / %i", fRunPos->GetChisquare(), fRunPos->GetNDF()));

  hRunNeg->GetXaxis()->SetTitle(corrVarName.Data());
  hRunNeg->GetYaxis()->SetTitle("Pol0 pull");
  hRunNeg->SetTitle(Form("Runwise: %s vs Pol0 Pull (sign = -1)", corrVarName.Data()));
  hRunNeg->Fit(fRunNeg, "Q", "", xLoRun + loFitFrac*(xHiRun - xLoRun), xLoRun + hiFitFrac*(xHiRun - xLoRun));
  ptRunNeg->SetFillColor(0); ptRunNeg->SetBorderSize(1);
  ptRunNeg->AddText("-------- Fit Params --------");
  ptRunNeg->AddText(Form("p0: %.4f +/- %.4f", fRunNeg->GetParameter(0), fRunNeg->GetParError(0)));
  ptRunNeg->AddText(Form("p1: %.4f +/- %.4f", fRunNeg->GetParameter(1), fRunNeg->GetParError(1)));
  ptRunNeg->AddText(Form("Chi2 / NDF: %.3f / %i", fRunNeg->GetChisquare(), fRunNeg->GetNDF()));

  TCanvas *cCorr = new TCanvas("cCorr", "Correlation Canvas", 1200, 800);
  cCorr->Divide(2, 2);
  cCorr->cd(1); hCycPos->Draw("prof"); ptCycPos->Draw("same");
  cCorr->cd(2); hRunPos->Draw("prof"); ptRunPos->Draw("same");
  cCorr->cd(3); hCycNeg->Draw("prof"); ptCycNeg->Draw("same");
  cCorr->cd(4); hRunNeg->Draw("prof"); ptRunNeg->Draw("same");
}


void globalPullPlot(Int_t prexOrCrex, Int_t snlLo=0, Int_t snlHi=1000){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *snl = (TTree *)f->Get("snl");
  TTree *run = (TTree *)f->Get("run");
  TTree *cyc = (TTree *)f->Get("cyc");

  const Int_t nSnails = (Int_t)snl->GetEntries();
  Int_t snailCycs[nSnails];
  Int_t snailRuns[nSnails];

  initTrees(snl, run, cyc);

  // Count cycles in snail range
  Int_t curRangeCycles = 0; Int_t rangeInd = 0;
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(snlSign == 0 || snlTreeNum == 175 || snlTreeNum == 176){
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

  // Count runs in snail range
  Int_t curRangeRuns = 0; rangeInd = 0;
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(snlSign == 0 || snlTreeNum == 175 || snlTreeNum == 176){
      snailRuns[i] = 0;
      continue;
    }
    if(snlTreeNum == rangeStarts[rangeInd]){
      curRangeRuns = 0;
    }
    snailRuns[i] = numRuns;
    curRangeRuns += numRuns;
    if(snlTreeNum == rangeEnds[rangeInd]){
      rangeRuns[rangeInd] = curRangeRuns;
      rangeInd++;
    }
  }
  
  for(Int_t i = 0; i < nFitRanges; i++){
    TString histName = Form("hRange_%i-%i", rangeStarts[i], rangeEnds[i]);
    TString histTitle = Form("Snails %i-%i: Cycle Pol0 Pulls", rangeStarts[i], rangeEnds[i]);
    TString pullName = Form("hPull_%i-%i", rangeStarts[i], rangeEnds[i]);
    TString pullTitle = Form("Snails %i-%i: Cycle Pol0 Pulls by Time", rangeStarts[i], rangeEnds[i]);
    TString histNameRuns = Form("hRangeRuns_%i-%i", rangeStarts[i], rangeEnds[i]);
    TString histTitleRuns = Form("Snails %i-%i: Run Pol0 Pulls", rangeStarts[i], rangeEnds[i]);
    TString pullNameRuns = Form("hPullRuns_%i-%i", rangeStarts[i], rangeEnds[i]);
    TString pullTitleRuns = Form("Snails %i-%i: Run Pol0 Pulls by Time", rangeStarts[i], rangeEnds[i]);
    rangeHists[i] = new TH1F(histName.Data(), histTitle.Data(), 100, -8.0, 8.0);
    rangeHists[i]->GetXaxis()->SetTitle("Cycle Pol0 Pull");
    rangePulls[i] = new TH1F(pullName.Data(), pullTitle.Data(), rangeCycs[i], 0, rangeCycs[i]);
    rangePulls[i]->GetYaxis()->SetTitle("Cycle Pol0 Pull");
    rangePulls[i]->SetLineColor(kGreen); 
    rangePulls[i]->SetFillColor(kGreen);
    rangeRunHists[i] = new TH1F(histNameRuns.Data(), histTitleRuns.Data(), 100, -8.0, 8.0);
    rangeRunHists[i]->GetXaxis()->SetTitle("Run Pol0 Pull");
    rangeRunPulls[i] = new TH1F(pullNameRuns.Data(), pullTitleRuns.Data(), rangeRuns[i], 0, rangeRuns[i]);
    rangeRunPulls[i]->GetYaxis()->SetTitle("Run Pol0 Pull");
    rangeRunPulls[i]->SetLineColor(kGreen); 
    rangeRunPulls[i]->SetFillColor(kGreen);
  }

  fillCycPullPlot(cyc, snlLo, snlHi);
  fillRunPullPlot(run, snlLo, snlHi);

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
  for(Int_t i = 0; i < nFitRanges; i++){
    TCanvas *c = new TCanvas(Form("cRangeRuns_%i", i), Form("Canvas for Range %i-%i (runwise)", rangeStarts[i], rangeEnds[i]), 1400, 800);
    c->Divide(1, 2);
    c->cd(1);
    rangeRunHists[i]->SetStats(220);
    rangeRunHists[i]->Draw();
    c->cd(2)->SetGridx();
    c->cd(2)->SetGridy();
    rangeRunPulls[i]->SetStats(0);
    rangeRunPulls[i]->Draw();
  }

  graphCanvas();
}
