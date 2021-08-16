#include "aggregate.h"
#include "../../online/runs.h"

using namespace std;

vector<Int_t> getSnailNums(Int_t grpNum){
  vector<int> snailList;

  printf("Reading snails from file esc%i.list...\n", grpNum);
  string runNumStr;
  ifstream infile(Form("%s/escargatoires/esc%i.list", getenv("COMPMON_SNAILS"), grpNum));
  if(infile.good()){
    while(getline(infile, runNumStr)){
      Int_t snailNum = atoi(runNumStr.c_str());
      if(snailNum == 0) continue;
      snailList.push_back(snailNum);
    }
  }
  else{
    printf("Could not find escargatoire file!\n");
  }
  return snailList;
}

vector<Int_t> snailBoundaries(vector<Int_t> snailList, vector<vector<Int_t>> runList){
  vector<Int_t> bounds;
  bounds.push_back(0);
  for(Int_t i = 1; i < snailList.size(); i++){
    Int_t snlNum = snailList[i];
    Int_t firstRun = runList[snlNum - 101][0];
    bounds.push_back(firstRun);
  }
  return bounds;
}

void aggregateSnailGroup(Int_t grpNum){
  TFile *grand = new TFile(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), exptName(grpNum).Data()), "READ");
  TTree *cyc = (TTree *)grand->Get("cyc");
  vector<Int_t> snailList = getSnailNums(grpNum);
  vector<vector<Int_t>> runList = productionRunList(2);
  vector<Int_t> snailBounds = snailBoundaries(snailList, runList);
  
  TString snlCut = "(";
  for(Int_t i = 0; i < snailList.size(); i++){
    if(i == snailList.size() - 1){
      snlCut += Form("snailNum==%i)", snailList[i]);
    }
    else{
      snlCut += Form("snailNum==%i || ", snailList[i]);
    }
  }
  printf("Snail Cut: %s\n", snlCut.Data());

  nCycles = (Int_t)cyc->GetEntries(snlCut.Data());
  nCyclesCut = (Int_t)cyc->GetEntries(Form("%s && CycleCut==0", snlCut.Data()));
  printf("Number of cycles: %i\n", nCycles);
  Int_t boundInd = 1;
  vector<Float_t> bounds;
  vector<Float_t> cutBounds;

  initHists(grpNum);
  initFits();
  initBranches(cyc);
  Int_t labelFreq = 5;
  
  Int_t nBin = 1; Int_t nBinCut = 1;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    Bool_t skip = kTRUE;
    for(Int_t j = 0; j < snailList.size(); j++){
      if(snailList[j] == snailNum) skip = kFALSE;
    }
    if(skip) continue;
    
    if(boundInd < snailBounds.size() && runNum >= snailBounds[boundInd]){
      bounds.push_back(1.0*(nBin - 1));
      cutBounds.push_back(1.0*(nBinCut - 1));
      boundInd++;
    }

    for(Int_t j = 0; j < 2*nPols; j += 2){
      Float_t binVal = pols0[(Int_t)(j/2)].mean*pol0Mults[(Int_t)(j/2)];
      Float_t binErr = pols0[(Int_t)(j/2)].meanErr*pol0Mults[(Int_t)(j/2)];
      pol0Hists[j+1]->SetBinContent(nBin, binVal);
      pol0Hists[j+1]->SetBinError(nBin, binErr);
      pol0Min[j+1] = (binVal - binErr < pol0Min[j+1]) ? binVal - binErr : pol0Min[j+1];
      pol0Max[j+1] = (binVal + binErr > pol0Max[j+1]) ? binVal + binErr : pol0Max[j+1];
      if(nBin % labelFreq == 1){
        pol0Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
      if(cycCut == 0){
        pol0Hists[j]->SetBinContent(nBinCut, binVal);
        pol0Hists[j]->SetBinError(nBinCut, binErr);
        pol0Min[j] = (binVal - binErr < pol0Min[j]) ? binVal - binErr : pol0Min[j];
        pol0Max[j] = (binVal + binErr > pol0Max[j]) ? binVal + binErr : pol0Max[j];
      }
      if(nBinCut % labelFreq == 1 && cycCut == 0){
        pol0Hists[j]->GetXaxis()->SetBinLabel(nBinCut, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nAccVars; j+=2){
      var0Hists[j]->SetBinContent(nBin, vars0[(Int_t)j/2].mean);
      var0Hists[j]->SetBinError(nBin, vars0[(Int_t)j/2].meanErr);
      var0Min[j] = (vars0[(Int_t)j/2].mean - vars0[(Int_t)j/2].meanErr < var0Min[j]) ? vars0[(Int_t)j/2].mean - vars0[(Int_t)j/2].meanErr : var0Min[j];
      var0Max[j] = (vars0[(Int_t)j/2].mean + vars0[(Int_t)j/2].meanErr > var0Max[j]) ? vars0[(Int_t)j/2].mean + vars0[(Int_t)j/2].meanErr : var0Max[j];
      var0Hists[j+1]->SetBinContent(nBin, vars0[(Int_t)j/2].rms);
      var0Hists[j+1]->SetBinError(nBin, vars0[(Int_t)j/2].rmsErr);
      var0Min[j+1] = (vars0[(Int_t)j/2].rms - vars0[(Int_t)j/2].rmsErr < var0Min[j+1]) ? vars0[(Int_t)j/2].rms - vars0[(Int_t)j/2].rmsErr : var0Min[j+1];
      var0Max[j+1] = (vars0[(Int_t)j/2].rms + vars0[(Int_t)j/2].rmsErr > var0Max[j+1]) ? vars0[(Int_t)j/2].rms + vars0[(Int_t)j/2].rmsErr : var0Max[j+1];
      if(nBin % labelFreq == 1){
        var0Hists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        var0Hists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nVars; j+=2){
      varHists[j]->SetBinContent(nBin, vars[(Int_t)j/2].mean);
      varHists[j]->SetBinError(nBin, vars[(Int_t)j/2].meanErr);
      varMin[j] = (vars[(Int_t)j/2].mean - vars[(Int_t)j/2].meanErr < varMin[j]) ? vars[(Int_t)j/2].mean - vars[(Int_t)j/2].meanErr : varMin[j];
      varMax[j] = (vars[(Int_t)j/2].mean + vars[(Int_t)j/2].meanErr > varMax[j]) ? vars[(Int_t)j/2].mean + vars[(Int_t)j/2].meanErr : varMax[j];
      varHists[j+1]->SetBinContent(nBin, vars[(Int_t)j/2].rms);
      varHists[j+1]->SetBinError(nBin, vars[(Int_t)j/2].rmsErr);
      varMin[j+1] = (vars[(Int_t)j/2].rms - vars[(Int_t)j/2].rmsErr < varMin[j+1]) ? vars[(Int_t)j/2].rms - vars[(Int_t)j/2].rmsErr : varMin[j+1];
      varMax[j+1] = (vars[(Int_t)j/2].rms + vars[(Int_t)j/2].rmsErr > varMax[j+1]) ? vars[(Int_t)j/2].rms + vars[(Int_t)j/2].rmsErr : varMax[j+1];
      if(nBin % labelFreq == 1){
        varHists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        varHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nBursts; j+=2){
      burstHists[j]->SetBinContent(nBin, burstVars[(Int_t)j/2].mean);
      burstHists[j]->SetBinError(nBin, burstVars[(Int_t)j/2].meanErr);
      burstMin[j] = (burstVars[(Int_t)j/2].mean - burstVars[(Int_t)j/2].meanErr < burstMin[j]) ? burstVars[(Int_t)j/2].mean - burstVars[(Int_t)j/2].meanErr : burstMin[j];
      burstMax[j] = (burstVars[(Int_t)j/2].mean + burstVars[(Int_t)j/2].meanErr > burstMax[j]) ? burstVars[(Int_t)j/2].mean + burstVars[(Int_t)j/2].meanErr : burstMax[j];
      if(burstVars[(Int_t)j/2].NDF == 0){burstHists[j+1]->SetBinContent(nBin, 0.0);}
      else{
        burstHists[j+1]->SetBinContent(nBin, burstVars[(Int_t)j/2].Chi2*1.0/burstVars[(Int_t)j/2].NDF);
        burstMin[j+1] = (burstVars[(Int_t)j/2].Chi2*1.0/burstVars[(Int_t)j/2].NDF < burstMin[j+1]) ? burstVars[(Int_t)j/2].Chi2*1.0/burstVars[(Int_t)j/2].NDF : burstMin[j+1];
        burstMax[j+1] = (burstVars[(Int_t)j/2].Chi2*1.0/burstVars[(Int_t)j/2].NDF > burstMax[j+1]) ? burstVars[(Int_t)j/2].Chi2*1.0/burstVars[(Int_t)j/2].NDF : burstMax[j+1];
      }
      if(nBin % labelFreq == 1){
        burstHists[j]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
        burstHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
    }
    for(Int_t j = 0; j < 2*nCombos; j+=2){
      Float_t binVal = comboVars[(Int_t)(j/2)].mean*comboMults[(Int_t)(j/2)];
      Float_t binErr = comboVars[(Int_t)(j/2)].meanErr*comboMults[(Int_t)(j/2)];
      comboHists[j+1]->SetBinContent(nBin, binVal);
      comboHists[j+1]->SetBinError(nBin, binErr);
      comboMin[j+1] = (binVal - binErr < comboMin[j+1]) ? binVal - binErr : comboMin[j+1];
      comboMax[j+1] = (binVal + binErr > comboMax[j+1]) ? binVal + binErr : comboMax[j+1];
      if(nBin % labelFreq == 1){
        comboHists[j+1]->GetXaxis()->SetBinLabel(nBin, Form("%i.%i", runNum, cycleNum));
      }
      if(cycCut == 0){
        comboHists[j]->SetBinContent(nBinCut, binVal);
        comboHists[j]->SetBinError(nBinCut, binErr);
        comboMin[j] = (binVal - binErr < comboMin[j]) ? binVal - binErr : comboMin[j];
        comboMax[j] = (binVal + binErr > comboMax[j]) ? binVal + binErr : comboMax[j];
      }
      if(nBinCut % labelFreq == 1 && cycCut == 0){
        comboHists[j]->GetXaxis()->SetBinLabel(nBinCut, Form("%i.%i", runNum, cycleNum));
      }
    }
    nBin++;
    if(cycCut == 0){nBinCut++;}
  }

  Int_t msmt = 0;
  //gStyle->SetOptStat(0);
  makePol0Lines(grpNum, bounds, cutBounds);
  makeVar0Lines(grpNum, bounds);
  makeVarLines(grpNum, bounds);
  makeBurstLines(grpNum, bounds);
  makeComboLines(grpNum, bounds, cutBounds);
  makePol0Plots(grpNum, &msmt);
  makeVar0Plots(grpNum, &msmt);
  makeVarPlots(grpNum,&msmt);
  makeBurstPlots(grpNum, &msmt);
  makeComboPlots(grpNum, &msmt);

  TString outputDir = Form("%s/snails/snail%i", getenv("COMPMON_WEB"), grpNum);
  gSystem->Exec(Form("pdfunite %s/agg_plots_*.pdf %s/snail%i_agg_plots.pdf", outputDir.Data(), outputDir.Data(), grpNum));
  gSystem->Exec(Form("rm -f %s/agg_plots_*.pdf", outputDir.Data()));
}
