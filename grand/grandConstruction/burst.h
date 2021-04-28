#ifndef BURST_H
#define BURST_H

#include "plot.h"
#include "vars.h"

const Int_t burstVars = 4;
const Int_t burstStates = 4;
TString burstTitles[burstVars] = {"BurstPosAcc0", "BurstNegAcc0", "BurstDiffAcc0", "BurstSummAcc0"};
TString burstLas[burstStates] = {"LasOn", "LasOff1", "LasOff2", "LasOff"};
//TString burstTitles[burstVars] = {"BurstPosAcc0LasOn",   "BurstNegAcc0LasOn",   "BurstDiffAcc0LasOn",   "BurstSumAcc0LasOn",
//                                  "BurstPosAcc0LasOff1", "BurstNegAcc0LasOff1", "BurstDiffAcc0LasOff1", "BurstSumAcc0LasOff1",
//                                  "BurstPosAcc0LasOff2", "BurstNegAcc0LasOff2", "BurstDiffAcc0LasOff2", "BurstSumAcc0LasOff2",
//                                  "BurstPosAcc0LasOff",  "BurstNegAcc0LasOff",  "BurstDiffAcc0LasOff",  "BurstSumAcc0LasOff"};
TString burstNames[burstVars] = {pos0, neg0, Form("%s - %s", pos0.Data(), neg0.Data()), Form("%s + %s", pos0.Data(), neg0.Data())};
TString burstCuts[burstStates] = {B1L1D0, B1L0D0, B1L0D0, B1L0D0};

const Int_t burstAsyms = 4;
TString burstAsymTitles[burstAsyms] = {"BurstAsym0LasOn", "BurstAsym0LasOff1", "BurstAsym0LasOff2", "BurstAsym0LasOff"};

const Int_t comboAsyms = 3;
TString comboAsymTitles[comboAsyms] = {"BurstAsym0NGC", "BurstAsym0", "BurstPol0"};

vector<vector<int>> getBurstList(int runNum){
  ifstream cycle_infile;   
  cycle_infile.open(Form("%s/bursts_%i.dat", getenv("COMPMON_BURSTS"), runNum));

  if(!cycle_infile.is_open()){
    printf("ERROR: couldn't open file!\n");
  }
  string read_str;
  vector<vector<int>> runCycles;
  while(getline(cycle_infile, read_str)){
    vector<int> limits; stringstream ss(read_str);
    for(int i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    runCycles.push_back(limits);
  }
  return runCycles;
}

vector<vector<vector<Int_t>>> sortedCycleBursts(Int_t rNum, Int_t cNum){
  vector<vector<Int_t>> allBursts = getBurstList(rNum);

  vector<vector<vector<vector<Int_t>>>> sortedBursts;
  Int_t curCycle = 0; Int_t curPeriod = 0;
  vector<vector<vector<Int_t>>> cycleBursts;
  vector<vector<Int_t>> periodBursts;
  for(Int_t i = 0; i < allBursts.size(); i++){
    vector<Int_t> burstLimits;
    if(i == 0){
      curCycle = allBursts[i][0]; curPeriod = allBursts[i][1];
    }
    if(curPeriod != allBursts[i][1]){
      cycleBursts.push_back(periodBursts);
      periodBursts.clear(); curPeriod = allBursts[i][1];
    }
    if(curCycle != allBursts[i][0]){
      sortedBursts.push_back(cycleBursts);
      cycleBursts.clear(); curCycle = allBursts[i][0];
    }

    burstLimits.push_back(allBursts[i][2]); burstLimits.push_back(allBursts[i][3]);
    periodBursts.push_back(burstLimits);
  }
  cycleBursts.push_back(periodBursts);
  sortedBursts.push_back(cycleBursts);

  return sortedBursts[cNum - 1];
}

void cycBurstPlots(TChain* quartetwise, Int_t runNum, Int_t cycNum, TFile *runOut){
  vector<vector<vector<Int_t>>> cycBursts = sortedCycleBursts(runNum, cycNum);

  TH1F *hists[burstVars][burstStates];
  TF1 *fits[burstVars][burstStates];

  Int_t per1Bursts = (Int_t)cycBursts[0].size();
  Int_t per2Bursts = (Int_t)cycBursts[1].size();
  Int_t per3Bursts = (Int_t)cycBursts[2].size();
  Int_t cycleBursts = per1Bursts + per2Bursts + per3Bursts;
  for(Int_t var = 0; var < burstVars; var++){
    for(Int_t state = 0; state < burstStates; state++){
      TString hName = Form("h%i.%i_%s%s", runNum, cycNum, burstTitles[var].Data(), burstLas[state].Data());
      TString fName = Form("f%i.%i_%s%s", runNum, cycNum, burstTitles[var].Data(), burstLas[state].Data());
      TString hTitle = Form("Cycle %i.%i: Multiplet %s", runNum, cycNum, burstTitles[var].Data());

      hists[var][state] = new TH1F(hName.Data(), hTitle.Data(), cycleBursts, 0, cycleBursts);
        
      fits[var][state] = new TF1(fName.Data(), "pol0");
    }
  }

  for(Int_t per = 0; per < cycBursts.size(); per++){
    for(Int_t burst = 0; burst < cycBursts[per].size(); burst++){
      for(Int_t msmt = 0; msmt < burstVars; msmt++){
        Int_t laserInd = 0; Int_t burstOffset = per1Bursts;
        if(per == 0){
          laserInd = 1; burstOffset = 0;
        }
        else if(per == 2){
          laserInd = 2; burstOffset = per1Bursts + per2Bursts;
        }
        TString cuts = Form("%s && firstMPSnumber>=%i && firstMPSnumber<=%i", 
                            burstCuts[laserInd].Data(), cycBursts[per][burst][0], cycBursts[per][burst][1]);
        TString hName = Form("h%i.%i_%s%s_period%i_burst%i", runNum, cycNum, burstTitles[msmt].Data(), burstLas[laserInd].Data(), per+1, burst+1);
        quartetwise->Project(hName.Data(), burstNames[msmt], cuts.Data(), "goff");
        TH1F *hData = (TH1F *)gDirectory->Get(hName.Data());
        if(per == 0){
          hists[msmt][1]->SetBinContent(burst + burstOffset + 1, hData->GetMean());
          hists[msmt][1]->SetBinError(burst + burstOffset + 1, hData->GetMeanError());
          hists[msmt][3]->SetBinContent(burst + burstOffset + 1, hData->GetMean());
          hists[msmt][3]->SetBinError(burst + burstOffset + 1, hData->GetMeanError());
        }
        else if(per == 1){
          hists[msmt][0]->SetBinContent(burst + burstOffset + 1, hData->GetMean());
          hists[msmt][0]->SetBinError(burst + burstOffset + 1, hData->GetMeanError());
        }
        else if(per == 2){
          hists[msmt][2]->SetBinContent(burst + burstOffset + 1, hData->GetMean());
          hists[msmt][2]->SetBinError(burst + burstOffset + 1, hData->GetMeanError());
          hists[msmt][3]->SetBinContent(burst + burstOffset + 1, hData->GetMean());
          hists[msmt][3]->SetBinError(burst + burstOffset + 1, hData->GetMeanError());
        }
      }
    }
  }

  for(Int_t burst = 0; burst < burstVars; burst++){
    for(Int_t state = 0; state < burstStates; state++){
      hists[burst][state]->Fit(fits[burst][state], "Q", "goff", 0, cycleBursts);
    }
  }

  Float_t sumOnMean = fits[3][0]->GetParameter(0);
  Float_t sumOffMean = fits[3][3]->GetParameter(0);
  TString asymLasOn = Form("(%s - %s)/(%s + %s - %f)", pos0.Data(), neg0.Data(), pos0.Data(), neg0.Data(), sumOffMean);
  TString asymLasOff = Form("(%s - %s)/(%f - %f)", pos0.Data(), neg0.Data(), sumOnMean, sumOffMean);
  TH1F *asymHists[burstStates];
  for(Int_t state = 0; state < burstStates; state++){
    TString hAsymName = Form("h%i.%i_%s", runNum, cycNum, burstAsymTitles[state].Data());
    TString hAsymTitle = Form("Cycle %i.%i: Multiplet Asym0", runNum, cycNum);

    asymHists[state] = new TH1F(hAsymName.Data(), hAsymTitle.Data(), cycleBursts, 0, cycleBursts);
  }

  for(Int_t per = 0; per < cycBursts.size(); per++){
    for(Int_t burst = 0; burst < cycBursts[per].size(); burst++){
      Int_t laserInd = 0; TString asymMsmt = asymLasOn;
      Int_t burstOffset = per1Bursts;
      if(per==0){laserInd = 1; asymMsmt = asymLasOff; burstOffset = 0;}
      else if(per==2){laserInd = 2; asymMsmt = asymLasOff; burstOffset = per1Bursts + per2Bursts;}
      TString cuts = Form("%s && firstMPSnumber>=%i && firstMPSnumber<=%i", 
                          burstCuts[laserInd].Data(), cycBursts[per][burst][0], cycBursts[per][burst][1]);
      TString hAsymName = Form("h%i.%i_Asym0%s_period%i_burst%i", runNum, cycNum, burstLas[laserInd].Data(), per+1, burst+1);
      quartetwise->Project(hAsymName.Data(), asymMsmt.Data(), cuts.Data(), "goff");
      TH1F *hAsymData = (TH1F *)gDirectory->Get(hAsymName.Data());
      if(per==0){
        asymHists[1]->SetBinContent(burst + burstOffset + 1, hAsymData->GetMean());
        asymHists[1]->SetBinError(burst + burstOffset + 1, hAsymData->GetMeanError());
        asymHists[3]->SetBinContent(burst + burstOffset + 1, hAsymData->GetMean());
        asymHists[3]->SetBinError(burst + burstOffset + 1, hAsymData->GetMeanError());
      }
      else if(per==1){
        asymHists[0]->SetBinContent(burst + burstOffset + 1, hAsymData->GetMean());
        asymHists[0]->SetBinError(burst + burstOffset + 1, hAsymData->GetMeanError());
      }
      else if(per==2){
        asymHists[2]->SetBinContent(burst + burstOffset + 1, hAsymData->GetMean());
        asymHists[2]->SetBinError(burst + burstOffset + 1, hAsymData->GetMeanError());
        asymHists[3]->SetBinContent(burst + burstOffset + 1, hAsymData->GetMean());
        asymHists[3]->SetBinError(burst + burstOffset + 1, hAsymData->GetMeanError());
      }
    }
  }

  for(Int_t var = 0; var < burstVars; var++){
    for(Int_t state = 0; state < burstStates; state++){
      runOut->cd();
      hists[var][state]->Write();
    }
  }

  for(Int_t state = 0; state < burstStates; state++){
    runOut->cd();
    asymHists[state]->Write();
  }
}


#endif
