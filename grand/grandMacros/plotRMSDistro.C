#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "../vars.h"

using namespace std;

void writeRootfile(vector<vector<vector<Int_t>>> cutInfo, Float_t facStart, Float_t facStep){
  TFile *outfile = new TFile("./cutCycles.root", "RECREATE");
  
  for(Int_t conf = 0; conf < cutInfo.size(); conf++){
    Int_t runNum, cycNum, flag;
    TTree *cut = new TTree(Form("cut_%.3f", facStart + facStep*conf), "a tree for cut cycles and their info");
    cut->Branch("runNum", &runNum, "runNum/I");
    cut->Branch("cycNum", &cycNum, "cycNum/I");
    cut->Branch("flag", &flag, "flag/I");

    for(Int_t cyc = 0; cyc < cutInfo[conf][0].size(); cyc++){
      runNum = cutInfo[conf][0][cyc];
      cycNum = cutInfo[conf][1][cyc];
      flag = cutInfo[conf][2][cyc];
      cut->Fill();
    }
    
    outfile->cd();
    cut->Write();
  }

  outfile->Close();
}

void printFactorCuts(Float_t factor){
  TFile *infile = new TFile("./cutCycles.root", "READ");
  TTree *cut = (TTree *)infile->Get(Form("cut_%.3f", factor));

  Int_t runNum, cycNum, flag;
  cut->SetBranchAddress("runNum", &runNum);
  cut->SetBranchAddress("cycNum", &cycNum);
  cut->SetBranchAddress("flag", &flag);

  for(Int_t i = 0; i < cut->GetEntries(); i++){
    cut->GetEntry(i);
    printf("%04i | Cycle %04i.%02i cut with flag %i\n", i+1, runNum, cycNum, flag);
  }
  infile->Close();
}

void printFalsePositiveRates(Float_t facStart, Float_t facEnd, Float_t facStep){
  TFile *infile = new TFile("./cutCycles.root", "READ");

  TCanvas *cRates = new TCanvas("cRates", "False Positive Rates", 1200, 800);
  cRates->SetGridx(); cRates->SetGridy();
  TGraph *g = new TGraph();
  Int_t nPts = 0;
  for(Float_t factor = facStart; factor <= facEnd; factor += facStep){
    TTree *cut = (TTree *)infile->Get(Form("cut_%.3f", factor));
    Float_t rate = cut->GetEntries("flag == 1 || flag ==2")*1.0/cut->GetEntries();
    g->SetPoint(nPts++, factor, rate);
  }

  cRates->cd();
  g->SetTitle("\"False Positive\" Rate by Factor");
  g->GetXaxis()->SetTitle("RMS Factor");
  g->GetYaxis()->SetTitle("#\"False Positive\" Rate");
  g->SetMarkerStyle(8);
  g->Draw("ap");
}

void cutDiffs(Float_t facLo, Float_t facHi){
  if(facHi <= facLo){
    printf("ERROR: Factors for comparison are in the wrong order.\n");
    return;
  }

  TFile *infile = new TFile("./cutCycles.root", "READ");
  TTree *cutLo = (TTree *)infile->Get(Form("cut_%.3f", facLo));
  TTree *cutHi = (TTree *)infile->Get(Form("cut_%.3f", facHi));

  Int_t runNumLo, cycNumLo, runNumHi, cycNumHi;

  cutLo->SetBranchAddress("runNum", &runNumLo);
  cutLo->SetBranchAddress("cycNum", &cycNumLo);
  cutHi->SetBranchAddress("runNum", &runNumHi);
  cutHi->SetBranchAddress("cycNum", &cycNumHi);

  Int_t lo = 0; Int_t hi = 0;
  while(lo < cutLo->GetEntries() && hi < cutHi->GetEntries()){
    cutLo->GetEntry(lo);
    cutHi->GetEntry(hi);

    if(runNumLo == runNumHi && cycNumLo == cycNumHi){
      lo++; hi++;
    }
    else{
      printf("Cycle %04i.%02i cut in %.3f but not in %.3f\n", runNumLo, cycNumLo, facLo, facHi);
      lo++;
    }
  }
}

void plotRMSDistro(Int_t prexOrCrex, Float_t specFactor, Float_t specFactor2){
  TString exptStr("");
  TString exptName("");
  if(prexOrCrex == 1){
    exptStr = "prex";
    exptName = "PREX-II";
  }
  else if(prexOrCrex == 2){
    exptStr = "crex";
    exptName = "CREX";
  }
  else{
    printf("Invalid Experiment Code\n");
    exit(1);
  }
  ofstream keys; keys.open(Form("%s/prex_rmsCut.key", getenv("COMPMON_MAPS")));
  
  TFile *f = TFile::Open(Form("%s/%sGrandCompton.root", getenv("COMPMON_GRAND"), exptStr.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  const Int_t nGroups = 14;
  Int_t runStarts[nGroups] = {4301, 4310, 4342, 4344, 4353, 4371, 4405,
                              4416, 4472, 4487, 4509, 4551, 4577, 4606};
  Int_t runStops[nGroups] =  {4309, 4341, 4343, 4352, 4370, 4404, 4415,
                              4471, 4486, 4508, 4550, 4576, 4605, 4622};
  Int_t groupCuts[nGroups] = {0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0};
  Int_t groupCycs[nGroups] = {0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0, 0};
  vector<vector<Float_t>> rmsOffLimits;
  vector<vector<Float_t>> rmsOnLimits;
  
  Float_t facStart = 0.02;
  Float_t facEnd = 0.2;
  Float_t facStep = 0.002;
  
  Int_t j = 0;
  for(Float_t factor = facStart; atof(Form("%.3f", factor)) <= facEnd; factor += facStep){
    vector<Float_t> rmsOffLimitsFac;
    vector<Float_t> rmsOnLimitsFac;

    for(Int_t i = 0; i < nGroups; i++){
      //TCanvas *c = new TCanvas(Form("cRMS_%i", i+1), Form("RMS Distribution %i", i+1), 1200, 600);
      //c->Divide(2, 1);
      TString hNameOff = Form("hOff_%i_%i", j, i+1);
      TString hNameOff1 = Form("hOff1_%i_%i", j, i+1);
      TString hNameOff2 = Form("hOff2_%i_%i", j, i+1);
      TString hNameOn = Form("hOn_%i_%i", j++, i+1);

      Int_t nBins = 200; Float_t xmin = 0.0; Float_t xmax = 0.7;
      TH1F *hOff = new TH1F(hNameOff.Data(), "Laser Off RMS's", nBins, xmin, xmax);
      TH1F *hOff1 = new TH1F(hNameOff1.Data(), "Laser Off RMS's", nBins, xmin, xmax);
      TH1F *hOff2 = new TH1F(hNameOff2.Data(), "Laser Off RMS's", nBins, xmin, xmax);
      TH1F *hOn = new TH1F(hNameOn.Data(), "Laser On RMS's", nBins, xmin, xmax);

      //c->cd(1); cyc->Draw(Form("Acc0LasOn.rms>>%s", hNameOff1.Data()), Form("runNum>=%i && runNum<=%i", runStarts[i], runStops[i]));
      //c->cd(2); cyc->Draw(Form("Acc0LasOff2.rms>>%s", hNameOff2.Data()), Form("runNum>=%i && runNum<=%i", runStarts[i], runStops[i]));
      cyc->Project(hNameOff1.Data(), "Acc0LasOff1.rms", Form("runNum>=%i && runNum<=%i", runStarts[i], runStops[i]));
      cyc->Project(hNameOff2.Data(), "Acc0LasOff2.rms", Form("runNum>=%i && runNum<=%i", runStarts[i], runStops[i]));
      cyc->Project(hNameOn.Data(), "Acc0LasOn.rms", Form("runNum>=%i && runNum<=%i", runStarts[i], runStops[i]));
    
      hOff->Add(hOff1, hOff2);

      //rmsOffLimitsFac.push_back(hOff->GetMean() + factor*hOff->GetRMS());  
      //rmsOnLimitsFac.push_back(hOn->GetMean() + factor*hOn->GetRMS());
      //rmsOffLimitsFac.push_back(factor*hOff->GetMean());  
      //rmsOnLimitsFac.push_back(factor*hOn->GetMean());
      if(atof(Form("%.3f", factor)) == 0.15){
        printf("For factor 0.05 and run range %i-%i off limit is: %.4f because mode is %.4f\n", 
                runStarts[i], runStops[i], hOff->GetBinCenter(hOff->GetMaximumBin()) + factor, hOff->GetBinCenter(hOff->GetMaximumBin()));
        keys<<Form("%i,%i,%.4f,%.4f,%.4f,%i,%.2f,%.2f\n", runStarts[i], runStops[i], 
                    hOff->GetBinCenter(hOff->GetMaximumBin()), hOn->GetBinCenter(hOn->GetMaximumBin()), factor, nBins, xmin, xmax);
      }
      rmsOffLimitsFac.push_back(hOff->GetBinCenter(hOff->GetMaximumBin()) + factor);
      rmsOnLimitsFac.push_back(hOn->GetBinCenter(hOn->GetMaximumBin()) + factor);
      //printf("Histogram Mode: %f\n", hOff->GetBinCenter(hOff->GetMaximumBin())); 
    }
    rmsOffLimits.push_back(rmsOffLimitsFac);
    rmsOnLimits.push_back(rmsOnLimitsFac);
  }
  keys.close();

  Int_t runNum, cycNum;
  DataVar acc0Off1, acc0Off2, acc0On;

  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("Acc0LasOn", &acc0On);
  cyc->SetBranchAddress("Acc0LasOff1", &acc0Off1);
  cyc->SetBranchAddress("Acc0LasOff2", &acc0Off2);

  TGraph *g = new TGraph();
  TString hTitle("RMS Cut Comparison (Las Off)"); Int_t nBins = 50;
  TH1F *hRMS_full = new TH1F("hRMS_full", hTitle.Data(), nBins, 0.0, 1.2);
  TH1F *hRMS_1 = new TH1F("hRMS_1", hTitle.Data(), nBins, 0.0, 1.2);
  TH1F *hRMS_2 = new TH1F("hRMS_2", hTitle.Data(), nBins, 0.0, 1.2);
  TH1F *hRMS_3 = new TH1F("hRMS_3", hTitle.Data(), nBins, 0.0, 1.2);
  TH1F *hRMS_4 = new TH1F("hRMS_4", hTitle.Data(), nBins, 0.0, 1.2);
  TH1F *hRMS_5 = new TH1F("hRMS_5", hTitle.Data(), nBins, 0.0, 1.2);
  TH1F *hRMS_spec = new TH1F("hRMS_spec", hTitle.Data(), nBins, 0.0, 1.2);

  vector<vector<vector<Int_t>>> allCutConfs;
  vector<vector<Int_t>> cutConf;
  vector<Int_t> runNums; vector<Int_t> cycNums; vector<Int_t> flags;
  Float_t selFacs[5] = {0.04, 0.08, 0.12, 0.16, 0.2};
  for(Float_t factor = facStart; atof(Form("%.3f", factor)) <= facEnd; factor += facStep){
    cutConf.clear(); runNums.clear(); cycNums.clear(); flags.clear();
    Int_t j = (Int_t)((factor - facStart)/facStep);
    for(Int_t i = 0; i < cyc->GetEntries(); i++){
      cyc->GetEntry(i);

      Int_t runRange = 0;
      while(!(runNum >= runStarts[runRange] && runNum <= runStops[runRange])){
        runRange++;
      }

      Int_t flagOff1 = (Int_t)(acc0Off1.rms > rmsOffLimits[j][runRange]);
      Int_t flagOff2 = (Int_t)(acc0Off2.rms > rmsOffLimits[j][runRange]);
      Int_t flagOn = (Int_t)(acc0On.rms > rmsOnLimits[j][runRange]);
      Int_t fullFlag = 0x1*flagOff1 + 0x2*flagOff2 + 0x4*flagOn;

      if(factor == facStart){hRMS_full->Fill(acc0Off1.rms); hRMS_full->Fill(acc0Off2.rms);}
      if(atof(Form("%.3f", factor)) == atof(Form("%.3f", specFactor))){
        groupCycs[runRange]++;
      }
      if(fullFlag != 0){
        runNums.push_back(runNum);
        cycNums.push_back(cycNum);
        flags.push_back(fullFlag);
        if(atof(Form("%.3f", factor)) == atof(Form("%.3f", specFactor))){
          groupCuts[runRange]++;
        }
      }
      else if(fullFlag == 0 && atof(Form("%.3f", factor)) == atof(Form("%.3f", selFacs[0]))){hRMS_1->Fill(acc0Off1.rms); hRMS_1->Fill(acc0Off2.rms);}
      else if(fullFlag == 0 && atof(Form("%.3f", factor)) == atof(Form("%.3f", selFacs[1]))){hRMS_2->Fill(acc0Off1.rms); hRMS_2->Fill(acc0Off2.rms);}
      else if(fullFlag == 0 && atof(Form("%.3f", factor)) == atof(Form("%.3f", selFacs[2]))){hRMS_3->Fill(acc0Off1.rms); hRMS_3->Fill(acc0Off2.rms);}
      else if(fullFlag == 0 && atof(Form("%.3f", factor)) == atof(Form("%.3f", selFacs[3]))){hRMS_4->Fill(acc0Off1.rms); hRMS_4->Fill(acc0Off2.rms);}
      else if(fullFlag == 0 && atof(Form("%.3f", factor)) == atof(Form("%.3f", selFacs[4]))){hRMS_5->Fill(acc0Off1.rms); hRMS_5->Fill(acc0Off2.rms);}
      else if(fullFlag == 0 && atof(Form("%.3f", factor)) == atof(Form("%.3f", specFactor))){hRMS_spec->Fill(acc0Off1.rms); hRMS_spec->Fill(acc0Off2.rms);}
    }
    cutConf.push_back(runNums);
    cutConf.push_back(cycNums);
    cutConf.push_back(flags);
    allCutConfs.push_back(cutConf);
    printf("RMS Factor %.3f cut %i cycles\n", factor, (Int_t)cutConf[0].size());
    g->SetPoint(j, factor, (Float_t)cutConf[0].size());
  }

  TCanvas *c = new TCanvas("c", "Cut Graph", 1200, 800);
  c->SetGridx(); c->SetGridy();
  g->SetTitle("# Cycles Cut by Mean Factor");
  g->GetXaxis()->SetTitle("RMS Factor");
  g->GetYaxis()->SetTitle("# Cycles Cut");
  g->SetMarkerStyle(8);
  g->Draw("ap");

  writeRootfile(allCutConfs, facStart, facStep);
  printFactorCuts(specFactor2);
  //printFalsePositiveRates(facStart, facEnd, facStep);
  cutDiffs(specFactor, specFactor2);

  TCanvas *cComp = new TCanvas("cComp", "RMS Cut Comparison", 1200, 800);
  THStack *hs = new THStack("hs", hTitle.Data());
  hRMS_full->SetLineColor(kBlack); hs->Add(hRMS_full);
  hRMS_5->SetLineColor(kRed); hs->Add(hRMS_5);
  hRMS_4->SetLineColor(kOrange); hs->Add(hRMS_4);
  hRMS_3->SetLineColor(kGreen); hs->Add(hRMS_3);
  hRMS_2->SetLineColor(kBlue); hs->Add(hRMS_2);
  hRMS_1->SetLineColor(kViolet-1); hs->Add(hRMS_1);
  hRMS_spec->SetLineColor(kMagenta); hs->Add(hRMS_spec);
  TLegend *leg = new TLegend(0.75, 0.75, 0.95, 0.9, "", "NDC");
  leg->AddEntry(hRMS_full, "No Cuts"); leg->AddEntry(hRMS_5, Form("Factor=%.3f", selFacs[4]));
  leg->AddEntry(hRMS_4, Form("Factor=%.3f", selFacs[3])); leg->AddEntry(hRMS_3, Form("Factor=%.3f", selFacs[2]));
  leg->AddEntry(hRMS_2, Form("Factor=%.3f", selFacs[1])); leg->AddEntry(hRMS_1, Form("Factor=%.3f", selFacs[0]));
  leg->AddEntry(hRMS_spec, Form("Factor=%.3f", specFactor));
  cComp->cd(); hs->Draw("nostack"); leg->Draw("same");
  
  for(Int_t i = 0; i < nGroups; i++){
    printf("Runs %i-%i: %i cycles cut out of %i (%.3f%%)\n", runStarts[i], runStops[i], groupCuts[i], groupCycs[i], groupCuts[i]*100.0/groupCycs[i]);
  }
}
