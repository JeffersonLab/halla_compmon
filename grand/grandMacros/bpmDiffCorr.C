#include "../grandOnline/makePlots.h"

using namespace std;

vector<vector<Float_t>> readCorrsFile(TString expt){
  vector<vector<Float_t>> corrs;
  ifstream keysfile(Form("%s/%s_corrSlopes.key", getenv("COMPMON_MAPS"), expt.Data()));
  if(keysfile.is_open()){
    string line;
    while(getline(keysfile, line)){
      vector<Float_t> keyRange; stringstream ss(line);
      for(Float_t i; ss >> i;){
        keyRange.push_back(i);
        if(ss.peek() == ',')
          ss.ignore();
      }
      corrs.push_back(keyRange);
    }
  }
  keysfile.close();
  return corrs;
}

void bpmDiffCorr(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  vector<vector<Float_t>> corrs = readCorrsFile(expt);
  
  DataVar diff_bpmAy, diff_bpmBy;
  Int_t snailNum, runNum, cycNum, sign, cycleCut;

  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("diff_bpmAy", &diff_bpmAy);
  cyc->SetBranchAddress("diff_bpmBy", &diff_bpmBy);
  cyc->SetBranchAddress("CycleCut", &cycleCut);

  TH1F *hDiffCorr = new TH1F("hDiffCorr", "Size of Sign-Corrected BPM Y-Difference Correction", 200, -5, 5);
  hDiffCorr->GetXaxis()->SetTitle("Asymmetry Correction [ppt]");

  Int_t rangeInd = 0;
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(runNum > corrs[rangeInd][1]){rangeInd++;}
    if(rangeInd >= corrs.size()){break;}
    if(runNum < corrs[rangeInd][0] || runNum > corrs[rangeInd][1] || cycleCut!=0){continue;}

    Float_t corr = 1000.0*corrs[rangeInd][2]*sign*(diff_bpmAy.mean + diff_bpmBy.mean)/2.0;
    hDiffCorr->Fill(corr);
  }

  TCanvas *c = new TCanvas("c", "Diff Correction Canvas", 1200, 800);
  c->cd();
  hDiffCorr->Draw();
}
