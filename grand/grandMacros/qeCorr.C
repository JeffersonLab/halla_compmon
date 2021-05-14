#include "../grandOnline/makePlots.h"

using namespace std;


vector<vector<Float_t>> readQEFile(){
  vector<vector<Float_t>> snlQE;
  //vector<Float_t> bpmMeans, bpmErrs;
  ifstream infile("qe.csv");
  string readStr;
  while(getline(infile, readStr)){
    stringstream ss(readStr);
    string token;
    vector<Float_t> qeLine;
    while(getline(ss, token, ',')){
      qeLine.push_back(atof(token.c_str()));
    }
    snlQE.push_back(qeLine);
  }

  return snlQE;
}


vector<Int_t> formSnailList(vector<vector<Float_t>> fullList){
  vector<Int_t> snailList;
  for(Int_t i = 0; i < fullList.size(); i++){
    snailList.push_back((Int_t)fullList[i][0]);
  }
  return snailList;
}


vector<Float_t> formQEList(vector<vector<Float_t>> fullList){
  vector<Float_t> qeList;
  for(Int_t i = 0; i < fullList.size(); i++){
    qeList.push_back(fullList[i][1]);
  }
  return qeList;
}


void qeCorr(){
  TFile *f = TFile::Open(Form("%s/aggregates/crexGrandCompton.root", getenv("COMPMON_WEB")));
  TTree *snl = (TTree *)f->Get("snl");

  FitPolVar pol0;
  Int_t snlNum, sign;
  Float_t ihwp;
  const Int_t rangeUpper = 137;
  vector<vector<Float_t>> fullList = readQEFile();
  vector<Int_t> snailList = formSnailList(fullList);
  vector<Float_t> qeList = formQEList(fullList);
  // for(Int_t i = 0; i < fullList.size(); i++){
  //   printf("Snail %i with QE %f\n", snailList[i], qeList[i]);
  // }
  // const Int_t rangeMsmts1 = 6;
  // Float_t qe1[rangeMsmts1] = {0.2296, 0.2656, 0.2913, 0.3468, 0.4171, 0.4628};
  // Float_t pol0Mean1[rangeMsmts1] = {86.2742, 86.73, 86.4, 86.8992, 87.86, 87.686};
  // Float_t pol0Err1[rangeMsmts1] = {0.0838, 0.217, 0.157, 0.0948, 0.7189, 0.1566};
  
  // const Int_t rangeMsmts2 = 13;
  // Float_t qe2[rangeMsmts2] = {0.5527, 0.5942, 0.6257, 0.601, 0.5917, 0.5807, 0.552, 0.5496, 0.5514, 0.5389, 0.481, 0.4479, 0.4335};
  // Float_t pol0Mean2[rangeMsmts2] = {87.73, 87.45, 87.7354, 87.685, 87.82, 87.725, 87.26, 87.58, 87.37, 87.7, 87.4042, 87.1, 87.1363};
  // Float_t pol0Err2[rangeMsmts2] = {0.1235, 0.1147, 0.0855, 0.0815, 0.1208, 0.0905, 0.1275, 0.1132, 0.1304, 0.1319, 0.1067, 0.1029, 0.0842};

  snl->SetBranchAddress("snailNum", &snlNum);
  snl->SetBranchAddress("sign", &sign);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("Pol0", &pol0);

  TCanvas *c1 = new TCanvas("c1", "Corr Before Spot Move", 1200, 800);
  TCanvas *c2 = new TCanvas("c2", "Corr After Spot Move", 1200, 800);

  TGraphErrors *g1 = new TGraphErrors();
  TGraphErrors *g2 = new TGraphErrors();
  g1->SetTitle("QE vs Pol0 Correlation Before Spot Move");
  g1->GetXaxis()->SetTitle("qe");
  // g1->GetYaxis()->SetTitle("Asym0 [ppt]");
  g1->GetYaxis()->SetTitle("Pol0 [pct]");
  g1->SetMarkerStyle(20);
  g2->SetTitle("QE vs Pol0 Correlation After Spot Move");
  g2->GetXaxis()->SetTitle("qe");
  // g2->GetYaxis()->SetTitle("Asym0 [ppt]");
  g2->GetYaxis()->SetTitle("Pol0 [pct]");
  g2->SetMarkerStyle(20);
  TF1 *fit1 = new TF1("f1", "pol1");
  TF1 *fit2 = new TF1("f2", "pol1");
  TPaveText *pt1 = new TPaveText(0.125, 0.75, 0.3, 0.9, "blNDC");
  TPaveText *pt2 = new TPaveText(0.125, 0.75, 0.3, 0.9, "blNDC");
  pt1->SetFillColor(0); pt1->SetBorderSize(1);
  pt2->SetFillColor(0); pt2->SetBorderSize(1);

  Int_t qeInd = 0; Int_t npts1 = 0; Int_t npts2 = 0;
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(snailList[qeInd] > snlNum || qeInd >= qeList.size()) continue;
    if(sign == 0 || snlNum == 150 || snlNum == 151 || snlNum == 159 || snlNum == 160 || snlNum == 220 || snlNum == 221){
      qeInd++;
      continue;
    }
    printf("Snail List Number: %i; Snail Tree Number: %i; QE: %.4f\n", snailList[qeInd], snlNum, qeList[qeInd]);
    if(snlNum <= rangeUpper){
      g1->SetPoint(npts1, qeList[qeInd], sign*1000*pol0.mean);
      g1->SetPointError(npts1++, 0.0, 1000*pol0.meanErr);
    }
    else{
      g2->SetPoint(npts2, qeList[qeInd], sign*1000*pol0.mean);
      g2->SetPointError(npts2++, 0.0, 1000*pol0.meanErr);
    }
    qeInd++;
  }

  // for(Int_t i = 0; i < rangeMsmts1; i++){
  //   g1->SetPoint(i, qe1[i], pol0Mean1[i]);
  //   g1->SetPointError(i, 0.0, pol0Err1[i]);
  // }
  // for(Int_t i = 0; i < rangeMsmts2; i++){
  //   g2->SetPoint(i, qe2[i], pol0Mean2[i]);
  //   g2->SetPointError(i, 0.0, pol0Err2[i]);
  // }

  g1->Fit(fit1, "", "Q");
  g2->Fit(fit2, "", "Q");
  pt1->AddText("--------Fit Parameters--------");
  pt1->AddText(Form("p0: %.4f +/- %.4f", fit1->GetParameter(0), fit1->GetParError(0)));
  pt1->AddText(Form("p1: %.4f +/- %.4f", fit1->GetParameter(1), fit1->GetParError(1)));
  pt1->AddText(Form("Chi2 / NDF: %.3f / %i", fit1->GetChisquare(), fit1->GetNDF()));
  pt2->AddText("--------Fit Parameters--------");
  pt2->AddText(Form("p0: %.4f +/- %.4f", fit2->GetParameter(0), fit2->GetParError(0)));
  pt2->AddText(Form("p1: %.4f +/- %.4f", fit2->GetParameter(1), fit2->GetParError(1)));
  pt2->AddText(Form("Chi2 / NDF: %.3f / %i", fit2->GetChisquare(), fit2->GetNDF()));

  c1->cd()->SetGridx();
  c1->cd()->SetGridy();
  g1->Draw("ap");
  pt1->Draw("same");
  c2->cd()->SetGridx();
  c2->cd()->SetGridy();
  g2->Draw("ap");
  pt2->Draw("same");
}
