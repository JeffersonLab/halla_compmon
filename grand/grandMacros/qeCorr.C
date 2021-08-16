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


void qeCorrGraphs(TF1 *fit1, TF1 *fit2, Bool_t residuals=kFALSE, Float_t yInt1=0.0, Float_t slope1=0.0, Float_t yInt2=0.0, Float_t slope2=0.0){
  TFile *f = TFile::Open(Form("%s/aggregates/crexGrandCompton.root", getenv("COMPMON_WEB")));
  TTree *snl = (TTree *)f->Get("snl");
  TTree *cyc = (TTree *)f->Get("cyc");

  FitPolVar pol0;
  Int_t snlNum, sign;
  Float_t ihwp;
  const Int_t rangeUpper = 137;
  vector<vector<Float_t>> fullList = readQEFile();
  vector<Int_t> snailList = formSnailList(fullList);
  vector<Float_t> qeList = formQEList(fullList);
  TString name   = residuals ? "Residuals" : "Correlations";
  TString title  = residuals ? "Pol0-QE Residuals" : "QE vs Pol0 Correlation";
  TString xTitle = residuals ? "Acc0LasOff/(Acc0LasOn - Acc0LasOff)" : "qe";
  TString yTitle = residuals ? "Pol0 Residuals" : "Pol0 [pct]";

  snl->SetBranchAddress("snailNum", &snlNum);
  snl->SetBranchAddress("sign", &sign);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("Pol0", &pol0);

  TCanvas *c1 = new TCanvas(Form("c1_%s", name.Data()), Form("%s Before Spot Move", name.Data()), 1200, 800);
  TCanvas *c2 = new TCanvas(Form("c2_%s", name.Data()), Form("%s After Spot Move", name.Data()), 1200, 800);

  TGraphErrors *g1 = new TGraphErrors();
  TGraphErrors *g2 = new TGraphErrors();
  g1->SetTitle(Form("%s Before Spot Move", title.Data()));
  g1->GetXaxis()->SetTitle(xTitle.Data());
  g1->GetYaxis()->SetTitle(yTitle.Data());
  g1->SetMarkerStyle(20);
  g2->SetTitle(Form("%s After Spot Move", title.Data()));
  g2->GetXaxis()->SetTitle(xTitle.Data());
  g2->GetYaxis()->SetTitle(yTitle.Data());
  g2->SetMarkerStyle(20);
  
  TPaveText *pt1 = new TPaveText(0.125, 0.75, 0.3, 0.9, "blNDC");
  TPaveText *pt2 = new TPaveText(0.125, 0.75, 0.3, 0.9, "blNDC");
  pt1->SetFillColor(0); pt1->SetBorderSize(1);
  pt2->SetFillColor(0); pt2->SetBorderSize(1);

  Int_t qeInd = 0; Int_t npts1 = 0; Int_t npts2 = 0;
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(snailList[qeInd] > snlNum || qeInd >= qeList.size()) continue;
    if(sign == 0){
      qeInd++;
      continue;
    }
    TString hNameOn = Form("h%sAcc0LasOn_snail%i", name.Data(), snlNum);
    TString hNameOff1 = Form("h%sAcc0LasOff1_snail%i", name.Data(), snlNum);
    TString hNameOff2 = Form("h%sAccoLasOff2_snail%i", name.Data(), snlNum);
    cyc->Project(hNameOn.Data(), "Acc0LasOn", Form("CycleCut==0 && snailNum==%i", snlNum), "goff");
    cyc->Project(hNameOff1.Data(), "Acc0LasOff1", Form("CycleCut==0 && snailNum==%i", snlNum), "goff");
    cyc->Project(hNameOff2.Data(), "Acc0LasOff2", Form("CycleCut==0 && snailNum==%i", snlNum), "goff");
    TH1F *hOn = (TH1F *)gDirectory->Get(hNameOn.Data());
    TH1F *hOff1 = (TH1F *)gDirectory->Get(hNameOff1.Data());
    TH1F *hOff2 = (TH1F *)gDirectory->Get(hNameOff2.Data());
    Float_t lasOn = hOn->GetMean();
    Float_t lasOff = (hOff1->GetMean() + hOff2->GetMean())/2.0;
    
    Float_t fitVal1 = slope1*qeList[qeInd] + yInt1;
    Float_t fitVal2 = slope2*qeList[qeInd] + yInt2;
    Float_t xVal  = residuals ? lasOff*1.0/(lasOn - lasOff) : qeList[qeInd];
    Float_t yVal1 = residuals ? (100.0*pol0.mean*sign - fitVal1) : sign*100.0*pol0.mean;
    Float_t yVal2 = residuals ? (100.0*pol0.mean*sign - fitVal2) : sign*100.0*pol0.mean;
    if(residuals){
      printf("Snail List Number: %i; Snail Tree Number: %i; Acc0LasOn: %.4f; Acc0LasOff: %.4f; X-Val: %.4f\n", snailList[qeInd], snlNum, lasOn, lasOff, xVal);
    }
    else{
      printf("Snail List Number: %i; Snail Tree Number: %i; QE: %.4f\n", snailList[qeInd], snlNum, qeList[qeInd]);
    }
    if(snlNum <= rangeUpper){
      g1->SetPoint(npts1, xVal, yVal1);
      g1->SetPointError(npts1, 0.0, 100.0*pol0.meanErr);
      npts1++;
    }
    else{
      g2->SetPoint(npts2, xVal, yVal2);
      g2->SetPointError(npts2, 0.0, 100.0*pol0.meanErr);
      npts2++;
    }
    qeInd++;
  }

  g1->Fit(fit1, "Q", "");
  g2->Fit(fit2, "Q", "");
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

void qeCorr(){
  TF1 *fit1 = new TF1("f1_residuals", "pol1");
  TF1 *fit2 = new TF1("f2_residuals", "pol1");
  TF1 *fit3 = new TF1("f3_correlations", "pol1");
  TF1 *fit4 = new TF1("f4_correlations", "pol1");

  qeCorrGraphs(fit1, fit2);
  printf("Fit Parameters: %.4f, %.4f, %.4f, %.4f\n", fit1->GetParameter(0), fit1->GetParameter(1), fit2->GetParameter(0), fit2->GetParameter(1));
  qeCorrGraphs(fit3, fit4, kTRUE, fit1->GetParameter(0), fit1->GetParameter(1), fit2->GetParameter(0), fit2->GetParameter(1));
}
