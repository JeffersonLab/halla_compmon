#include "../grandOnline/makePlots.h"
#include <ctime>


struct std::tm start = {0, 0, 0, 0, 0, 119};

Int_t daysInMonth(Int_t monthNum, Int_t yearNum){
  if(monthNum == 1 || monthNum == 3 || monthNum == 5 || monthNum == 7 || monthNum == 8 || monthNum == 10 || monthNum == 12){
    return 31;
  }
  else if(monthNum == 4 || monthNum == 6 || monthNum == 9 || monthNum == 11){
    return 30;
  }
  else if(monthNum == 2 && yearNum % 4 == 0 && yearNum % 100 != 0){
    return 29;
  }
  else if(monthNum == 2 && yearNum % 400 == 0){
    return 29;
  }
  else if(monthNum == 2){
    return 28;
  }
  else{
    printf("ERROR: Invalid month entered!");
    return 0;
  }
}


void anPowTime(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *run = (TTree *)f->Get("run");

  TCanvas *c = new TCanvas("c", "Analyzing Power Canvas", 1500, 600);
  c->SetGridx();
  c->SetGridy();
  TGraph *g = new TGraph();
  g->SetTitle("Analyzing Power vs Time");
  g->GetXaxis()->SetTitle("Days After 2019-01-01");
  g->GetYaxis()->SetTitle("Analyzing Power [ppt]");
  g->SetMarkerStyle(3);

  PolVar anPow, collOff;
  FitPolVar pol0;
  Int_t year, month, day, hour, minute, second;
  Int_t runNum, nCycles;

  run->SetBranchAddress("AnalyzingPower", &anPow);
  run->SetBranchAddress("CollimatorOffset", &collOff);
  run->SetBranchAddress("year", &year);
  run->SetBranchAddress("month", &month);
  run->SetBranchAddress("day", &day);
  run->SetBranchAddress("hour", &hour);
  run->SetBranchAddress("minute", &minute);
  run->SetBranchAddress("second", &second);
  run->SetBranchAddress("runNum", &runNum);
  run->SetBranchAddress("numCyclesAcc", &nCycles);
  run->SetBranchAddress("Pol0", &pol0);

  Float_t avgTerm = 0.0;
  Float_t weightSum = 0.0;
  Int_t nPts = 0;
  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    if(pol0.meanErr == 0){continue;}
    struct std::tm runTime = {second, minute, hour, day-1, month-1, year-1900};
    double diff = std::difftime(std::mktime(&runTime), std::mktime(&start))/(60 * 60 * 24);

    g->SetPoint(nPts++, diff, 1000*anPow.mean);

    avgTerm += TMath::Power(1.0/pol0.meanErr, 2)*anPow.mean;
    weightSum += TMath::Power(1.0/pol0.meanErr, 2);
  }

  Float_t errTerm = 0.0;
  Float_t p0 = 16.5612/1000.0;
  Float_t p1 = -0.0102/1000.0;
  Float_t p2 = 0.0159/1000.0;
  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    if(pol0.meanErr == 0){continue;}
    Float_t anPowLo = p0 + p1*(collOff.mean - collOff.meanErr) + p2*TMath::Power(collOff.mean - collOff.meanErr, 2);
    Float_t anPowHi = p0 + p1*(collOff.mean + collOff.meanErr) + p2*TMath::Power(collOff.mean + collOff.meanErr, 2);

    Float_t diff1 = anPow.mean - anPowLo;
    Float_t diff2 = anPowHi - anPow.mean;

    Float_t anPowErr = (diff1 > diff2) ? diff1 : diff2;
    Float_t relWeight = (1.0/TMath::Power(pol0.meanErr, 2))/weightSum;

    errTerm += anPowErr*relWeight;
  }

  printf("Analyzing Power Average: %.4f +/- %.4f\n", 1000*avgTerm/weightSum, 1000.0*errTerm);

  c->cd();
  g->Draw("ap");
}
