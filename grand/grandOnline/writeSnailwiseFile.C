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

#include "../grandConstruction/vars.h"

using namespace std;

TString ihwpStateName(Float_t ihwp){
  if(ihwp == 0){return "OUT";}
  else if(ihwp == 1){return "IN";}
  else{return "UNK";}
}

TString signString(Int_t sign){
  if(sign == 1){return " 1";}
  else if(sign == -1){return "-1";}
  else{return " 0";}
}

void writeSnailwiseFile(Int_t prexOrCrex){
  printf("Writing summary textfile...\n");
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
  
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), exptStr.Data()));
  TTree *snl = (TTree *)f->Get("snl");

  FitPolVar pol0, asymOff;
  Int_t sign, snlNum, year, month, day, hour, minute, second;
  Float_t qw1, hw1, ihwp, vWien, hWien, solWien;

  snl->SetBranchAddress("Pol0", &pol0);
  snl->SetBranchAddress("Asym0LasOff", &asymOff);
  snl->SetBranchAddress("qw1", &qw1);
  snl->SetBranchAddress("hw1", &hw1);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("sign", &sign);
  snl->SetBranchAddress("VWienAngle", &vWien);
  snl->SetBranchAddress("HWienAngle", &hWien);
  snl->SetBranchAddress("PhiFG", &solWien);
  snl->SetBranchAddress("snailNum", &snlNum);
  snl->SetBranchAddress("year", &year);
  snl->SetBranchAddress("month", &month);
  snl->SetBranchAddress("day", &day);
  snl->SetBranchAddress("hour", &hour);
  snl->SetBranchAddress("minute", &minute);
  snl->SetBranchAddress("second", &second);

  ofstream outfile(Form("%s/aggregates/%sCompton.csv", getenv("COMPMON_WEB"), exptStr.Data()));
  outfile<<Form("//Pol0 (pct), Pol0 Err (pct), Chi2, NDF, ihwp state, vert Wien Angle, horiz Wien Angle, solenoid Wien Angle, calculated sign, qw1 setting, hw1 setting\n");
  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    TString polMsmt = Form("%.2f,%.2f,%3.2f,%03i", 100*pol0.mean, 100*pol0.meanErr, pol0.Chi2, pol0.NDF);
    TString signStr = Form("%s,%3.4f,%3.4f,%3.4f,%i", ihwpStateName(ihwp).Data(), vWien, hWien, solWien, sign);
    TString laser = Form("%05i,%05i", (Int_t)qw1, (Int_t)hw1);
    TString date = Form("%04i-%02i-%02i", year, month, day);
    TString time = Form("%02i:%02i:%02i", hour, minute, second);
    outfile<<Form("%03i,%s,%s,%s,%s,%s\n", snlNum, polMsmt.Data(), signStr.Data(), laser.Data(), date.Data(), time.Data());
  }
  outfile.close();
  printf("...Done!\n");
}
