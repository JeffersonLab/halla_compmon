#include "../grandOnline/makePlots.h"

using namespace std;

/**
vector<vector<Float_t>> readBPMsFile(Int_t runNum){
  vector<vector<Float_t>> cycBPMs;
  //vector<Float_t> bpmMeans, bpmErrs;
  ifstream infile(Form("%s/Run%i_bpms.csv", getenv("COMPMON_RUNPLOTS"), runNum));
  string readStr;
  while(getline(infile, readStr)){
    stringstream ss(readStr);
    string token;
    vector<Float_t> bpmLine;
    while(getline(ss, token, ',')){
      bpmLine.push_back(atof(token.c_str()));
    }
    cycBPMs.push_back(bpmLine);
  }

  return cycBPMs;
}
**/

vector<Float_t> calcCollOffset(Int_t prexOrCrex, vector<vector<Float_t>> bpms){
  Float_t collCentX = 0.0; Float_t collCentXErr = 0.0;
  Float_t collCentY = 0.0; Float_t collCentYErr = 0.0;
  if(prexOrCrex == 1){
    collCentX = 2.160; collCentXErr = 0.080;
    collCentY = 0.709; collCentYErr = 0.098;
  }
  else if(prexOrCrex == 2){
    collCentX = -0.133; collCentXErr = 0.038;
    collCentY = 0.747; collCentYErr = 0.021;
  }

  Float_t collPosX = bpms[0][0] + 6*(bpms[2][0] - bpms[0][0]); 
  Float_t collPosY = bpms[1][0] + 6*(bpms[3][0] - bpms[1][0]);

  Float_t xDiffErr = 6*TMath::Sqrt(TMath::Power(bpms[0][1], 2) + TMath::Power(bpms[2][1], 2));
  Float_t yDiffErr = 6*TMath::Sqrt(TMath::Power(bpms[1][1], 2) + TMath::Power(bpms[3][1], 2));
  Float_t xProjErr = TMath::Sqrt(TMath::Power(bpms[0][1], 2) + TMath::Power(xDiffErr, 2));
  Float_t yProjErr = TMath::Sqrt(TMath::Power(bpms[1][1], 2) + TMath::Power(yDiffErr, 2));

  Float_t collOffX = collPosX - collCentX;
  Float_t collOffY = collPosY - collCentY;
  Float_t collOffXErr = TMath::Sqrt(TMath::Power(xProjErr, 2) + TMath::Power(collCentXErr, 2));
  Float_t collOffYErr = TMath::Sqrt(TMath::Power(yProjErr, 2) + TMath::Power(collCentYErr, 2));

  Float_t collOffset = TMath::Sqrt(TMath::Power(collOffX, 2) + TMath::Power(collOffY, 2));
  Float_t offXsqErr = TMath::Abs(TMath::Power(collOffX, 2))*2*collOffXErr/collOffX;
  Float_t offYsqErr = TMath::Abs(TMath::Power(collOffY, 2))*2*collOffYErr/collOffY;
  Float_t offSqSumErr = TMath::Sqrt(TMath::Power(offXsqErr, 2) + TMath::Power(offYsqErr, 2));
  Float_t collOffsetErr = TMath::Abs(collOffset)*0.5*offSqSumErr/(TMath::Power(collOffX, 2) + TMath::Power(collOffY, 2));

  vector<Float_t> offset;
  offset.push_back(collOffset);
  offset.push_back(collOffsetErr);
  offset.push_back(collOffX);
  offset.push_back(collOffXErr);
  offset.push_back(collOffY);
  offset.push_back(collOffYErr);
  return offset;
}

void bpmCalibrationCheck(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *run = (TTree *)f->Get("run");

  Int_t runNum;
  const Int_t nBPMs = 4;
  DataVar bpms[nBPMs];
  PolVar epicsBPMs[nBPMs];
  TString bpmNames[nBPMs] = {"bpmAx", "bpmAy", "bpmBx", "bpmBy"};
  TString epicsNames[nBPMs] = {"epics_bpmAx", "epics_bpmAy", "epics_bpmBx", "epics_bpmBy"};

  run->SetBranchAddress("runNum", &runNum);
  for(Int_t i = 0; i < nBPMs; i++){
    run->SetBranchAddress(bpmNames[i].Data(), &bpms[i]);
    run->SetBranchAddress(epicsNames[i].Data(), &epicsBPMs[i]);
  }
  TH1F *daqOffset = new TH1F("daqOffset", "DAQ vs EPICS: Collimator Offset (runwise)", 100, 0.0, 8.0);
  TH1F *epicsOffset = new TH1F("epicsOffset", "DAQ vs EPICS: Collimator Offset (runwise)", 100, 0.0, 8.0);
  TH1F *daqOffsetX = new TH1F("daqOffsetX", "DAQ vs EPICS: Collimator X Position (runwise)", 100, -8.0, 8.0);
  TH1F *epicsOffsetX = new TH1F("epicsOffsetX", "DAQ vs EPICS: Collimator X Position (runwise)", 100, -8.0, 8.0);
  TH1F *daqOffsetY = new TH1F("daqOffsetY", "DAQ vs EPICS: Collimator Y Position (runwise)", 100, -8.0, 8.0);
  TH1F *epicsOffsetY = new TH1F("epicsOffsetY", "DAQ vs EPICS: Collimator Y Position (runwise)", 100, -8.0, 8.0);

  vector<vector<Float_t>> daqBPMs;
  vector<vector<Float_t>> epicsBPMInfo;
  for(Int_t i = 0; i < nBPMs; i++){
    vector<Float_t> bpmData;
    vector<Float_t> epicsData;
    for(Int_t j = 0; j < 2; j++){
      bpmData.push_back(0.0);
      epicsData.push_back(0.0);
    }
    daqBPMs.push_back(bpmData);
    epicsBPMInfo.push_back(epicsData);
  }

  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    for(Int_t j = 0; j < nBPMs; j++){
      daqBPMs[j][0] = bpms[j].mean;
      daqBPMs[j][1] = bpms[j].meanErr;
      epicsBPMInfo[j][0] = epicsBPMs[j].mean;
      epicsBPMInfo[j][1] = epicsBPMs[j].meanErr;
    }

    vector<Float_t> daqOffData = calcCollOffset(prexOrCrex, daqBPMs);
    vector<Float_t> epicsOffData = calcCollOffset(prexOrCrex, epicsBPMInfo);
    daqOffset->Fill(daqOffData[0]);
    epicsOffset->Fill(epicsOffData[0]);
    daqOffsetX->Fill(daqOffData[2]);
    epicsOffsetX->Fill(epicsOffData[2]);
    daqOffsetY->Fill(daqOffData[4]);
    epicsOffsetY->Fill(epicsOffData[4]);
  }

  TH1F *daqHists[3] = {daqOffset, daqOffsetX, daqOffsetY};
  TH1F *epicsHists[3] = {epicsOffset, epicsOffsetX, epicsOffsetY};

  for(Int_t i = 0; i < 3; i++){
  TCanvas *cOff = new TCanvas(Form("cOff_%i", i + 1), "Offset Canvas", 1200, 800);
    cOff->cd()->SetGridx();
    cOff->cd()->SetGridy();
    daqHists[i]->SetStats(0); daqHists[i]->SetLineColor(kBlue);
    daqHists[i]->GetXaxis()->SetTitle("Coll. Offset [mm]");
    daqHists[i]->Draw();
    daqHists[i]->SetStats(0); epicsHists[i]->SetLineColor(kRed);
    epicsHists[i]->Draw("same");

    TLegend *leg = new TLegend(0.75, 0.75, 0.9, 0.9, "", "NDC");
    leg->AddEntry(daqHists[i], "DAQ BPM Data");
    leg->AddEntry(epicsHists[i], "EPICS BPM Data");
    leg->Draw("same");
  }

}
