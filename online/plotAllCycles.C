#include <TMath.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <THStack.h>
#include <TPad.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "utils.h"

using namespace std;

vector<vector<int>> completeCycleList(int runNum){
  ifstream cycle_infile(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), runNum));
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

void plotAllCycles(int runNum){
  vector<vector<int>> allCycles = completeCycleList(runNum);
  //TFile *f = new TFile(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum), "READ");
  //TTree *mpswise = (TTree *)f->Get("mpswise");
  //TTree *quartetwise = (TTree *)f->Get("quartetwise");
  vector<TChain *> runChains = loadChain(runNum);
  Int_t mpsInd = 0 ; Int_t quartetInd = 1;
  TChain *mpswise = runChains[mpsInd]; TChain *quartetwise = runChains[quartetInd];

  gStyle->SetOptStat(0);
  Float_t sumOnMean = 0; Float_t sumOffMean = 0;
  for(int i = 0; i < allCycles.size(); i++){
    vector<int> cycle = allCycles[i];
    printf("Cycle %i Limits: %i, %i, %i, %i, %i, %i\n", i + 1, cycle[0], cycle[1], cycle[2], cycle[3], cycle[4], cycle[5]);
    TCanvas *currentCan = new TCanvas(Form("cycle_%i_pad", i), Form("Run %i, Laser Cycle %i", runNum, i + 1), 1200, 800);
    currentCan->Divide(2, 2);
    
    TString pos("PosHelAcc0/PosHelNSamples0");
    TString neg("NegHelAcc0/NegHelNSamples0");
    TString msmts[3] = {Form("%s - %s", pos.Data(), neg.Data()), Form("%s + %s", pos.Data(), neg.Data()),
                      Form("(%s - %s)/(%s + %s)", pos.Data(), neg.Data(), pos.Data(), neg.Data())};
    TString ids[3] = {"Diffs", "Sums", "Asyms"};
    TString cycleParts[3] = {"preOff", "on", "postOff"};
    for(int j = 0; j < 3; j++){
      TString hPreOff_name = Form("h_cycle_%i_%s_%s", i, ids[j].Data(), cycleParts[0].Data());
      TString hOn_name = Form("h_cycle_%i_%s_%s", i, ids[j].Data(), cycleParts[1].Data());
      TString hPostOff_name = Form("h_cycle_%i_%s_%s", i, ids[j].Data(), cycleParts[2].Data());
      TString onCut("(laserState==0 || laserState==1) && beamState==1");
      TString offCut("(laserState==2 || laserState==3) && beamState==1");
      //printf("Cut1: %s; Cut 2: %s\n", onCut.Data(), offCut.Data());
      if(j != 2){
        quartetwise->Project(hPreOff_name.Data(), msmts[j].Data(), Form("%s && firstMPSnumber >= %i && firstMPSnumber <= %i", offCut.Data(), cycle[0], cycle[1]));
        quartetwise->Project(hOn_name.Data(), msmts[j].Data(), Form("%s && firstMPSnumber >= %i && firstMPSnumber <= %i", onCut.Data(), cycle[2], cycle[3]));
        quartetwise->Project(hPostOff_name.Data(), msmts[j].Data(), Form("%s && firstMPSnumber >= %i && firstMPSnumber <= %i", offCut.Data(), cycle[4], cycle[5]));
      }
      else{
        TString newMsmtOn  = Form("(%s - %s)/(%s + %s - %f)", pos.Data(), neg.Data(), pos.Data(), neg.Data(), sumOffMean);
        TString newMsmtOff = Form("(%s - %s)/(%f - %f)", pos.Data(), neg.Data(), sumOnMean, sumOffMean);
        quartetwise->Project(hPreOff_name.Data(), newMsmtOff.Data(), Form("%s && firstMPSnumber >= %i && firstMPSnumber <= %i", offCut.Data(), cycle[0], cycle[1]));
        quartetwise->Project(hOn_name.Data(), newMsmtOn.Data(), Form("%s && firstMPSnumber >= %i && firstMPSnumber <= %i", onCut.Data(), cycle[2], cycle[3]));
        quartetwise->Project(hPostOff_name.Data(), newMsmtOff.Data(), Form("%s && firstMPSnumber >= %i && firstMPSnumber <= %i", offCut.Data(), cycle[4], cycle[5]));
      }
      //quartetwise->Draw(Form("%s>>%s", msmts[j].Data(), hPreOff_name.Data()),
      //                  Form("(laserState==2 || laserState==3) && beamState==1 && firstMPSnumber>=%i && firstMPSnumber<%i", cycle[0], cycle[1]), "goff");
      //quartetwise->Draw(Form("%s>>%s", msmts[i].Data(), hOn_name.Data()),
      //                  Form("(laserState==0 || laserState==1) && beamState==1 && firstMPSnumber>=%i && firstMPSnumber<%i",  cycle[2], cycle[3]), "goff");
      //quartetwise->Draw(Form("%s>>%s", msmts[i].Data(), hPostOff_name.Data()),
      //                  Form("(laserState==2 || laserState==3) && beamState==1 && firstMPSnumber>=%i && firstMPSnumber<%i", cycle[4], cycle[5]), "goff");

      TH1F *hPreOff = (TH1F *)gDirectory->Get(hPreOff_name.Data());
      TH1F *hOn = (TH1F *)gDirectory->Get(hOn_name.Data());
      TH1F *hPostOff = (TH1F *)gDirectory->Get(hPostOff_name.Data());
      hPreOff->SetLineColor(kRed); hOn->SetLineColor(kGreen + 2); hPostOff->SetLineColor(kRed + 2);

      if(j == 1){
        sumOnMean = hOn->GetMean();
        sumOffMean = (hPreOff->GetMean() +hPostOff->GetMean())/2.0;
      }

      currentCan->cd(j + 1);
      TPad *pName  = new TPad(Form("p_%s_div%i_name", msmts[j].Data(), i),  "Hist Title", 0.0, 0.9, 1.0, 1.0);
      TPad *pStat1 = new TPad(Form("p_%s_div%i_stat1", msmts[j].Data(), i), "Stats 1", 0.0, 0.7, 0.33, 0.9);
      TPad *pStat2 = new TPad(Form("p_%s_div%i_stat2", msmts[j].Data(), i), "Stats 2", 0.33, 0.7, 0.67, 0.9);
      TPad *pStat3 = new TPad(Form("p_%s_div%i_stat3", msmts[j].Data(), i), "Stats 3", 0.67, 0.7, 1.0, 0.9);
      TPad *pHist  = new TPad(Form("p_%s_div%i_hist", msmts[j].Data(), i),  "Histogram", 0.0, 0.0, 1.0, 0.7);
      pName->Draw(); pStat1->Draw(); pStat2->Draw(); pStat3->Draw(); pHist->Draw();
      
      TPaveText *ptPreOff = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptPreOff->SetBorderSize(1); ptPreOff->SetFillColor(0);
      TPaveText *ptOn   = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptOn->SetBorderSize(1);   ptOn->SetFillColor(0);
      TPaveText *ptPostOff = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptPostOff->SetBorderSize(1); ptPostOff->SetFillColor(0);
      vector<TString> sPreOff = hist_stats(hPreOff); vector<TString> sOn = hist_stats(hOn); vector<TString> sPostOff = hist_stats(hPostOff);
      for(int k = 0; k < sPreOff.size(); k++){ptPreOff->AddText(sPreOff[k].Data())->SetTextColor(kRed);}
      for(int k = 0; k < sOn.size(); k++){ptOn->AddText(sOn[k].Data())->SetTextColor(kGreen + 2);}
      for(int k = 0; k < sPostOff.size(); k++){ptPostOff->AddText(sPostOff[k].Data())->SetTextColor(kRed + 2);}
      pStat1->cd(); ptPreOff->Draw(); 
      pStat2->cd(); ptOn->Draw(); 
      pStat3->cd(); ptPostOff->Draw();

      TPaveText *ptName = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptName->SetBorderSize(0); ptName->SetFillColor(16);
      ptName->AddText(Form("Cycle %i, %s", i + 1, ids[j].Data()));
      pName->cd(); ptName->Draw();

      THStack *hs = new THStack(Form("%s_cycle%i_stack", ids[j].Data(), i), Form("Run %i quartetwise, %s", runNum, ids[j].Data()));
      hs->Add(hOn); hs->Add(hPreOff); hs->Add(hPostOff);
      pHist->cd(); hs->Draw("nostack");
    }
    currentCan->cd(4);
    mpswise->Draw("Acc0/NAcc0:mpsCoda", Form("mpsCoda>=%i && mpsCoda<=%i && beamState==1", cycle[0], cycle[5]));
    currentCan->Print(Form("%s/runs/Run%i/cycle_%03i_plots.png", getenv("COMPMON_WEB"), runNum, i + 1), "png");
    
    TCanvas *graphCan = new TCanvas(Form("graph_cycle%i_pad", i), "Graphs Canvas", 1200, 800);
    graphCan->Divide(2, 2);
    TString cycCut = Form("mpsCoda>=%i && mpsCoda<=%i", cycle[0], cycle[5]);
    graphCan->cd(1); mpswise->Draw("laserState:mpsCoda", cycCut.Data());
    graphCan->cd(2); mpswise->Draw("beamState + 1:mpsCoda", cycCut.Data());
    graphCan->cd(3); mpswise->Draw("cavPowerCalibrated:mpsCoda", cycCut.Data());
    graphCan->cd(4); mpswise->Draw("bcm:mpsCoda", cycCut.Data());
    graphCan->Print(Form("%s/runs/Run%i/cycle_%03i.png", getenv("COMPMON_WEB"), runNum, i + 1), "png");
  }
  
  gSystem->Exec(Form("convert %s/runs/Run%i/cycle_*.png %s/runs/Run%i/cycle_qVars.pdf", getenv("COMPMON_WEB"), runNum, getenv("COMPMON_WEB"), runNum));
  gSystem->Exec(Form("rm -rf %s/runs/Run%i/*.png", getenv("COMPMON_WEB"), runNum));
}
