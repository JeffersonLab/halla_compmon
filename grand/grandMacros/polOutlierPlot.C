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
#include "makePlots.h"

void polOutlierPlot(Int_t prexOrCrex){
  //Open files and get trees
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  TTree *run = (TTree *)f->Get("run");
  TTree *snl = (TTree *)f->Get("snl");

  PolVar pol0;
  FitPolVar snlPol, runPol;
  Int_t snail, snailNum, runNum, cycNum, cycCut;
  Float_t hWien, vWien, solWien, ihwp;
  Int_t runSnail, runTreeNum, nCycles, runSign;
  Float_t hWienRun, vWienRun, solWienRun, ihwpRun;

  //Note which trees are being set to which variable
  snl->SetBranchAddress("snailNum", &snail);
  snl->SetBranchAddress("Pol0", &snlPol);
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("Pol0", &pol0);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  run->SetBranchAddress("snailNum", &runSnail);
  run->SetBranchAddress("runNum", &runTreeNum);
  run->SetBranchAddress("HWienAngle", &hWienRun);
  run->SetBranchAddress("VWienAngle", &vWienRun);
  run->SetBranchAddress("PhiFG", &solWienRun);
  run->SetBranchAddress("ihwp", &ihwpRun);
  run->SetBranchAddress("Pol0", &runPol);
  run->SetBranchAddress("sign", &runSign);
  run->SetBranchAddress("numCyclesAcc", &nCycles);
  snl->SetBranchAddress("HWienAngle", &hWien);
  snl->SetBranchAddress("VWienAngle", &vWien);
  snl->SetBranchAddress("PhiFG", &solWien);
  snl->SetBranchAddress("ihwp", &ihwp);

  Int_t startSnail = prexOrCrex == 1 ? 0 : 100;
  //Float_t sectionPols[4] = {0.89461, -0.88599, -0.88837, 0.87826};
  //Float_t sectionErrs[4] = {0.00182,  0.00178, 0.00426,  0.00397};
  Float_t sectionPols[4] = {0.89149, -0.88261, -0.88260, 0.87420};
  Float_t sectionErrs[4] = {0.00182,  0.00178,  0.00426, 0.00397};

  //hPull == Plot for cycles pulled against the polarization of their snail
  //hSecPull == Plot for cycles pulled against the average polarization of the wien/ihwp region they're in
  //hSnlPull == Plot for snails pulled against the average polarization of the wien/ihwp region they're in
  TH1F *hPull = new TH1F("hPull", Form("Polarization Global Pull Plot (%s)", expt.Data()), 100, -10, 10);
  TH1F *hSecPull = new TH1F("hSecPull", Form("Polarization Pull Plot (%s, by wien avg)", expt.Data()), 100, -10, 10);
  TH1F *hRunPull = new TH1F("hRunPull", Form("Polarization Run Pull Plot (%s)", expt.Data()), 320, 4300, 4620);
  TH1F *hSnlPull = new TH1F("hSnlPull", Form("Polarization Snail Pull Plot (%s)", expt.Data()), 100, -6, 6);

  //Wien and IHWP Indices tell us where to look for average polarization
  Int_t curSnail = 0; Float_t curPol = 0.0; Float_t curPolErr = 0.0;
  Int_t wienInd = 0; Int_t ihwpInd = 0; Int_t maxFillSnl = 0;
  //Loop thru cycles
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(snailNum != curSnail){
      //If you're starting a new snail, loop thru snail tree to find snail's polarization
      for(Int_t j = 0; j < snl->GetEntries(); j++){
        snl->GetEntry(j);
        if(snail == snailNum){
          curPol = snlPol.mean; curPolErr = snlPol.meanErr;
          curSnail = snail;
          wienInd = isFlipLeft(snailNum, hWien, vWien, solWien, true) == true ? 0 : 2;
          ihwpInd = (Int_t)(ihwp > 0.98);
          if(snailNum > maxFillSnl){
            Float_t snlErr = TMath::Sqrt(TMath::Power(snlPol.meanErr, 2) + TMath::Power(sectionErrs[wienInd + ihwpInd], 2));
            hSnlPull->Fill((snlPol.mean - sectionPols[wienInd + ihwpInd])/(snlErr));
            maxFillSnl = snailNum;
          }
          break;
        }
      }
    }
    //Actual calculation of pulls
    Float_t polDiffErr = TMath::Sqrt(TMath::Power(pol0.meanErr, 2) + TMath::Power(curPolErr, 2));
    Float_t polDiff = pol0.mean - curPol;
    Float_t sectionDiffErr = TMath::Sqrt(TMath::Power(pol0.meanErr, 2) + TMath::Power(sectionErrs[wienInd + ihwpInd], 2));
    Float_t sectionDiff = pol0.mean - sectionPols[wienInd + ihwpInd];
    if(cycCut==0){
      hPull->Fill(polDiff/polDiffErr);
      hSecPull->Fill(sectionDiff/sectionDiffErr);
    }
  }
  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    wienInd = isFlipLeft(runTreeNum, hWienRun, vWienRun, solWienRun, false) == true ? 0 : 2;
    ihwpInd = (Int_t)(ihwpRun > 0.98);
    //Float_t runPull = runPol.mean - sectionPols[wienInd + ihwpInd];
    //Float_t runPullErr = TMath::Sqrt(TMath::Power(runPol.meanErr, 2) + TMath::Power(sectionErrs[wienInd + ihwpInd], 2));
    Float_t runPull = runPol.mean*runSign - 0.889;
    Float_t runPullErr = runPol.meanErr;
    if(nCycles > 0 && TMath::Power(runPull/runPullErr, 2) > 8){
      //hRunPull->Fill(runPull/runPullErr);
      //hRunPull->Fill(TMath::Power(runPull/runPullErr, 2));
      hRunPull->Fill(runTreeNum);
    }
  }

  //Build canvas
  TCanvas *cPull = new TCanvas("cPull", "Pull Plot", 1200, 1200);
  cPull->Divide(2, 2);
  cPull->cd(1); hPull->SetStats(220111); hPull->Draw();
  cPull->cd(2); hSecPull->SetStats(220111); hSecPull->Draw();
  cPull->cd(3); hRunPull->SetStats(220111); hRunPull->Draw();
  cPull->cd(4); hSnlPull->SetStats(220111); hSnlPull->Draw();
}

