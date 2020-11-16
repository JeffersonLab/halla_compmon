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

void cycleCutPct(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  TTree *snl = (TTree *)f->Get("snl");

  Int_t snailNum, runNum, cycNum, numCycles;

  snl->SetBranchAddress("snailNum", &snailNum);
  snl->SetBranchAddress("numCycles", &numCycles);

  Int_t startSnail = (prexOrCrex == 1) ? 1 : 101;

  TH1F *hCut = new TH1F("hCut", Form("Pct Cycles Cut by Snail (%s)", expt.Data()),
                        snl->GetEntries(), startSnail, startSnail + snl->GetEntries());

  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    Int_t nRejCycles = (Int_t)cyc->GetEntries(Form("snailNum==%i && CycleCut==1", snailNum));
    hCut->SetBinContent(i+1, nRejCycles*1.0/numCycles);
  }

  TCanvas *cCut = new TCanvas("cCut", "Cycle Cuts", 1200, 800);
  cCut->cd();
  hCut->SetStats(0);
  hCut->Draw("h");
}
