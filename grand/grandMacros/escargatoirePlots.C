#include "makePlots.h"
#include "../vars.h"

void escargatoirePlots(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  if(prexOrCrex != 1){
    printf("Configuration not implemented for CREX yet, please run PREX-only.\n");
    exit(1);
  }

  const Int_t nPairs = 19;
  const Int_t nTriples = 1;
  Int_t pairs[nPairs][3] = {{ 0,  4,  6}, { 0,  5,  7}, { 0,  8, 10}, { 0,  9, 11}, 
                            { 0, 12, 14}, { 0, 13, 15}, { 0, 16, 18}, { 0, 17, 19},
                            { 0, 20, 22}, { 0, 21, 23}, {24, 26, 31}, { 0, 25, 28},
                            {27, 29, 43}, {30, 32, 42}, { 0, 33, 38}, { 0, 34, 36}, 
                            { 0, 35, 37}, { 0, 39, 41}, { 0, 40, 44}};
  //Int_t triples[nTriples][3] = {{27, 29, 31}};

  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  TTree *snl = (TTree *)f->Get("snl");

  const Int_t nGraphs = 4;
  const Int_t nData = 1;
  const Int_t nCans = 2*nData;
  TCanvas *cans[nCans];
  TGraphErrors *graphs[nCans][nGraphs];
  Int_t counts[nCans][nGraphs];
  TPaveText *texts[nData][nGraphs];
  TF1 *fits[nData];
  TF1 *escFits[nData][nGraphs];
  FitPolVar vars[nData];
  TString names[nData] = {"Pol0"};
  TString titles[nData] = {"Polarization (Sign Corrected)"};
  Float_t ymins[nCans]; Float_t ymaxs[nCans];
  Float_t xmins[nGraphs]; Float_t xmaxs[nGraphs];
  Float_t xLo[nGraphs] = {0.10, 0.25, 0.60, 0.75};
  Float_t xHi[nGraphs] = {0.25, 0.40, 0.75, 0.90};
  Float_t yLo[nGraphs] = {0.80, 0.80, 0.80, 0.80};
  Float_t yHi[nGraphs] = {0.90, 0.90, 0.90, 0.90};

  PolVar polVars[nData];
  Int_t cycNum, runNum, snailNum, cycCut;
  Int_t snailTreeNum, nCycles, sign;
  Float_t hWien, vWien, solWien, ihwp, laserPol;

  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  for(Int_t i = 0; i < nData; i++){
    cyc->SetBranchAddress(names[i].Data(), &polVars[i]);
  }
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  snl->SetBranchAddress("snailNum", &snailTreeNum);
  snl->SetBranchAddress("numCyclesAcc", &nCycles);
  snl->SetBranchAddress("HWienAngle", &hWien);
  snl->SetBranchAddress("VWienAngle", &vWien);
  snl->SetBranchAddress("PhiFG", &solWien);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("LaserPolarization", &laserPol);
  snl->SetBranchAddress("sign", &sign);

  
  for(Int_t i = 0; i < nData; i++){
    cans[i] = new TCanvas(Form("c%s", names[i].Data()), Form("%s Canvas", titles[i].Data()), 1200, 800);
    cans[i]->SetGridx(); cans[i]->SetGridy();
    cans[i+1] = new TCanvas(Form("c%s_Chi2", names[i].Data()), Form("%s Chi2/NDF Canvas", titles[i].Data()), 1200, 800);
    cans[i+1]->SetGridx(); cans[i+1]->SetGridy();
    ymins[i] = 1e16; ymaxs[i] = -1e16;
    ymins[i+1] = 1e16; ymaxs[i] = -1e16;
    fits[i] = new TF1(Form("fit_%s", names[i].Data()), "pol0");
    for(Int_t j = 0; j < nGraphs; j++){
      graphs[i][j] = new TGraphErrors();
      graphs[i][j]->SetTitle(Form("%s vs Escargatoire", titles[i].Data()));
      graphs[i][j]->GetXaxis()->SetTitle("Escargatoire Num");
      graphs[i][j]->GetYaxis()->SetTitle(titles[i].Data());      
      graphs[i][j]->SetMarkerColor(getColor(j));
      graphs[i][j]->SetMarkerStyle(getMarkerStyle(j));
      counts[i][j] = 0;
      texts[i][j] = new TPaveText(xLo[j], yLo[j], xHi[j], yHi[j], "blNDC");
      texts[i][j]->SetFillColor(0); texts[i][j]->SetBorderSize(1);
      escFits[i][j] = new TF1(Form("escFit_%s_graph%i", names[i].Data(), j), "pol0");
      escFits[i][j]->SetLineColor(getColor(j));
      xmins[j] = 1e10; xmaxs[j] = -1e10;

      graphs[i+1][j] = new TGraphErrors();
      graphs[i+1][j]->SetTitle(Form("%s Chi2/NDF vs Escargatoire", titles[i].Data()));
      graphs[i+1][j]->GetXaxis()->SetTitle("Escargatoire Num");
      graphs[i+1][j]->GetYaxis()->SetTitle("Chi2/NDF");
      graphs[i+1][j]->GetXaxis()->SetLimits(0, nPairs+1);
      graphs[i+1][j]->SetMarkerColor(getColor(j));
      graphs[i+1][j]->SetMarkerStyle(getMarkerStyle(j));
      counts[i][j] = 0;
    }
  }

  for(Int_t i = 0; i < nPairs; i++){
    Float_t escHWien = 0.0; Float_t escVWien = 0.0; Float_t escSolWien = 0.0; 
    Float_t escIHWP = -1.0; Float_t escLaserPol = 0.0; 
    Int_t escCycles = 0; Int_t escSign = 0;
    for(Int_t j = 0; j < snl->GetEntries(); j++){
      snl->GetEntry(j);
      if(snailTreeNum == pairs[i][0] || snailTreeNum == pairs[i][1] || snailTreeNum == pairs[i][2]){
        escHWien = hWien; escVWien = vWien; escSolWien = solWien;
        escIHWP = ihwp; escLaserPol = laserPol; escCycles += nCycles;
        escSign = sign;
      }
    }

    TH1F *hists[nData];
    Int_t cycCounts[nData];
    for(Int_t j = 0; j < nData; j++){
      hists[j] = new TH1F(Form("h%s_esc%i", names[j].Data(), i+1), "", escCycles, 0, escCycles);
      cycCounts[j] = 1;
    }

    for(Int_t j = 0; j < cyc->GetEntries(); j++){
      cyc->GetEntry(j);
      if((snailNum == pairs[i][0] || snailNum == pairs[i][1] || snailNum == pairs[i][2]) && cycCut == 0){
        for(Int_t k = 0; k < nData; k++){
          hists[k]->SetBinContent(cycCounts[k], polVars[k].mean);
          hists[k]->SetBinError(cycCounts[k]++, polVars[k].meanErr);
        }
      }
    }
  
    for(Int_t j = 0; j < nData; j++){
      hists[j]->Fit(Form("fit_%s", names[j].Data()), "Q", "", 0, escCycles);
      //Int_t ind = getGraphInd(snailTreeNum, escHWien, escVWien, escSolWien, escIHWP);
      Int_t ind = getGraphInd(snailTreeNum, -13.0, -89.9008, 86.9014, escIHWP);
      graphs[j][ind]->SetPoint(counts[j][ind], i+1, 100*escSign*fits[j]->GetParameter(0));
      graphs[j][ind]->SetPointError(counts[j][ind]++, 0.0, 100*fits[j]->GetParError(0));

      graphs[j+1][ind]->SetPoint(counts[j+1][ind]++, i+1, fits[j]->GetChisquare()*1.0/fits[j]->GetNDF());
      if(i+1 < xmins[ind]) xmins[ind] = i+1;
      if(i+1 > xmaxs[ind]) xmaxs[ind] = i+1;
      if(100*escSign*fits[j]->GetParameter(0) - fits[j]->GetParError(0) < ymins[j]) 
        ymins[j] = 100*escSign*fits[j]->GetParameter(0) - fits[j]->GetParError(0);
      if(100*escSign*fits[j]->GetParameter(0) + fits[j]->GetParError(0) > ymaxs[j]) 
        ymaxs[j] = 100*escSign*fits[j]->GetParameter(0) + fits[j]->GetParError(0);
      if(fits[j]->GetChisquare()*1.0/fits[j]->GetNDF() < ymins[j+1]) ymins[j+1] = fits[j]->GetChisquare()*1.0/fits[j]->GetNDF();
      if(fits[j]->GetChisquare()*1.0/fits[j]->GetNDF() > ymaxs[j+1]) ymaxs[j+1] = fits[j]->GetChisquare()*1.0/fits[j]->GetNDF();
    }
  }

  for(Int_t i = 0; i < nData; i++){
    cans[i]->cd();
    for(Int_t j = 0; j < nGraphs; j++){
      graphs[i][j]->GetXaxis()->SetLimits(0, nPairs+1);
      graphs[i][j]->GetXaxis()->SetRangeUser(0, nPairs+1);
      graphs[i][j]->GetYaxis()->SetLimits(0.95*ymins[i], 1.05*ymaxs[i]);
      graphs[i][j]->GetYaxis()->SetRangeUser(0.95*ymins[i], 1.05*ymaxs[i]);
      if(j == 0) graphs[i][j]->Draw("ap");
      else graphs[i][j]->Draw("p && same");
      graphs[i][j]->Fit(Form("escFit_%s_graph%i", names[i].Data(), j), "Q", "", xmins[j], xmaxs[j]);
      TText *line1 = texts[i][j]->AddText(Form("------%s------\n", getLegendEntry(j).Data()));
      line1->SetTextColor(getColor(j));
      TText *line2 = texts[i][j]->AddText(Form("Pol0: %.3f +/- %.3f\n", escFits[i][j]->GetParameter(0), escFits[i][j]->GetParError(0)));
      line2->SetTextColor(getColor(j));
      TText *line3 = texts[i][j]->AddText(Form("Chi2 / NDF: %.3f / %i\n", escFits[i][j]->GetChisquare(), escFits[i][j]->GetNDF()));
      line3->SetTextColor(getColor(j));
      texts[i][j]->Draw("same");
    }

    cans[i+1]->cd();
    for(Int_t j = 0; j < nGraphs; j++){
      graphs[i+1][j]->GetXaxis()->SetLimits(0, nPairs+1);
      graphs[i+1][j]->GetXaxis()->SetRangeUser(0, nPairs+1);
      graphs[i+1][j]->GetYaxis()->SetRangeUser(0.0, 1.1*ymaxs[i+1]);
      if(j == 0) graphs[i+1][j]->Draw("ap");
      else graphs[i+1][j]->Draw("p && same");
    }
  }
}
