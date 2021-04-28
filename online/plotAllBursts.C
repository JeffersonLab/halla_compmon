#include "utils.h"

using namespace std;

void plotAllBursts_alt(Int_t runNum){
  vector<vector<int>> burstList = getCycleList(runNum, true);
  vector<TChain *> chains = loadChain(runNum);
  TTree *quartetwise = chains[1];
  TString posH("PosHelAcc0/PosHelNSamples0");
  TString negH("NegHelAcc0/NegHelNSamples0");
  
  TCanvas *c = new TCanvas("c", "Burst Plots Canvas", 1200, 800);
  c->Divide(2, 2);

  Int_t burstCount = 0; Int_t currentState = 0;
  printf("Burst list size: %i\n", (Int_t)burstList.size());
  for(Int_t i = 0; i < burstList.size(); i++){
    if(burstList[i][0] % 2 == 1) continue;
    TString laserCut("(laserState==2 || laserState==3)");
    Int_t lineColor = 0;
    if(burstList[i][1] == 1){
      lineColor = kRed;
    }
    else if(burstList[i][1] == 2){
      laserCut = "(laserState==0 || laserState==1)";
      lineColor = kGreen + 2;
    }
    else{
      lineColor = kRed + 2;
    }

    if(currentState != burstList[i][1]){
      burstCount = 1; currentState = burstList[i][1];
    }
    else{
      burstCount++;
    }
    TString cuts = Form("%s && beamState==1 && dithering==0 && firstMPSnumber>=%i && firstMPSnumber<=%i",
                        laserCut.Data(), burstList[i][2], burstList[i][3]);
    TString hPosName = Form("hPos_cycle%i_period%i_burst%i", burstList[i][0], burstList[i][1], burstCount);
    TString hNegName = Form("hNeg_cycle%i_period%i_burst%i", burstList[i][0], burstList[i][1], burstCount);
    TString hDiffName = Form("hDiff_cycle%i_period%i_burst%i", burstList[i][0], burstList[i][1], burstCount);
    TString hSummName = Form("hSumm_cycle%i_period%i_burst%i", burstList[i][0], burstList[i][1], burstCount);
    quartetwise->Project(hPosName.Data(), posH.Data(), cuts.Data(), "goff");
    quartetwise->Project(hNegName.Data(), negH.Data(), cuts.Data(), "goff");
    quartetwise->Project(hDiffName.Data(), Form("%s - %s", posH.Data(), negH.Data()), cuts.Data(), "goff");
    quartetwise->Project(hSummName.Data(), Form("%s + %s", posH.Data(), negH.Data()), cuts.Data(), "goff");

    TH1F *hPos = (TH1F *)gDirectory->Get(hPosName.Data());
    TH1F *hNeg = (TH1F *)gDirectory->Get(hNegName.Data());
    TH1F *hDiff = (TH1F *)gDirectory->Get(hDiffName.Data());
    TH1F *hSumm = (TH1F *)gDirectory->Get(hSummName.Data());

    hPos->SetTitle(Form("Burst %i.%i.%i.%i: PosHelAcc0", runNum, burstList[i][0], burstList[i][1], burstCount));
    hNeg->SetTitle(Form("Burst %i.%i.%i.%i: NegHelAcc0", runNum, burstList[i][0], burstList[i][1], burstCount));
    hDiff->SetTitle(Form("Burst %i.%i.%i.%i: Diff Acc0", runNum, burstList[i][0], burstList[i][1], burstCount));
    hSumm->SetTitle(Form("Burst %i.%i.%i.%i: Sum Acc0",  runNum, burstList[i][0], burstList[i][1], burstCount));
    hPos->SetLineColor(lineColor); hDiff->SetLineColor(lineColor);
    hNeg->SetLineColor(lineColor); hSumm->SetLineColor(lineColor);

    c->cd(1); hPos->Draw();
    c->cd(2); hNeg->Draw();
    c->cd(3); hDiff->Draw();
    c->cd(4); hSumm->Draw();

    c->Print(Form("%s/runs/Run%i/cycle%02i_period%i_burst%02i.pdf", getenv("COMPMON_WEB"), 
                  runNum, burstList[i][0], burstList[i][1], burstCount), "pdf");
  }

  gSystem->Exec(Form("pdfunite %s/runs/Run%i/cycle*_period*_burst*.pdf %s/runs/Run%i/burst_diagnostics.pdf", 
                getenv("COMPMON_WEB"), runNum, getenv("COMPMON_WEB"), runNum));
  gSystem->Exec(Form("rm -rf %s/runs/Run%i/cycle*_period*_burst*.pdf", getenv("COMPMON_WEB"), runNum));
}

vector<vector<vector<vector<Int_t>>>> sortBurstsByCycle(vector<vector<int>> allBursts){
  vector<vector<vector<vector<Int_t>>>> sortedBursts;
  if(allBursts.size() == 0) return sortedBursts;
  
  Int_t curCycle = 0; Int_t curPeriod = 0;
  vector<vector<vector<Int_t>>> cycleBursts;
  vector<vector<Int_t>> periodBursts;
  for(Int_t i = 0; i < allBursts.size(); i++){
    vector<Int_t> burstLimits;
    if(i == 0){
      curCycle = allBursts[i][0]; curPeriod = allBursts[i][1];
    }
    if(curPeriod != allBursts[i][1]){
      cycleBursts.push_back(periodBursts);
      periodBursts.clear(); curPeriod = allBursts[i][1];
    }
    if(curCycle != allBursts[i][0]){
      sortedBursts.push_back(cycleBursts);
      cycleBursts.clear(); curCycle = allBursts[i][0];
    }

    burstLimits.push_back(allBursts[i][2]); burstLimits.push_back(allBursts[i][3]);
    periodBursts.push_back(burstLimits);
  }
  cycleBursts.push_back(periodBursts);
  sortedBursts.push_back(cycleBursts);

  return sortedBursts;
}

void plotAllBursts(Int_t runNum){
  vector<vector<int>> burstList = getCycleList(runNum, true);
  vector<TChain *> chains = loadChain(runNum);
  TTree *quartetwise = chains[1];
  TString posH("PosHelAcc0/PosHelNSamples0");
  TString negH("NegHelAcc0/NegHelNSamples0");
  TString diff(Form("%s - %s", posH.Data(), negH.Data()));
  TString summ(Form("%s + %s", posH.Data(), negH.Data()));

  const Int_t nPlots = 3;
  const Int_t nHists = 2;
  TString msmts[nPlots] = {"Diffs", "Sums", "Asyms"};
  TString codes[nPlots] = {diff, summ, ""};
  TString states[nHists] = {"On", "Off"};

  vector<vector<vector<vector<Int_t>>>> runBursts = sortBurstsByCycle(burstList);
  TH1F *hists[nPlots][nHists];
  TF1 *fits[nPlots][nHists];
  TPaveText *stats[nPlots][nHists];
  Int_t colors[nHists] = {kGreen + 2, kRed};
  Int_t marker[nHists] = {20, 21};

  for(Int_t i = 0; i < runBursts.size(); i++){
    printf("Starting cycle %i...\n", i+1);
    Int_t per1Bursts = (Int_t)runBursts[i][0].size();
    Int_t per2Bursts = (Int_t)runBursts[i][1].size();
    Int_t per3Bursts = (Int_t)runBursts[i][2].size();
    Int_t cycleBursts = per1Bursts + per2Bursts + per3Bursts;

    Float_t ymins[nPlots];
    Float_t ymaxs[nPlots];

    for(Int_t plot = 0; plot < nPlots; plot++){
      for(Int_t state = 0; state < nHists; state++){
        TString hName = Form("h%s%s_cycle%i", msmts[plot].Data(), states[state].Data(), i+1);
        TString fName = Form("f%s%s_cycle%i", msmts[plot].Data(), states[state].Data(), i+1);
        TString hTitle = Form("Cycle %i.%i: Multiplet %s", runNum, i+1, msmts[plot].Data());

        hists[plot][state] = new TH1F(hName.Data(), hTitle.Data(), cycleBursts, 0, cycleBursts);
        hists[plot][state]->SetMarkerStyle(marker[state]); hists[plot][state]->SetMarkerColor(colors[state]);
        hists[plot][state]->SetLineColor(colors[state]); hists[plot][state]->SetStats(0);
        
        fits[plot][state] = new TF1(fName.Data(), "pol0");
        fits[plot][state]->SetLineColor(colors[state]);

        Float_t xmin = (state == 0) ? 0.1 : 0.7;
        Float_t xmax = (state == 0) ? 0.3 : 0.9;
        stats[plot][state] = new TPaveText(xmin, 0.75, xmax, 0.92, "blNDC");
        stats[plot][state]->SetFillColor(0); stats[plot][state]->SetBorderSize(1);
      }
      ymins[plot] = 1e16;
      ymaxs[plot] = -1e16;
    }
    
    for(Int_t j = 0; j < runBursts[i].size(); j++){
      printf("  Starting period %i sums and diffs...\n", j+1);
      for(Int_t k = 0; k < runBursts[i][j].size(); k++){
        Int_t laserInd = (j == 1) ? 0 : 1;
        TString laserCut = (j == 1) ? "(laserState==0 || laserState==1)" : "(laserState==2 || laserState==3)";
        Int_t base1 = (j > 0) ? per1Bursts : 0;
        Int_t base2 = (j > 1) ? per2Bursts : 0;
        TString cuts = Form("%s && beamState==1 && dithering==0 && firstMPSnumber>=%i && firstMPSnumber<=%i",
                            laserCut.Data(), runBursts[i][j][k][0], runBursts[i][j][k][1]);
        for(Int_t l = 0; l < nPlots - 1; l++){
          TString hName = Form("h%s%s_cycle%i_period%i_burst%i", msmts[l].Data(), states[laserInd].Data(), i+1, j+1, k+1);
          quartetwise->Project(hName.Data(), codes[l].Data(), cuts.Data(), "goff");
          TH1F *hData = (TH1F *)gDirectory->Get(hName.Data());
          hists[l][laserInd]->SetBinContent(k + base1 + base2 + 1, hData->GetMean());
          hists[l][laserInd]->SetBinError(k + base1 + base2 + 1, hData->GetMeanError());
          ymins[l] = (hData->GetMean() - hData->GetMeanError() < ymins[l]) ? hData->GetMean() - hData->GetMeanError() : ymins[l];
          ymaxs[l] = (hData->GetMean() + hData->GetMeanError() > ymaxs[l]) ? hData->GetMean() + hData->GetMeanError() : ymaxs[l];
        }
      }
    }

    for(Int_t plot = 0; plot < nPlots - 1; plot++){
      for(Int_t state = 0; state < nHists; state++){
        hists[plot][state]->Fit(fits[plot][state], "Q", "goff", 0, cycleBursts);
        stats[plot][state]->AddText(Form("========%s Stats========", states[state].Data()));
        stats[plot][state]->AddText(Form("Mean: %.4f +/- %.4f", fits[plot][state]->GetParameter(0), fits[plot][state]->GetParError(0)));
        stats[plot][state]->AddText(Form("Chi2 / NDF: %.4f / %i", fits[plot][state]->GetChisquare(), fits[plot][state]->GetNDF()));
      }
    }

    Float_t sumOnMean = fits[1][0]->GetParameter(0); Float_t sumOnMeanErr = fits[1][0]->GetParError(0);
    Float_t sumOffMean = fits[1][1]->GetParameter(0); Float_t sumOffMeanErr = fits[1][1]->GetParError(0);
    TString onAsym = Form("(%s)/(%s - %f)", diff.Data(), summ.Data(), sumOffMean);
    TString offAsym = Form("(%s)/(%f - %f)", diff.Data(), sumOnMean, sumOffMean);
    TString asyms[nHists] = {onAsym, offAsym};

    for(Int_t j = 0; j < runBursts[i].size(); j++){
      for(Int_t k = 0; k < runBursts[i][j].size(); k++){
        Int_t laserInd = (j == 1) ? 0 : 1;
        TString laserCut = (j == 1) ? "(laserState==0 || laserState==1)" : "(laserState==2 || laserState==3)";
        Int_t base1 = (j > 0) ? per1Bursts : 0;
        Int_t base2 = (j > 1) ? per2Bursts : 0;
        TString cuts = Form("%s && beamState==1 && dithering==0 && firstMPSnumber>=%i && firstMPSnumber<=%i",
                            laserCut.Data(), runBursts[i][j][k][0], runBursts[i][j][k][1]);
        TString hName = Form("h%s%s_cycle%i_period%i_burst%i", msmts[2].Data(), states[laserInd].Data(), i+1, j+1, k+1);
        quartetwise->Project(hName.Data(), asyms[laserInd], cuts.Data(), "goff");
        TH1F *hData = (TH1F *)gDirectory->Get(hName.Data());
        hists[2][laserInd]->SetBinContent(k + base1 + base2 + 1, hData->GetMean());
        hists[2][laserInd]->SetBinError(k + base1 + base2 + 1, hData->GetMeanError());
        ymins[2] = (hData->GetMean() - hData->GetMeanError() < ymins[2]) ? hData->GetMean() - hData->GetMeanError() : ymins[2];
        ymaxs[2] = (hData->GetMean() + hData->GetMeanError() > ymaxs[2]) ? hData->GetMean() + hData->GetMeanError() : ymaxs[2];
      }
    }

    for(Int_t state = 0; state < nHists; state++){
      hists[2][state]->Fit(fits[2][state], "Q", "goff", 0, cycleBursts);
      stats[2][state]->AddText(Form("========%s Stats========", states[state].Data()));
      stats[2][state]->AddText(Form("Mean: %.4f +/- %.4f", fits[2][state]->GetParameter(0), fits[2][state]->GetParError(0)));
      stats[2][state]->AddText(Form("Chi2 / NDF: %.4f / %i", fits[2][state]->GetChisquare(), fits[2][state]->GetNDF()));
    }

    TCanvas *can = new TCanvas(Form("can_cycle%i", i+1), "Burst Plots Canvas", 1200, 800);
    can->Divide(2, 2);

    for(Int_t plot = 0; plot < nPlots; plot++){
      //TString hTitle = Form("Cycle %i.%i: Multiplet %s", runNum, i+1, msmts[plot].Data());
      //THStack *hs = new THStack(); hs->SetTitle(hTitle.Data());
      //for(Int_t state = 0; state < nHists; state++){hs->Add(hists[plot][state]);}
      can->cd(plot + 1);
      //hs->Draw();
      for(Int_t state = 0; state < nHists; state++){
        hists[plot][state]->GetYaxis()->SetLimits(ymins[plot], ymaxs[plot]);
        hists[plot][state]->GetYaxis()->SetRangeUser(ymins[plot], ymaxs[plot]);
        if(state == 0) hists[plot][state]->Draw();
        else hists[plot][state]->Draw("same");
        stats[plot][state]->Draw("same");
      }
    }

    can->cd(0);
    can->Print(Form("%s/runs/Run%i/cycleBursts_%02i.pdf", getenv("COMPMON_WEB"), runNum, i+1), "pdf");
  }

  gSystem->Exec(Form("pdfunite %s/runs/Run%i/cycleBursts_*.pdf %s/runs/Run%i/burst_qVars.pdf", 
                getenv("COMPMON_WEB"), runNum, getenv("COMPMON_WEB"), runNum));
  gSystem->Exec(Form("rm -rf %s/runs/Run%i/cycleBursts_*.pdf", getenv("COMPMON_WEB"), runNum));
}
