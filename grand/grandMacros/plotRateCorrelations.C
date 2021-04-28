#include "../grandConstruction/plot.h"
#include "../grandOnline/makePlots.h"

using namespace std;

vector<vector<Int_t>> readKeysFile(TString expt){
  vector<vector<Int_t>> keys;
  ifstream keysfile(Form("%s/%s_rateRegions.key", getenv("COMPMON_MAPS"), expt.Data()));
  if(keysfile.is_open()){
    string line;
    while(getline(keysfile, line)){
      vector<Int_t> keyRange; stringstream ss(line);
      for(Int_t i; ss >> i;){
        keyRange.push_back(i);
        if(ss.peek() == ',')
          ss.ignore();
      }
      keys.push_back(keyRange);
    }
  }
  keysfile.close();

  return keys;
}

void plotRateCorrelations(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  vector<vector<Int_t>> ranges = readKeysFile(expt);
  const Int_t nRanges = (Int_t)ranges.size();

  TGraphErrors *graphs[nRanges];
  Int_t nPts[nRanges];
  TF1 *fits[nRanges];
  TH1F *rateHists[nRanges];
  TH1F *xPosHists[nRanges];
  TH1F *yPosHists[nRanges];
  TH1F *xDiffHists[nRanges];
  TH1F *yDiffHists[nRanges];
  TGraph *gAllRates = new TGraph();
  for(Int_t i = 0; i < nRanges; i++){
    graphs[i] = new TGraphErrors();
    graphs[i]->SetName(Form("g%i-%i", ranges[i][0], ranges[i][1]));
    nPts[i] = 0;
    fits[i] = new TF1(Form("f%i-%i", ranges[i][0], ranges[i][1]), "pol1");

    TString rateName = Form("rateHist_%i-%i", ranges[i][0], ranges[i][1]);
    TString xPosName = Form("xPosHist_%i-%i", ranges[i][0], ranges[i][1]);
    TString yPosName = Form("yPosHist_%i-%i", ranges[i][0], ranges[i][1]);
    TString xDiffName = Form("xDiffHist_%i-%i", ranges[i][0], ranges[i][1]);
    TString yDiffName = Form("yDiffHist_%i-%i", ranges[i][0], ranges[i][1]);
    TString cut = Form("CycleCut==0 && runNum>=%i && runNum<=%i", ranges[i][0], ranges[i][1]);

    cyc->Project(rateName.Data(), "120*(CentralRateLasOn - CentralRateLasOff)/(BeamCurrent*LaserPower)", cut.Data());
    cyc->Project(xPosName.Data(), "(bpmAx + bpmBx)/2.0", cut.Data());
    cyc->Project(yPosName.Data(), "(bpmAy + bpmBy)/2.0", cut.Data());
    cyc->Project(xDiffName.Data(), "(diff_bpmAx + diff_bpmBx)/2.0", cut.Data());
    cyc->Project(yDiffName.Data(), "(diff_bpmAy + diff_bpmBy)/2.0", cut.Data());
    rateHists[i] = (TH1F *)gDirectory->Get(rateName.Data());
    xPosHists[i] = (TH1F *)gDirectory->Get(xPosName.Data());
    yPosHists[i] = (TH1F *)gDirectory->Get(yPosName.Data());
    xDiffHists[i] = (TH1F *)gDirectory->Get(xDiffName.Data());
    yDiffHists[i] = (TH1F *)gDirectory->Get(yDiffName.Data());
  }
  

  Int_t runNum, cycleCut, sign;
  DataVar bpmAy, bpmBy, diffAx, diffBx, diffAy, diffBy, rateOn, rateOff, las, bcm;
  PolVar anPow;
  PolVar asym0;
  Int_t ptNum = 0;

  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("bpmAy", &bpmAy);
  cyc->SetBranchAddress("bpmBy", &bpmBy);
  cyc->SetBranchAddress("diff_bpmAx", &diffAx);
  cyc->SetBranchAddress("diff_bpmAy", &diffAy);
  cyc->SetBranchAddress("diff_bpmBx", &diffBx);
  cyc->SetBranchAddress("diff_bpmBy", &diffBy);
  cyc->SetBranchAddress("CentralRateLasOn", &rateOn);
  cyc->SetBranchAddress("CentralRateLasOff", &rateOff);
  cyc->SetBranchAddress("BeamCurrent", &bcm);
  cyc->SetBranchAddress("LaserPower", &las);
  cyc->SetBranchAddress("Asym0", &asym0);
  cyc->SetBranchAddress("CycleCut", &cycleCut);
  cyc->SetBranchAddress("AnalyzingPower", &anPow);
  cyc->SetBranchAddress("sign", &sign);

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(cycleCut != 0 or sign == 0) continue;
    Int_t curRange = -1;
    for(Int_t j = 0; j < nRanges; j++){
      if(runNum >= ranges[j][0] && runNum <= ranges[j][1]){
        curRange = j;
        break;
      }
    }
    if(curRange == -1) continue;

    Float_t asymMean = sign*asym0.mean;
    Float_t asymErr = asym0.meanErr;
    graphs[curRange]->SetPoint(nPts[curRange], sign*(diffAy.mean + diffBy.mean)/2.0, asymMean);
    graphs[curRange]->SetPointError(nPts[curRange], 0.0, asymErr);
    nPts[curRange] = nPts[curRange] + 1;
    gAllRates->SetPoint(ptNum++, runNum, 120*(rateOn.mean - rateOff.mean)/(bcm.mean*las.mean));
  }

  TGraphErrors *g = new TGraphErrors();
  TGraph *gPulls = new TGraph();
  TGraph *gChi = new TGraph();
  TH1F *hPulls = new TH1F("hPulls", "CREX Correlation Slope Pulls", 100, -8.0, 8.0);
  TGraphErrors *gRate = new TGraphErrors();
  TGraphErrors *gXPos = new TGraphErrors();
  TGraphErrors *gYPos = new TGraphErrors();
  TGraphErrors *gXDiff = new TGraphErrors();
  TGraphErrors *gYDiff = new TGraphErrors();

  for(Int_t i = 0; i < nRanges; i++){
    graphs[i]->Fit(fits[i], "Q", "");
    g->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, fits[i]->GetParameter(1));
    g->SetPointError(i, 0.0, fits[i]->GetParError(1));
    gPulls->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, TMath::Abs(fits[i]->GetParameter(1))*1.0/fits[i]->GetParError(1));
    hPulls->Fill(fits[i]->GetParameter(1)*1.0/fits[i]->GetParError(1));
    gChi->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, fits[i]->GetChisquare()*1.0/fits[i]->GetNDF());
    gRate->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, rateHists[i]->GetMean());
    gRate->SetPointError(i, 0.0, rateHists[i]->GetMeanError());
    gXPos->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, xPosHists[i]->GetMean());
    gXPos->SetPointError(i, 0.0, xPosHists[i]->GetMeanError());
    gYPos->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, yPosHists[i]->GetMean());
    gYPos->SetPointError(i, 0.0, yPosHists[i]->GetMeanError());
    gXDiff->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, xDiffHists[i]->GetMean());
    gXDiff->SetPointError(i, 0.0, xDiffHists[i]->GetMeanError());
    gYDiff->SetPoint(i, (ranges[i][0] + ranges[i][1])*1.0/2.0, yDiffHists[i]->GetMean());
    gYDiff->SetPointError(i, 0.0, yDiffHists[i]->GetMeanError());
    printf("%i,%i,%.4f,%.4f,%.4f\n", ranges[i][0], ranges[i][1], fits[i]->GetParameter(1), fits[i]->GetParError(1), fits[i]->GetParameter(1)*1.0/fits[i]->GetParError(1));
  }

  TCanvas *c = new TCanvas("c", "Range Fit Summary Canvas", 1400, 800);
  c->Divide(1, 2);
  c->cd(1)->SetGridx();
  c->cd(1)->SetGridy();
  g->SetTitle("CREX Correlation Slopes by Run Range");
  g->GetXaxis()->SetTitle("Median Run # of Range");
  g->GetYaxis()->SetTitle("Slope of Correlation Fit [pct um^-1]");
  g->SetMarkerStyle(8);
  g->Draw("ap");
  // c->cd(2)->SetGridx();
  // c->cd(2)->SetGridy();
  // gPulls->SetTitle("CREX Correlation Slope Pulls");
  // gPulls->GetXaxis()->SetTitle("Median Run # of Range");
  // gPulls->GetYaxis()->SetTitle("Correlation Slope Pull");
  // gPulls->SetMarkerStyle(8);
  // gPulls->Draw("ap");
  c->cd(2)->SetGridx();
  c->cd(2)->SetGridy();
  gChi->SetTitle("CREX Correlation Slope Chi2 / NDF");
  gChi->GetXaxis()->SetTitle("Median Run # of Range");
  gChi->GetYaxis()->SetTitle("Correlation Fit Chi2 / NDF");
  gChi->SetMarkerStyle(8);
  gChi->Draw("ap");

  TCanvas *cPulls = new TCanvas("cPulls", "Pulls Canvas", 1200, 800);
  cPulls->cd();
  hPulls->GetXaxis()->SetTitle("Correlation Slope Pull");
  hPulls->Draw();

  TCanvas *cRate = new TCanvas("cRate", "Rate Canvas", 1200, 800);
  cRate->cd();
  gRate->SetTitle("CREX Background Subtracted Central Rate");
  gRate->GetXaxis()->SetTitle("Median Run # of Range");
  gRate->GetYaxis()->SetTitle("Central Rate [Hz W^-1 uA^-1]");
  gRate->SetMarkerStyle(8);
  gRate->Draw("ap");

  TCanvas *cPos = new TCanvas("cPos", "Position Canvas", 1400, 800);
  cPos->Divide(1, 2);
  cPos->cd(1);
  gXPos->SetTitle("CREX Average X-Position");
  gXPos->GetXaxis()->SetTitle("Median Run # of Range");
  gXPos->GetYaxis()->SetTitle("Avg X Pos (mm)");
  gXPos->SetMarkerStyle(8);
  gXPos->Draw("ap");
  cPos->cd(2);
  gYPos->SetTitle("CREX Average Y-Position");
  gYPos->GetXaxis()->SetTitle("Median Run # of Range");
  gYPos->GetYaxis()->SetTitle("Avg Y Pos (mm)");
  gYPos->SetMarkerStyle(8);
  gYPos->Draw("ap");

  TCanvas *cAllRates = new TCanvas("cAllRates", "All Rates Canvas", 1200, 800);
  cAllRates->cd();
  gAllRates->SetTitle("CREX All Background Subtracted Rates");
  gAllRates->GetXaxis()->SetTitle("runNum");
  gAllRates->GetYaxis()->SetTitle("Central Rate [Hz W^-1 uA^-1]");
  gAllRates->GetYaxis()->SetLimits(0.0, 1.0);
  gAllRates->GetYaxis()->SetRangeUser(0.0, 1.0);
  gAllRates->Draw("ap");
  for(Int_t i = 0; i < nRanges; i++){
    TLine *lStart = new TLine(ranges[i][0], 0.0, ranges[i][0], 1.0);
    TLine *lEnd = new TLine(ranges[i][1], 0.0, ranges[i][1], 1.0);
    lStart->SetLineColor(kRed);
    lEnd->SetLineColor(kRed);
    lStart->Draw("same");
    // lEnd->Draw("same");
  }

  TCanvas *cDiff = new TCanvas("cDiff", "BPM Diff Canvas", 1400, 800);
  cDiff->Divide(1, 2);
  cDiff->cd(1);
  gXDiff->SetTitle("CREX Average X Difference");
  gXDiff->GetXaxis()->SetTitle("Median Run # of Range");
  gXDiff->GetYaxis()->SetTitle("Avg X Difference [mm]");
  gXDiff->SetMarkerStyle(8);
  gXDiff->Draw("ap");
  cDiff->cd(2);
  gYDiff->SetTitle("CREX Average Y Difference");
  gYDiff->GetXaxis()->SetTitle("Median Run # of Range");
  gYDiff->GetYaxis()->SetTitle("Avg Y Difference [mm]");
  gYDiff->SetMarkerStyle(8);
  gYDiff->Draw("ap");
}
