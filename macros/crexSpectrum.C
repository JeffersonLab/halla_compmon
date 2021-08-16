#include "../online/utils.h"
#include "simulationHists.h"

using namespace std;


Float_t findComptonEdge(Int_t run, TH1F *hSum, Int_t nBins){
  printf("  Finding Compton Edge...\n");
  Float_t max = -1e16;
  Int_t maxBin = 0;
  for(Int_t bin = 1; bin <= nBins; bin++){
    Float_t center = hSum->GetBinCenter(bin);
    if(center < 10000) continue;
    Bool_t replaceMax = (hSum->GetBinContent(bin) > max);
    max = replaceMax ? hSum->GetBinContent(bin) : max;
    maxBin = replaceMax ? bin : maxBin;
  }

  printf("  Maximum and bin found! Maximum: %.2f, Max Bin: %i Max SRAU: %.2f\n", max, maxBin, hSum->GetBinCenter(maxBin));

  Float_t ce;
  for(Int_t bin = maxBin; bin <= nBins; bin++){
    if(hSum->GetBinContent(bin) < 0.5*max){
      ce = hSum->GetBinCenter(bin);
      break;
    }
  }

  printf("  Found compton edge: %.2f\n", ce);

  return ce;
}

Float_t getClockCutoff(TTree *trig, Int_t run){
  printf("  Finding sumClock cutoff...\n");
  const Int_t nHists = 16;
  TH1F *hists[nHists];
  Float_t refCount = 0;
  Float_t limit = 0;
  for(Int_t i = 0; i < nHists; i++){
    printf("    Trying clock value: %f...\n", i*100e3);
    TString clockCut = Form("sumIsRandom==0 && sum>0 && (laserState==0 || laserState==1) && beamState==1 && abs(3790 - sumPre)<10 && sumClock>=%f && sumClock<%f", i*100e3, (i+1)*100e3);
    trig->Project(Form("hClock%i_run%i", i, run), "sum", clockCut.Data());
    hists[i] = (TH1F *)gDirectory->Get(Form("hClock%i_run%i", i, run));
    if(i == 0){
      refCount = hists[i]->GetEntries();
    }
    limit = i*100e3;
    if(hists[i]->GetEntries() < 0.5*refCount){
      break;
    }
  }
  
  printf("  Found sumClock upper bound: %f\n", limit);
  return limit;
  // return 1.6e16;
}

Float_t getNormFactor(TH1F *hSpec, Int_t nBins){
  Float_t max = -1e16;
  Int_t maxBin = 0;
  for(Int_t bin = 1; bin <= nBins; bin++){
    Float_t center = hSpec->GetBinCenter(bin);
    if(center < 0.5) continue;
    Bool_t replaceMax = (hSpec->GetBinContent(bin) > max);
    max = replaceMax ? hSpec->GetBinContent(bin) : max;
    maxBin = replaceMax ? bin : maxBin;
  }

  return max;
}

Float_t getSubtractionFactor(TTree *trig, Int_t run, Float_t ce, Float_t clock, Float_t extraFactor){
  printf("  Finding subtraction factor...\n");
  Float_t energyLo = 2.4*ce;
  Float_t energyHi = 90e3;
  if(ce > 30e3){
    energyLo = 1.5*ce;
  }
  TString cuts = Form("beamState==1 && sumIsRandom==0 && abs(sumPre - 3790)<10 && sum>%f && sum<%f", energyLo, energyHi);
  TString lasOn = Form("(laserState==0 || laserState==1) && sumClock<%f", clock);
  TString lasOff = "(laserState==2 || laserState==3) && sumClock<1.6e6";
  
  TString hName1("hBackOn_run%i");
  TString hName2("hBackOff_run%i");
  TH1F *hBackOn = new TH1F(hName1.Data(), "", 100, energyLo, energyHi);
  TH1F *hBackOff = new TH1F(hName2.Data(), "", 100, energyLo, energyHi);
  trig->Project(hName1.Data(), "sum", Form("%s && %s", cuts.Data(), lasOn.Data()));
  trig->Project(hName2.Data(), "sum", Form("%s && %s", cuts.Data(), lasOff.Data()));
  Float_t fac = hBackOn->Integral()*1.0/hBackOff->Integral();
  hBackOn->SetLineColor(kGreen + 2);
  hBackOff->SetLineColor(kRed);
  hBackOff->Scale(fac);

  TCanvas *cSubNorm = new TCanvas(Form("cSubNorm_run%i_fac%.2f", run, extraFactor), Form("Run %i Subtraction Normalization Canvas", run), 1200, 800);
  cSubNorm->cd();
  THStack *hsBackSub = new THStack(Form("hsBackSub_run%i_fac%.2f", run, extraFactor), Form("Run %i: Normalization Range Comparison ;Pulse Intergal [SRAU];", run));
  hsBackSub->Add(hBackOn);
  hsBackSub->Add(hBackOff);
  hsBackSub->Draw("nostack");

  TLegend *leg = new TLegend(0.70, 0.80, 0.90, 0.90);
  leg->AddEntry(hBackOn, Form("Las On, %.0f Entries", hBackOn->GetEntries()));
  leg->AddEntry(hBackOff, Form("Las Off, %.0f Entries", hBackOff->GetEntries()));
  leg->Draw("same");

  printf("  Found subtraction factor: %.4f\n", fac);
  return extraFactor*fac;
}

vector<TH1F *> getHistograms(Int_t nBins){
  TString hTitle = "Spectrum Comparison";
  const Int_t nHists = 3;
  TString histLabels[nHists] = {"OnSub", "OnNoSub", "Off"};
  Int_t colors[nHists] = {kRed, kBlue, kRed};
  TString axisTitle("Pulse Intergal [% of CE]");
  vector<TH1F *> hists;
  for(Int_t i = 0; i < nHists; i++){
    TH1F *hSpec = new TH1F(Form("hSpec%s", histLabels[i].Data()), hTitle.Data(), nBins, 0, 2);
    // hSpec->SetLineColor(colors[i]);
    hSpec->GetXaxis()->SetTitle(axisTitle.Data());
    hists.push_back(hSpec);
  }
  
  return hists;
}

vector<TH1F *> fillHistogram(TTree *trig, vector<TH1F *> hists, Int_t run, Float_t ce, Float_t clock){
  Float_t sum, sumPre, sumClock;
  Int_t beamState, laserState;
  Bool_t random;
  trig->SetBranchAddress("sum", &sum);
  trig->SetBranchAddress("sumPre", &sumPre);
  trig->SetBranchAddress("sumIsRandom", &random);
  trig->SetBranchAddress("beamState", &beamState);
  trig->SetBranchAddress("laserState", &laserState);
  trig->SetBranchAddress("sumClock", &sumClock);

  printf("  Plotting run %i...\n", run);

  Float_t pct = 0.1;
  for(Int_t i = 0; i < trig->GetEntries(); i++){
    trig->GetEntry(i);
    if(i > pct*trig->GetEntries()){
      printf("    %.2f%% done with trig...\n", 100*pct);
      pct += 0.1;
    }
    if(!(beamState==1 && laserState<4 && !random && abs(sumPre - 3790)<10 && sum<75e3 && sumClock<1.6e6))
      continue;
    if((laserState == 0 || laserState == 1) && sumClock<clock){
      hists[0]->Fill(sum*1.0/ce);
      if(sum>3e3){
        hists[1]->Fill(sum*1.0/ce);
      }
    }
    else if(laserState == 2 || laserState == 3){
      hists[2]->Fill(sum*1.0/ce);
    }
  }

  printf("  Plotted run %i!\n", run);
  return hists;
}


// *********************************************
// A one-run variant of crexSpectrumOverlay.C  * 
// *********************************************
void crexSpectrum(Int_t runNum){
  const Int_t nRuns = 4;
  Int_t runs[nRuns] = {runNum, runNum, runNum, 5718};
  Float_t bkSubFactor[nRuns] = {1.0, 0.5, 2.0, 1.0};
  Int_t color[nRuns] = {kBlue, kRed, kGreen + 2, kViolet};
  TString legTitles[nRuns] = {Form("Run %i (Reg. Subtraction)", runNum), Form("Run %i (2x Undersubtraction)", runNum), 
                              Form("Run %i (2x Oversubtraction)", runNum), "Run 5718 (5.88 mm offset)"};
  TLegend *leg = new TLegend(0.70, 0.70, 0.90, 0.90);
  TLegend *legSub = new TLegend(0.70, 0.70, 0.90, 0.90);
  vector<vector<TH1F *>> allHists;
  THStack *hs = new THStack("hs", "Spectrum Comparison (Background Subtracted) ;Pulse Intergal [% of CE];");
  THStack *hsSub = new THStack("hsSub", "Background Subtraction Amount ;Pulse Intergal [% of CE]; Size of Subtraction Relative to CE bin size");

  Float_t regScaleFactor = 0.0;
  for(Int_t i = 0; i < nRuns; i++){
    printf("Starting on run %i...\n", runs[i]);
    TFile *file = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runs[i]));
    TTree *trig = (TTree *)file->Get("triggerwise");
    Int_t nBins = 500;
    Float_t clock = getClockCutoff(trig, runs[i]);
    TString cuts = Form("beamState==1 && (laserState==0 || laserState==1) && sumIsRandom==0 && abs(sumPre - 3790)<10 && sum>4e3 && sum<75e3 && sumClock<%f", clock);
    TH1F *hSum = new TH1F(Form("hSum_run%i", runs[i]), Form("Run %i Spectrum", runs[i]), nBins, 0, 75e3);
    trig->Project(Form("hSum_run%i", runs[i]), "sum", cuts.Data(), "goff");

    Float_t ce = findComptonEdge(runs[i], hSum, nBins);
    Float_t fac = getSubtractionFactor(trig, runs[i], ce, clock, bkSubFactor[i]);

    nBins = 1000;
    vector<TH1F *> hists = getHistograms(nBins);
    hists = fillHistogram(trig, hists, runs[i], ce, clock);
    hists[0]->SetLineColor(color[i]);
    hists[2]->SetLineColor(color[i]);
    hists[0]->Add(hists[2], -fac);

    Float_t scale = getNormFactor(hists[0], nBins);
    if(i == 0){regScaleFactor = scale;}
    Float_t curScale = (i == 0 || i == 3) ? scale : regScaleFactor;
    hists[0]->Scale(1.0/curScale);
    hists[1]->Scale(1.0/curScale);
    hists[2]->Scale(fac/curScale);
    
    leg->AddEntry(hists[0], legTitles[i]);
    legSub->AddEntry(hists[2], legTitles[i]);
    // leg->AddEntry(hists[1], Form("Run %i No Sub", runs[i]));
    // TLegend *legSub = new TLegend(0.7, 0.75, 0.9, 0.9);
    // legSub->AddEntry(hists[0][2], Form("Run %i", run1));
    // legSub->AddEntry(hists[1][2], Form("Run %i", run2));
    allHists.push_back(hists);
    hs->Add(hists[0]);
    // hs->Add(hists[1]);
    hsSub->Add(hists[2]);
    printf("Finished run %i!\n", runs[i]);
    // break;
  }

  // hs->GetXaxis()->SetTitle("Pulse Intergal [% of CE]");
  fillSimHist_6mm(hs, leg, kBlack);

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", "Spectrum Comparison", 1200, 800);
  c->Divide(1, 2);
  c->cd(1)->SetGridx();
  c->cd(1)->SetGridy();
  hs->Draw("nostack");
  hs->GetYaxis()->SetRangeUser(0, 2);
  hs->Draw("nostack");
  leg->Draw("same");

  //TCanvas *cSub = new TCanvas("cSub", "Spectrum Subtraction Amount", 1200, 800);
  c->cd(2)->SetGridx();
  c->cd(2)->SetGridy();
  c->cd(2)->SetLogy();
  hsSub->Draw("nostack");
  legSub->Draw("same");
}
