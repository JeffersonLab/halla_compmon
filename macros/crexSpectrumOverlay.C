#include "../online/utils.h"

using namespace std;

Float_t findComptonEdge(Int_t run, TH1F *hSum, Int_t nBins){
  printf("Finding Compton Edge...\n");
  Float_t max = -1e16;
  Int_t maxBin = 0;
  for(Int_t bin = 1; bin <= nBins; bin++){
    Float_t center = hSum->GetBinCenter(bin);
    if(center < 10000) continue;
    Bool_t replaceMax = (hSum->GetBinContent(bin) > max);
    max = replaceMax ? hSum->GetBinContent(bin) : max;
    maxBin = replaceMax ? bin : maxBin;
  }

  printf("Maximum and bin found! Maximum: %.2f, Max Bin: %i Max SRAU: %.2f\n", max, maxBin, hSum->GetBinCenter(maxBin));

  Float_t ce;
  for(Int_t bin = maxBin; bin <= nBins; bin++){
    if(hSum->GetBinContent(bin) < 0.5*max){
      ce = hSum->GetBinCenter(bin);
      break;
    }
  }

  printf("Found compton edge: %.2f\n", ce);

  return ce;
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

Float_t getSubtractionFactor(TTree *trig, Float_t ce){
  printf("Finding subtraction factor...\n");
  TString cuts = Form("beamState==1 && sumIsRandom==0 && abs(sumPre - 3790)<10 && sum>%f && sum<%f", 2.4*ce, 90e3);
  TString lasOn = "(laserState==0 || laserState==1) && sumClock<1.2e6";
  TString lasOff = "(laserState==2 || laserState==3) && sumClock<1.6e6";
  
  TString hName1("hBackOn");
  TString hName2("hBackOff");
  trig->Project(hName1.Data(), "sum", Form("%s && %s", cuts.Data(), lasOn.Data()));
  trig->Project(hName2.Data(), "sum", Form("%s && %s", cuts.Data(), lasOff.Data()));
  TH1F *hBackOn = (TH1F *)gDirectory->Get(hName1.Data());
  TH1F *hBackOff = (TH1F *)gDirectory->Get(hName2.Data());
  Float_t fac = hBackOn->Integral()*1.0/hBackOff->Integral();

  printf("Found subtraction factor: %.4f\n", fac);
  return hBackOn->Integral()*1.0/hBackOff->Integral();
}

vector<vector<TH1F *>> getHistograms(Int_t nBins){
  TString hTitle = "Spectrum Comparison";
  const Int_t nRuns = 2;
  const Int_t nHists = 3;
  TString histLabels[nHists] = {"On", "Off"};
  Int_t colors[nRuns][nHists] = {{kRed, kGreen + 2, kRed}, {kBlue, kViolet, kBlue}};
  TString axisTitle("Pulse Intergal [% of CE]");
  vector<vector<TH1F *>> hists;
  for(Int_t i = 0; i < nRuns; i++){
    vector<TH1F *> runHists;
    for(Int_t j = 0; j < nHists; j++){
      TH1F *hSpec = new TH1F(Form("hSpec%s%i", histLabels[j].Data(), i+1), hTitle.Data(), nBins, 0, 2);
      hSpec->SetLineColor(colors[i][j]);
      hSpec->GetXaxis()->SetTitle(axisTitle.Data());
      runHists.push_back(hSpec);
    }
    hists.push_back(runHists);
  }
  
  return hists;
}

vector<vector<TH1F *>> fillHistogram(TTree *trig, vector<vector<TH1F *>> hists, Int_t ind, Int_t run, Float_t ce){
  Float_t sum, sumPre;
  Int_t beamState, laserState;
  Bool_t random;
  trig->SetBranchAddress("sum", &sum);
  trig->SetBranchAddress("sumPre", &sumPre);
  trig->SetBranchAddress("sumIsRandom", &random);
  trig->SetBranchAddress("beamState", &beamState);
  trig->SetBranchAddress("laserState", &laserState);

  printf("Plotting run %i...\n", run);

  Float_t pct = 0.1;
  for(Int_t i = 0; i < trig->GetEntries(); i++){
    trig->GetEntry(i);
    if(i > pct*trig->GetEntries()){
      printf("  %.2f%% done with trig...\n", 100*pct);
      pct += 0.1;
    }
    if(!(beamState==1 && laserState<4 && !random && abs(sumPre - 3790)<10 && sum<75e3))
      continue;
    if(laserState == 0 || laserState == 1){
      hists[ind][0]->Fill(sum*1.0/ce);
      if(sum>3e3){
        hists[ind][1]->Fill(sum*1.0/ce);
      }
    }
    else if(laserState == 2 || laserState == 3){
      hists[ind][2]->Fill(sum*1.0/ce);
    }
  }

  printf("Plotted run %i!\n", run);
  return hists;
}

void crexSpectrumOverlay(Int_t run1, Int_t run2){
  TFile *file1 = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), run1));
  TFile *file2 = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), run2));
  TTree *trig1 = (TTree *)file1->Get("triggerwise");
  TTree *trig2 = (TTree *)file2->Get("triggerwise");

  Int_t nBins = 500;
  TString cuts("beamState==1 && (laserState==0 || laserState==1) && sumIsRandom==0 && abs(sumPre - 3790)<10 && sum>4e3 && sum<75e3");
  TH1F *hSum1 = new TH1F("hSum1", Form("Run %i Spectrum", run1), nBins, 0, 75e3);
  TH1F *hSum2 = new TH1F("hSum2", Form("Run %i Spectrum", run2), nBins, 0, 75e3);
  trig1->Project("hSum1", "sum", cuts.Data(), "goff");
  trig2->Project("hSum2", "sum", cuts.Data(), "goff");

  Float_t ce1 = findComptonEdge(run1, hSum1, nBins);
  Float_t ce2 = findComptonEdge(run2, hSum2, nBins);

  Float_t fac1 = getSubtractionFactor(trig1, ce1);
  Float_t fac2 = getSubtractionFactor(trig2, ce2);

  nBins = 250;
  vector<vector<TH1F *>> hists = getHistograms(nBins);
  hists = fillHistogram(trig1, hists, 0, run1, ce1);
  hists = fillHistogram(trig2, hists, 1, run2, ce2);
  hists[0][0]->Add(hists[0][2], -fac1);
  hists[1][0]->Add(hists[1][2], -fac2);

  Float_t scale1 = getNormFactor(hists[0][0], nBins);
  Float_t scale2 = getNormFactor(hists[1][0], nBins);
  hists[0][0]->Scale(1.0/scale1);
  hists[1][0]->Scale(1.0/scale2);
  hists[0][1]->Scale(1.0/scale1);
  hists[1][1]->Scale(1.0/scale2);
  hists[0][2]->Scale(fac1/scale1);
  hists[1][2]->Scale(fac2/scale2);
  TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9);
  leg->AddEntry(hists[0][0], Form("Run %i Bk Sub", run1));
  leg->AddEntry(hists[1][0], Form("Run %i Bk Sub", run2));
  leg->AddEntry(hists[0][1], Form("Run %i No Sub", run1));
  leg->AddEntry(hists[1][1], Form("Run %i No Sub", run2));
  TLegend *legSub = new TLegend(0.7, 0.75, 0.9, 0.9);
  legSub->AddEntry(hists[0][2], Form("Run %i", run1));
  legSub->AddEntry(hists[1][2], Form("Run %i", run2));

  THStack *hs = new THStack("hs", "Spectrum Comparison ;Pulse Intergal [% of CE];");
  THStack *hsSub = new THStack("hsSub", "Background Subtraction Amount ;Pulse Intergal [% of CE]; Size of Subtraction Relative to CE bin size");
  // hs->GetXaxis()->SetTitle("Pulse Intergal [% of CE]");
  for(Int_t i = 0; i < hists.size(); i++){
    hs->Add(hists[i][0]);
    hs->Add(hists[i][1]);
    hsSub->Add(hists[i][2]);
  }

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", "Spectrum Comparison", 1200, 800);
  c->Divide(1, 2);
  c->cd(1)->SetGridx();
  c->cd(1)->SetGridy();
  hs->Draw("nostack");
  leg->Draw("same");

  //TCanvas *cSub = new TCanvas("cSub", "Spectrum Subtraction Amount", 1200, 800);
  c->cd(2)->SetGridx();
  c->cd(2)->SetGridy();
  hsSub->Draw("nostack");
  legSub->Draw("same");
}
