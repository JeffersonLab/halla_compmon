#include "../grandOnline/makePlots.h"

using namespace std;

Int_t getMarkerStyle(Int_t inOrOut, Int_t index){
  if(inOrOut == 0){
    if(index == 1 || index == 2){return 23;}
    else if(index == 0 || index == 3){return 20;}
    else{
      printf("Invalid graph index (marker)!\n");
      return 0;
    }
  }
  else if(inOrOut == 1){
    if(index == 1 || index == 2){return 22;}
    else if(index == 0 || index == 3){return 21;}
    else{
      printf("Invalid graph index (marker)!\n");
      return 0;
    }
  }
  else{
    printf("Invalid IHWP state! (marker)\n");
    return 0;
  }
}

Int_t getColor(Int_t inOrOut, Int_t index){
  if(inOrOut == 0){
    if(index == 1 || index == 2){return kViolet;}
    else if(index == 0 || index == 3){return kRed;}
    else{
      printf("Invalid graph index (color)!\n");
      return 0;
    }
  }
  else if(inOrOut == 1){
    if(index == 1 || index == 2){return kOrange;}
    else if(index == 0 || index == 3){return kBlue;}
    else{
      printf("Invalid graph index (color)!\n");
      return 0;
    }
  }
  else{
    printf("Invalid IHWP state! (color)\n");
    return 0;
  }
}

TString getLegendEntry(Int_t inOrOut, Int_t index){
  TString ihwp = (inOrOut == 0) ? "Out" : "In";
  TString wien = (index == 1 || index == 2) ? "Left" : "Right";
  if(!(inOrOut == 0 || inOrOut == 1) || index < 0 || index > 3){
    printf("Invalid parameters (legend)!\n");
    return "";
  }
  return Form("%s %s", wien.Data(), ihwp.Data());
}

void fitSnailPol(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *snl = (TTree *)f->Get("snl");

  FitPolVar pol0;
  Int_t snlNum, sign;
  Float_t ihwp;

  snl->SetBranchAddress("snailNum", &snlNum);
  snl->SetBranchAddress("Pol0", &pol0);
  snl->SetBranchAddress("sign", &sign);
  snl->SetBranchAddress("ihwp", &ihwp);

  const Int_t nFitRanges = 4;
  Int_t rangeStarts[nFitRanges] = {101, 122, 138, 177};
  Int_t rangeEnds[nFitRanges] = {121, 137, 174, 221};
  // Int_t rangeStarts[nFitRanges] = {101, 138, 177};
  // Int_t rangeEnds[nFitRanges] = {137, 174, 219};
  TF1 *fitsOUT[nFitRanges];
  TF1 *fitsIN[nFitRanges];
  TPaveText *texts[nFitRanges];
  Float_t xLo[nFitRanges] = {0.11, 0.25, 0.45, 0.75};
  Float_t yLo[nFitRanges] = {0.75, 0.20, 0.60, 0.70};
  Float_t xHi[nFitRanges] = {0.21, 0.35, 0.55, 0.85};
  Float_t yHi[nFitRanges] = {0.65, 0.30, 0.70, 0.80};
  // Float_t xLo[nFitRanges] = {0.25, 0.45, 0.75};
  // Float_t yLo[nFitRanges] = {0.20, 0.60, 0.70};
  // Float_t xHi[nFitRanges] = {0.35, 0.55, 0.85};
  // Float_t yHi[nFitRanges] = {0.30, 0.70, 0.80};

  TGraphErrors *gOUT[nFitRanges];
  TGraphErrors *gIN[nFitRanges];
  Float_t xmin = 1e16; Float_t xmax = -1e16;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Int_t nPtsOUT[nFitRanges]; Int_t nPtsIN[nFitRanges];
  for(Int_t i = 0; i < nFitRanges; i++){
    gOUT[i] = new TGraphErrors();
    gOUT[i]->SetMarkerStyle(getMarkerStyle(0, i));
    gOUT[i]->SetMarkerColor(getColor(0, i));
    gIN[i] = new TGraphErrors();
    gIN[i]->SetMarkerStyle(getMarkerStyle(1, i));
    gIN[i]->SetMarkerColor(getColor(1, i));
    TString fitStyle = (i == 1 || i == 2) ? "pol1" : "pol0";

    fitsOUT[i] = new TF1(Form("fit%i-%i_OUT", rangeStarts[i], rangeEnds[i]), fitStyle.Data());
    fitsOUT[i]->SetLineColor(getColor(0, i));
    fitsIN[i] = new TF1(Form("fit%i-%i_IN", rangeStarts[i], rangeEnds[i]), fitStyle.Data());
    fitsIN[i]->SetLineColor(getColor(1, i));

    nPtsOUT[i] = 0; 
    nPtsIN[i] = 0;

    texts[i] = new TPaveText(xLo[i], yLo[i], xHi[i], yHi[i], "blNDC");
  }

  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    //if(sign == 0 || snlNum == 150 || snlNum == 151 || snlNum == 159 || snlNum == 160 || snlNum == 220 || snlNum == 221) continue;
    if(sign == 0) continue;
    Int_t rng = 0;
    for(Int_t j = 0; j < nFitRanges; j++){
      if(snlNum >= rangeStarts[j] && snlNum <= rangeEnds[j])
        rng = j;
    }
    if(ihwp < 0.5){
      gOUT[rng]->SetPoint(nPtsOUT[rng], snlNum, 100*sign*pol0.mean);
      gOUT[rng]->SetPointError(nPtsOUT[rng], 0.0, 100*pol0.meanErr);
      nPtsOUT[rng] = nPtsOUT[rng] + 1;
    }
    else{
      gIN[rng]->SetPoint(nPtsIN[rng], snlNum, 100*sign*pol0.mean);
      gIN[rng]->SetPointError(nPtsIN[rng], 0.0, 100*pol0.meanErr);
      nPtsIN[rng] = nPtsIN[rng] + 1;
    }
    xmin = (snlNum - 1 < xmin) ? snlNum - 1 : xmin;
    xmax = (snlNum + 1 > xmax) ? snlNum + 1 : xmax;
    ymin = (0.95*100*(sign*pol0.mean - pol0.meanErr) < ymin) ? 0.95*100*(sign*pol0.mean - pol0.meanErr) : ymin;
    ymax = (1.05*100*(sign*pol0.mean + pol0.meanErr) > ymax) ? 1.05*100*(sign*pol0.mean + pol0.meanErr) : ymax;
    printf("%i,%.4f,%i\n", snlNum, pol0.Chi2, pol0.NDF);
  }

  TCanvas *c = new TCanvas("c", "Pol0 Fits", 1200, 800);
  c->cd()->SetGridx();
  c->cd()->SetGridy();
  TMultiGraph *mg = new TMultiGraph();
  for(Int_t i = 0; i < nFitRanges; i++){
    gOUT[i]->Fit(fitsOUT[i], "Q", "", rangeStarts[i], rangeEnds[i]);
    gIN[i]->Fit(fitsIN[i], "Q", "", rangeStarts[i], rangeEnds[i]);
    mg->Add(gOUT[i], "p");
    mg->Add(gIN[i], "p");
    texts[i]->AddText(Form("Chi2: %.2f", fitsOUT[i]->GetChisquare()*1.0/fitsOUT[i]->GetNDF()));
    ((TText*)texts[i]->GetListOfLines()->Last())->SetTextColor(getColor(0, i));
    texts[i]->AddText(Form("Chi2: %.2f", fitsIN[i]->GetChisquare()*1.0/fitsIN[i]->GetNDF()));
    ((TText*)texts[i]->GetListOfLines()->Last())->SetTextColor(getColor(1, i));
    printf("Range %i OUT: %.4f +/- %.4f\n", i+1, fitsOUT[i]->GetParameter(0), fitsOUT[i]->GetParError(0));
    printf("Range %i IN: %.4f +/- %.4f\n", i+1, fitsIN[i]->GetParameter(0), fitsIN[i]->GetParError(0));
    if(i == 1 || i == 2){
      printf("Range %i OUT: %.4f +/- %.4f\n", i+1, fitsOUT[i]->GetParameter(1), fitsOUT[i]->GetParError(1));
      printf("Range %i IN: %.4f +/- %.4f\n", i+1, fitsIN[i]->GetParameter(1), fitsIN[i]->GetParError(1));
    }
  }

  mg->SetTitle("Polarization Fits (Excluding Snails taken with <100% DOCP)");
  mg->GetXaxis()->SetTitle("snailNum");
  mg->GetYaxis()->SetTitle("Pol0 [pct]");
  mg->GetXaxis()->SetLimits(xmin, xmax); mg->GetXaxis()->SetRangeUser(xmin, xmax);
  mg->GetYaxis()->SetLimits(84, 90); mg->GetYaxis()->SetRangeUser(84, 90);
  mg->Draw("a");
  
  TLegend *leg = new TLegend(0.75, 0.35, 0.9, 0.2);
  leg->AddEntry(gOUT[0], getLegendEntry(0, 0));
  leg->AddEntry(gIN[0], getLegendEntry(1, 0));
  leg->AddEntry(gOUT[1], getLegendEntry(0, 1));
  leg->AddEntry(gIN[1], getLegendEntry(1, 1));
  leg->Draw("same");

  for(Int_t i = 0; i < nFitRanges; i++){
    texts[i]->SetFillColor(0);
    texts[i]->SetFillStyle(0);
    texts[i]->SetLineColor(0);
    texts[i]->SetLineWidth(0);
    texts[i]->Draw("same");
  }
}
