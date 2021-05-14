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

void snailFitRange(Int_t prexOrCrex, Int_t startSnail, Int_t endSnail){
  if(endSnail <= startSnail){
    printf("Invalid snail combination. List snail range from start to finish, including more than one snail.");
  }
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

  TF1 *fit = new TF1("fit", "pol0");

  TGraphErrors *g = new TGraphErrors();
  g->SetMarkerStyle(20);
  // Float_t xmin = 1e16; Float_t xmax = -1e16;
  // Float_t ymin = 1e16; Float_t ymax = -1e16;
  Int_t nPts = 0;

  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(sign == 0 || snlNum == 150 || snlNum == 151 || snlNum == 159 || snlNum == 160 || snlNum == 220 || snlNum == 221) continue;
    g->SetPoint(nPts, snlNum, 100*sign*pol0.mean);
    g->SetPointError(nPts++, 0.0, 100*pol0.meanErr);

    // xmin = (snlNum - 1 < xmin) ? snlNum - 1 : xmin;
    // xmax = (snlNum + 1 > xmax) ? snlNum + 1 : xmax;
    // ymin = (0.95*100*(sign*pol0.mean - pol0.meanErr) < ymin) ? 0.95*100*(sign*pol0.mean - pol0.meanErr) : ymin;
    // ymax = (1.05*100*(sign*pol0.mean + pol0.meanErr) > ymax) ? 1.05*100*(sign*pol0.mean + pol0.meanErr) : ymax;
    // printf("%i,%.4f,%i\n", snlNum, pol0.Chi2, pol0.NDF);
  }

  TCanvas *c = new TCanvas("c", "Pol0 Fits", 1200, 800);
  c->cd()->SetGridx();
  c->cd()->SetGridy();
  // TMultiGraph *mg = new TMultiGraph();
  g->Fit(fit, "Q", "", startSnail, endSnail);
  printf("Range %i-%i: %.4f +/- %.4f\n", startSnail, endSnail, fit->GetParameter(0), fit->GetParError(0));

  g->SetTitle("Polarization Fits (Excluding Snails taken with <100% DOCP)");
  g->GetXaxis()->SetTitle("snailNum");
  g->GetYaxis()->SetTitle("Pol0 [pct]");
  g->Draw("ap");
  
  // TLegend *leg = new TLegend(0.75, 0.35, 0.9, 0.2);
  // leg->AddEntry(gOUT[0], getLegendEntry(0, 0));
  // leg->AddEntry(gIN[0], getLegendEntry(1, 0));
  // leg->AddEntry(gOUT[1], getLegendEntry(0, 1));
  // leg->AddEntry(gIN[1], getLegendEntry(1, 1));
  // leg->Draw("same");

  // for(Int_t i = 0; i < nFitRanges; i++){
  //   texts[i]->SetFillColor(0);
  //   texts[i]->SetFillStyle(0);
  //   texts[i]->SetLineColor(0);
  //   texts[i]->SetLineWidth(0);
  //   texts[i]->Draw("same");
  // }
}
