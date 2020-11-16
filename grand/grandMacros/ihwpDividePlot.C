#include "makePlots.h"

vector<Float_t> legendCoords(Int_t i){
  vector<Float_t> coords;
  if(i == 0){
    coords.push_back(0.10); coords.push_back(0.80);
    coords.push_back(0.25); coords.push_back(0.90);
  }
  else if(i == 1){
    coords.push_back(0.25); coords.push_back(0.80);
    coords.push_back(0.40); coords.push_back(0.90);
  }
  else if(i == 2){
    coords.push_back(0.75); coords.push_back(0.80);
    coords.push_back(0.90); coords.push_back(0.90);
  }
  else if(i == 3){
    coords.push_back(0.60); coords.push_back(0.80);
    coords.push_back(0.75); coords.push_back(0.90);
  }
  else{
    printf("Invalid graph index (legend)!\n");
  }
  return coords;
}

void ihwpDividePlot(Int_t snlStart, Int_t snlEnd){
  TString expt("");
  if(snlStart <= 99 && snlEnd >= 100){
    printf("Please no mixed-experiment snails.\n");
    exit(1);
  }
  else if(snlStart <= 99){
    expt = "prex";
  }
  else{
    expt = "crex";
  }

  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *snl = (TTree *)f->Get("snl");
  vector<TGraphErrors *> graphs;
  vector<Int_t> graphCounts;
  vector<TF1 *> fits;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;
  Float_t xmins[4] = {1e10, 1e10, 1e10, 1e10};
  Float_t xmaxs[4] = {-1e10, -1e10, -1e10, -1e10};

  FitPolVar pol0;
  Int_t snailNum, sign, nCycles;
  Float_t hWien, vWien, solWien, ihwp;
  snl->SetBranchAddress("snailNum", &snailNum);
  snl->SetBranchAddress("numCyclesAcc", &nCycles);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("VWienAngle", &vWien);
  snl->SetBranchAddress("HWienAngle", &hWien);
  snl->SetBranchAddress("PhiFG", &solWien);
  snl->SetBranchAddress("sign", &sign);
  snl->SetBranchAddress("Pol0", &pol0);

  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("Pol0 vs Snail (Snails %i-%i By IHWP)", snlStart, snlEnd));
    g->GetXaxis()->SetTitle("snailNum"); g->GetYaxis()->SetTitle("Pol0 (pct)");
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
    TF1 *fit = new TF1(Form("fit%i", i), "pol0");
    fit->SetLineColor(getColor(i));
    fits.push_back(fit);
  }

  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    if(sign == 0 || snailNum < snlStart || snailNum > snlEnd || nCycles == 0) continue;
    Int_t ind = getGraphInd(snailNum, hWien, vWien, solWien, ihwp, true);
    Float_t polVar = 100.*pol0.mean*sign; Float_t polErr = 100.*pol0.meanErr;
    graphs[ind]->SetPoint(graphCounts[ind], snailNum, polVar);
    graphs[ind]->SetPointError(graphCounts[ind], 0.0, TMath::Abs(polErr));
    graphCounts[ind] += 1;
    if(polVar - polErr < ymin){ymin = polVar - polErr;}
    if(polVar + polErr > ymax){ymax = polVar + polErr;}
    if(snailNum - 1 < xmin){xmin = snailNum - 1;}
    if(snailNum + 1 > xmax){xmax = snailNum + 1;}
    if(snailNum < xmins[ind]){xmins[ind] = snailNum;}
    if(snailNum > xmaxs[ind]){xmaxs[ind] = snailNum;}
  }

  adjustGraphLimits(graphs, ymin, ymax, 0.99, 1.01, xmin, xmax);

  TCanvas *c = new TCanvas("c", "Pol0 Graph Comparison Canvas", 1200, 800);
  c->cd();

  Bool_t drawn = false;
  for(Int_t i = 0; i < 4; i++){
    if(graphs[i]->GetN() > 0 && !drawn){
      graphs[i]->Draw("ap");
      drawn = true;
    }
    else if(graphs[i]->GetN() > 0 && drawn){
      graphs[i]->Draw("p && same");
    }
  }
  for(Int_t i = 0; i < 4; i++){
    graphs[i]->Fit(Form("fit%i", i), "Q", "", xmins[i], xmaxs[i]);
    vector<Float_t> coords = legendCoords(i);
    TPaveText *pt = new TPaveText(coords[0], coords[1], coords[2], coords[3], "blNDC");
    TText *line1 = pt->AddText(Form("------%s------\n", getLegendEntry(i).Data()));
    line1->SetTextColor(getColor(i));
    TText *line2 = pt->AddText(Form("Pol0: %.3f +/- %.3f\n", fits[i]->GetParameter(0), fits[i]->GetParError(0)));
    line2->SetTextColor(getColor(i));
    TText *line3 = pt->AddText(Form("Chi2 / NDF: %.3f / %i\n", fits[i]->GetChisquare(), fits[i]->GetNDF()));
    line3->SetTextColor(getColor(i));
    pt->SetFillColor(0); pt->SetBorderSize(1);
    pt->Draw("same");
  }
}
