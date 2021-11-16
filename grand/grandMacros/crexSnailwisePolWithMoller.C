#include "../grandOnline/makePlots.h"

using namespace std;

const Int_t nWiens = 2;
const Int_t nIHWPs = 2;

TGraphErrors *gCompton[nWiens][nIHWPs];
TGraphErrors *gMoller[nWiens][nIHWPs];
Int_t nPtsCompton[nWiens][nIHWPs]; 
Int_t nPtsMoller[nWiens][nIHWPs];

Int_t markerStyles[nWiens][nIHWPs] = {{23, 22}, {20, 21}};
Int_t colors[nWiens][nIHWPs] = {{kViolet, kOrange}, {kRed, kBlue}};
TString legEntry[nWiens][nIHWPs] = {{"Left Out", "Left In"}, {"Right Out", "Right In"}};

Int_t mollerMarkerStyles[nWiens][nIHWPs] = {{47, 29}, {33, 34}};
Int_t mollerColors[nWiens][nIHWPs] = {{kTeal, kOrange+2}, {kGreen+1, kMagenta-4}};


void fillMollerGraphs(){
  gMoller[1][0]->SetPoint(0, 80.0, 86.85);
  gMoller[1][0]->SetPointError(0, 0.0, 0.19);

  gMoller[1][1]->SetPoint(0, 78.0, 87.01);
  gMoller[1][1]->SetPointError(0, 0.0, 0.19);

  gMoller[1][1]->SetPoint(1, 106.5, 86.39);
  gMoller[1][1]->SetPointError(1, 0.0, 0.17);

  gMoller[1][0]->SetPoint(1, 106.5, 87.12);
  gMoller[1][0]->SetPointError(1, 0.0, 0.51);

  gMoller[1][0]->SetPoint(2, 116.5, 86.98);
  gMoller[1][0]->SetPointError(2, 0.0, 0.20);
  
  gMoller[1][1]->SetPoint(2, 116.5, 86.83);
  gMoller[1][1]->SetPointError(2, 0.0, 0.20);

  gMoller[0][1]->SetPoint(0, 127, 86.23);
  gMoller[0][1]->SetPointError(0, 0.0, 0.16);

  gMoller[0][0]->SetPoint(0, 127, 86.79);
  gMoller[0][0]->SetPointError(0, 0.0, 0.16);

  gMoller[0][0]->SetPoint(1, 161.5, 87.08);
  gMoller[0][0]->SetPointError(1, 0.0, 0.23);

  gMoller[1][1]->SetPoint(3, 182.5, 87.43);
  gMoller[1][1]->SetPointError(3, 0.0, 0.20);

  gMoller[1][0]->SetPoint(3, 182.5, 87.68);
  gMoller[1][0]->SetPointError(3, 0.0, 0.20);

  gMoller[1][0]->SetPoint(4, 207.5, 87.45);
  gMoller[1][0]->SetPointError(4, 0.0, 0.16);

  gMoller[1][1]->SetPoint(4, 207.5, 87.49);
  gMoller[1][1]->SetPointError(4, 0.0, 0.17);

  gMoller[1][0]->SetPoint(5, 222.0, 87.62);
  gMoller[1][0]->SetPointError(5, 0.0, 0.18);

  gMoller[1][1]->SetPoint(5, 222.0, 87.11);
  gMoller[1][1]->SetPointError(5, 0.0, 0.18);
}


void crexSnailwisePolWithMoller(){
  TString expt = experimentCode(2);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *snl = (TTree *)f->Get("snl");

  FitPolVar pol0;
  Int_t snlNum, sign;
  Float_t ihwp, vWien;

  snl->SetBranchAddress("snailNum", &snlNum);
  snl->SetBranchAddress("Pol0", &pol0);
  snl->SetBranchAddress("sign", &sign);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("VWienAngle", &vWien);

  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      gCompton[i][j] = new TGraphErrors();
      gCompton[i][j]->SetMarkerStyle(markerStyles[i][j]);
      gCompton[i][j]->SetMarkerSize(1.5);
      gCompton[i][j]->SetMarkerColor(colors[i][j]);
      gMoller[i][j] = new TGraphErrors();
      gMoller[i][j]->SetMarkerStyle(mollerMarkerStyles[i][j]);
      gMoller[i][j]->SetMarkerSize(1.5);
      gMoller[i][j]->SetMarkerColor(mollerColors[i][j]);

      nPtsCompton[i][j] = 0; 
      nPtsMoller[i][j] = 0;
    }
  }

  for(Int_t i = 0; i < snl->GetEntries(); i++){
    snl->GetEntry(i);
    //if(sign == 0 || snlNum == 150 || snlNum == 151 || snlNum == 159 || snlNum == 160 || snlNum == 220 || snlNum == 221) continue;
    if(sign == 0) continue;
    Int_t wien = (vWien < 0) ? 0 : 1;
    Int_t ihwpInd = (ihwp < 0.5) ? 0 : 1;
    gCompton[wien][ihwpInd]->SetPoint(nPtsCompton[wien][ihwpInd], snlNum, 100*sign*pol0.mean);
    gCompton[wien][ihwpInd]->SetPointError(nPtsCompton[wien][ihwpInd], 0.0, 100*pol0.meanErr);
    nPtsCompton[wien][ihwpInd]++;

    printf("%i,%.4f,%i\n", snlNum, pol0.Chi2, pol0.NDF);
  }

  fillMollerGraphs();

  TCanvas *c = new TCanvas("c", "Pol0 Fits", 1500, 500);
  c->cd()->SetGridx();
  c->cd()->SetGridy();
  TLegend *leg = new TLegend(0.75, 0.35, 0.9, 0.1);
  TMultiGraph *mg = new TMultiGraph();
  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      mg->Add(gCompton[i][j], "p");
      leg->AddEntry(gCompton[i][j], Form("%s (Compton)", legEntry[i][j].Data()));
    }
  }

  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      mg->Add(gMoller[i][j], "p");
      leg->AddEntry(gMoller[i][j], Form("%s (Moller)", legEntry[i][j].Data()));
    }
  }

  mg->SetTitle("CREX Polarization Measurements (Compton & Moller)");
  mg->GetXaxis()->SetTitle("Approx Time of Measurement");
  mg->GetYaxis()->SetTitle("Beam Polarization [pct]");
  // mg->GetXaxis()->SetLimits(xmin, xmax); mg->GetXaxis()->SetRangeUser(xmin, xmax);
  // mg->GetYaxis()->SetLimits(84, 90); mg->GetYaxis()->SetRangeUser(84, 90);
  mg->Draw("a");
  leg->Draw("same");
}
