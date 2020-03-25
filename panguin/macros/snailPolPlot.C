#include <TCanvas.h>
#include <TPad.h>
#include <TObject.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TPad.h>
#include <TString.h>
#include <TMath.h>
#include <TLegend.h>
#include <vector>

using namespace std;

void snailPolPlot(){
  gStyle->SetOptStat(0);

  float hwpInPol[20] =  {85.177, 87.806, 87.374, 86.687, 85.996, 90.422, 87.659, 88.374, 87.358, 87.138, 
                         87.928, 88.161, 89.045, 88.748, 88.604, 87.430, 86.971, 82.983, 86.270, 88.931};
  float hwpInErr[20] =  {1.086, 1.406, 1.173, 2.456, 1.264, 1.427, 2.000, 2.152, 1.763, 1.238, 
                         1.337, 1.191, 0.567, 0.534, 0.553, 0.875, 0.729, 0.690, 0.702, 0.631};
  int hwpInNum[20] =    {1, 3, 5, 7, 9, 11, 12, 14, 16, 18, 20, 22, 23, 25, 27, 30, 32, 35, 37, 39};
  float hwpOutPol[18] = {88.665, 90.209, 86.462, 86.676, 90.059, 91.141, 87.040, 91.558, 88.305, 89.472, 
                         89.498, 86.885, 87.263, 88.972, 87.151, 86.435, 88.281, 88.186};
  float hwpOutErr[18] = {1.270, 1.093, 1.400, 1.324, 1.383, 1.556, 1.441, 1.658, 1.220, 0.550, 
                         0.539, 0.878, 0.774, 0.895, 0.850, 0.864, 0.724, 1.242};
  int hwpOutNum[18] =   {2, 4, 6, 8, 10, 13, 15, 17, 19, 24, 26, 28, 29, 31, 33, 36, 38, 40};

  TCanvas *cPol = new TCanvas("cPol", "Polarizations by Snail", 1200, 600);
  cPol->Divide(2, 1);
  TH1F *hIN = new TH1F("hIN", "Polarizations by Snail (IHWP IN)", 20, 0, 20);
  TH1F *hOUT = new TH1F("hOUT", "Polarization by Snail (IHWP OUT)", 18, 0, 18);
  TF1 *fconstIN  = new TF1("fconstIN",  "pol0");
  TF1 *fconstOUT = new TF1("fconstOUT", "pol0");
  TPaveText *ptIN =  new TPaveText(0.60, 0.75, 0.95, 0.90, "blNDC");
  TPaveText *ptOUT = new TPaveText(0.60, 0.75, 0.95, 0.90, "blNDC");
  
  for(int i = 0; i < 20; i++){
    hIN->SetBinContent(i + 1, hwpInPol[i]);
    hIN->SetBinError(i + 1, hwpInErr[i]);
    hIN->GetXaxis()->SetBinLabel(i + 1, Form("%i", hwpInNum[i]));
  }
  for(int i = 0; i < 18; i++){
    hOUT->SetBinContent(i + 1, hwpOutPol[i]);
    hOUT->SetBinError(i + 1, hwpOutErr[i]);
    hOUT->GetXaxis()->SetBinLabel(i + 1, Form("%i", hwpOutNum[i]));
  }

  cPol->cd(1);
  hIN->Fit("fconstIN", "Q", "", 0, 20);
  fconstIN->SetLineColor(4);
  ptIN->AddText(Form("--------Polarization (IHWP IN)--------"));
  ptIN->AddText(Form("%.3f%% +/- %.3f",fconstIN->GetParameter(0), fconstIN->GetParError(0)));
  ptIN->AddText(Form("Chi^2 / ndf: %f / %d", fconstIN->GetChisquare(), fconstIN->GetNDF()));
  ptIN->SetBorderSize(1); ptIN->SetFillColor(0);
  hIN->Draw("P");
  ptIN->Draw();

  cPol->cd(2);
  hOUT->Fit("fconstOUT", "Q", "", 0, 18);
  fconstOUT->SetLineColor(2);
  ptOUT->AddText(Form("--------Polarization (IHWP OUT)--------"));
  ptOUT->AddText(Form("%.3f%% +/- %.3f",fconstOUT->GetParameter(0), fconstOUT->GetParError(0)));
  ptOUT->AddText(Form("Chi^2 / ndf: %f / %d", fconstOUT->GetChisquare(), fconstOUT->GetNDF()));
  ptOUT->SetBorderSize(1); ptOUT->SetFillColor(0);
  hOUT->Draw("P");
  ptOUT->Draw();
  
}
