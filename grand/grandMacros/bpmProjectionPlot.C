#include "makePlots.h"

using namespace std;

void bpmProjectionPlot(Int_t prexOrCrex){
  TFile *f;
  if(prexOrCrex == 1){
    f = TFile::Open(Form("%s/aggregates/prexGrandCompton.root", getenv("COMPMON_WEB")));
  }
  else if(prexOrCrex == 2){
    f = TFile::Open(Form("%s/aggregates/crexGrandCompton.root", getenv("COMPMON_WEB")));
  }
  else{
    printf("ERROR: invalid experiment code\n");
    exit(1);
  }
  TTree *tree = (TTree *)f->Get("run");

  TCanvas *c = new TCanvas("cBPM", "BPM Canvas", 1200, 800);
  TGraphErrors *gLeftIn = new TGraphErrors();
  TGraphErrors *gLeftOut = new TGraphErrors();
  TGraphErrors *gRightIn = new TGraphErrors();
  TGraphErrors *gRightOut = new TGraphErrors();
  TGraphErrors *gCenter = new TGraphErrors();

  DataVar Ax, Ay, Bx, By;
  Int_t runNum, sign;
  Float_t vWien, hWien, solWien, ihwp;
  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress("bpmAx", &Ax);
  tree->SetBranchAddress("bpmAy", &Ay);
  tree->SetBranchAddress("bpmBx", &Bx);
  tree->SetBranchAddress("bpmBy", &By);

  Int_t nLeftIn = 0; Int_t nLeftOut = 0; Int_t nRightIn = 0; Int_t nRightOut = 0;
  Float_t xmin = 1e16; Float_t xmax = -1e16;
  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    Float_t xDiff = Bx.mean - Ax.mean;
    Float_t xDiffErr = TMath::Sqrt(TMath::Power(Bx.meanErr, 2) + TMath::Power(Ax.meanErr, 2));
    Float_t yDiff = By.mean - Ay.mean;
    Float_t yDiffErr = TMath::Sqrt(TMath::Power(By.meanErr, 2) + TMath::Power(Ay.meanErr, 2));
    Float_t projX = 6.0*xDiff; Float_t projXErr = 6.0*xDiffErr;
    Float_t projY = 6.0*yDiff; Float_t projYErr = 6.0*yDiffErr;
    if(projX - projXErr < xmin) xmin = projX - projXErr;
    if(projX + projXErr > xmax) xmax = projX + projXErr;
    if(projY - projYErr < xmin) xmin = projY - projYErr;
    if(projY + projYErr > xmax) xmax = projY + projYErr;

    if(isFlipLeft(runNum, hWien, vWien, solWien, false) && ihwp > 0.5){
      gLeftIn->SetPoint(nLeftIn, projX, projY);
      gLeftIn->SetPointError(nLeftIn++, projXErr, projYErr);
    }
    else if(isFlipLeft(runNum, hWien, vWien, solWien, false) && ihwp < 0.5){
      gLeftOut->SetPoint(nLeftOut, projX, projY);
      gLeftOut->SetPointError(nLeftOut++, projXErr, projYErr);
    }
    else if(!isFlipLeft(runNum, hWien, vWien, solWien, false) && ihwp > 0.5){
      gRightIn->SetPoint(nRightIn, projX, projY);
      gRightIn->SetPointError(nRightIn++, projXErr, projYErr);
    }
    else if(!isFlipLeft(runNum, hWien, vWien, solWien, false) && ihwp < 0.5){
      gRightOut->SetPoint(nRightOut, projX, projY);
      gRightOut->SetPointError(nRightOut++, projXErr, projYErr);
    }
  }

  gCenter->SetPoint(0, 2.160, 0.709); gCenter->SetPointError(0, 0.080, 0.098);
  gLeftIn->SetTitle("Projected Collimator Position");
  gLeftIn->GetXaxis()->SetTitle("x pos (mm)"); gLeftIn->GetYaxis()->SetTitle("y pos (mm)");
  gLeftIn->GetXaxis()->SetLimits(xmin, xmax); gLeftIn->GetYaxis()->SetLimits(xmin, xmax);
  gLeftOut->SetTitle("Projected Collimator Position");
  gLeftOut->GetXaxis()->SetTitle("x pos (mm)"); gLeftOut->GetYaxis()->SetTitle("y pos (mm)");
  gLeftOut->GetXaxis()->SetLimits(xmin, xmax); gLeftOut->GetYaxis()->SetLimits(xmin, xmax);
  gRightIn->SetTitle("Projected Collimator Position");
  gRightIn->GetXaxis()->SetTitle("x pos (mm)"); gRightIn->GetYaxis()->SetTitle("y pos (mm)");
  gRightIn->GetXaxis()->SetLimits(xmin, xmax); gRightIn->GetYaxis()->SetLimits(xmin, xmax);
  gRightOut->SetTitle("Projected Collimator Position");
  gRightOut->GetXaxis()->SetTitle("x pos (mm)"); gRightOut->GetYaxis()->SetTitle("y pos (mm)");
  gRightOut->GetXaxis()->SetLimits(xmin, xmax); gRightOut->GetYaxis()->SetLimits(xmin, xmax);
  gCenter->SetTitle("Projected Collimator Position");
  gCenter->GetXaxis()->SetTitle("x pos (mm)"); gCenter->GetYaxis()->SetTitle("y pos (mm)");
  gCenter->GetXaxis()->SetLimits(xmin, xmax); gCenter->GetYaxis()->SetLimits(xmin, xmax);
  gLeftIn->SetMarkerStyle(getMarkerStyle(1)); gLeftIn->SetMarkerColor(getColor(1));
  gLeftOut->SetMarkerStyle(getMarkerStyle(0)); gLeftOut->SetMarkerColor(getColor(0));
  gRightIn->SetMarkerStyle(getMarkerStyle(3)); gRightIn->SetMarkerColor(getColor(3));
  gRightOut->SetMarkerStyle(getMarkerStyle(2)); gRightOut->SetMarkerColor(getColor(2));
  gCenter->SetMarkerStyle(47); gCenter->SetMarkerColor(1);
  

  gLeftIn->Draw("ap");
  gLeftOut->Draw("p && same");
  gRightIn->Draw("p && same");
  gRightOut->Draw("p && same");
  gCenter->Draw("p && same");
}
