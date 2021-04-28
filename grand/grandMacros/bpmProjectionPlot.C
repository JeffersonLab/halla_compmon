#include "../grandOnline/makePlots.h"

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
  TTree *tree = (TTree *)f->Get("cyc");

  TCanvas *c = new TCanvas("cBPM", "BPM Canvas", 1200, 800);
  TGraphErrors *gLeftIn = new TGraphErrors();
  TGraphErrors *gLeftOut = new TGraphErrors();
  TGraphErrors *gRightIn = new TGraphErrors();
  TGraphErrors *gRightOut = new TGraphErrors();
  TGraphErrors *gCenter = new TGraphErrors();

  DataVar Ax, Ay, Bx, By;
  Int_t runNum, sign, cycCut;
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
  tree->SetBranchAddress("CycleCut", &cycCut);

  Int_t nLeftIn = 0; Int_t nLeftOut = 0; Int_t nRightIn = 0; Int_t nRightOut = 0;
  Float_t xmin = 1e16; Float_t xmax = -1e16;
  Float_t maxOffset = 0.0;
  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(cycCut != 0){continue;}
    Float_t projX = Ax.mean + 6.0*(Bx.mean - Ax.mean);
    Float_t projY = Ay.mean + 6.0*(By.mean - Ay.mean);
    Float_t xDiffErr = 6.0*TMath::Sqrt(TMath::Power(Bx.meanErr, 2) + TMath::Power(Ax.meanErr, 2));
    Float_t yDiffErr = 6.0*TMath::Sqrt(TMath::Power(By.meanErr, 2) + TMath::Power(Ay.meanErr, 2));
    Float_t projXErr = TMath::Sqrt(TMath::Power(Ax.meanErr, 2) + TMath::Power(xDiffErr, 2));
    Float_t projYErr = TMath::Sqrt(TMath::Power(Ay.meanErr, 2) + TMath::Power(yDiffErr, 2));
    if(1.05*(projX - projXErr) < xmin) xmin = 1.05*(projX - projXErr);
    if(1.05*(projX + projXErr) > xmax) xmax = 1.05*(projX + projXErr);
    if(1.05*(projY - projYErr) < xmin) xmin = 1.05*(projY - projYErr);
    if(1.05*(projY + projYErr) > xmax) xmax = 1.05*(projY + projYErr);
    maxOffset = (TMath::Sqrt(TMath::Power(projX, 2) + TMath::Power(projY, 2)) > maxOffset) ? TMath::Sqrt(TMath::Power(projX, 2) + TMath::Power(projY, 2)) : maxOffset;

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
  printf("Graph Counts: %i, %i, %i, %i\n", nRightOut, nLeftOut, nRightIn, nLeftIn);
  printf("Graph Limits: %.4f, %.4f\n", xmin, xmax);
  printf("Max Offset: %.4f\n", maxOffset);

  if(prexOrCrex==2){
    gCenter->SetPoint(0, -0.133, 0.747); gCenter->SetPointError(0, 0.038, 0.021);
  }
  else if(prexOrCrex==1){
    gCenter->SetPoint(0, 2.160, 0.709); gCenter->SetPointError(0, 0.080, 0.098);
  }
  gLeftIn->SetTitle("Projected Collimator Position");
  gLeftIn->GetXaxis()->SetTitle("x pos (mm)"); gLeftIn->GetYaxis()->SetTitle("y pos (mm)");
  gLeftIn->GetXaxis()->SetLimits(xmin, xmax); gLeftIn->GetYaxis()->SetLimits(xmin, xmax);
  gLeftIn->GetYaxis()->SetRangeUser(xmin, xmax);
  gLeftOut->SetTitle("Projected Collimator Position");
  gLeftOut->GetXaxis()->SetTitle("x pos (mm)"); gLeftOut->GetYaxis()->SetTitle("y pos (mm)");
  gLeftOut->GetXaxis()->SetLimits(xmin, xmax); gLeftOut->GetYaxis()->SetLimits(xmin, xmax);
  gLeftOut->GetYaxis()->SetRangeUser(xmin, xmax);
  gRightIn->SetTitle("Projected Collimator Position");
  gRightIn->GetXaxis()->SetTitle("x pos (mm)"); gRightIn->GetYaxis()->SetTitle("y pos (mm)");
  gRightIn->GetXaxis()->SetLimits(xmin, xmax); gRightIn->GetYaxis()->SetLimits(xmin, xmax);
  gRightIn->GetYaxis()->SetRangeUser(xmin, xmax);
  gRightOut->SetTitle("Projected Collimator Position");
  gRightOut->GetXaxis()->SetTitle("x pos (mm)"); gRightOut->GetYaxis()->SetTitle("y pos (mm)");
  gRightOut->GetXaxis()->SetLimits(xmin, xmax); gRightOut->GetYaxis()->SetLimits(xmin, xmax);
  gRightOut->GetYaxis()->SetRangeUser(xmin, xmax);
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
