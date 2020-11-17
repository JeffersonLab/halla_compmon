#ifndef makePlots_h
#define makePlots_h

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "../grandConstruction/vars.h"

TString experimentCode(Int_t prexOrCrex){
  TString expt("prex");
  if(prexOrCrex == 2){
    expt = "crex";
  }
  else if(prexOrCrex != 1 && prexOrCrex != 2){
    printf("Please enter correct parameter for experiment:\n");
    printf("    PREX == 1\n");
    printf("    CREX == 2\n\n");
    exit(1);
  }
  return expt;
}

Int_t getMarkerStyle(Int_t index){
  if(index == 0){return 23;}
  else if(index == 1){return 22;}
  else if(index == 2){return 20;}
  else if(index == 3){return 21;}
  else{
    printf("Invalid graph index (marker)!\n");
    return 0;
  }
}

Int_t getColor(Int_t index){
  if(index == 0){return kViolet;}
  else if(index == 1){return kOrange;}
  else if(index == 2){return kRed;}
  else if(index == 3){return kBlue;}
  else{
    printf("Invalid graph index (color)!\n");
    return 0;
  }
}

TString getLegendEntry(Int_t index){
  if(index == 0){return "Left Out";}
  else if(index == 1){return "Left In";}
  else if(index == 2){return "Right Out";}
  else if(index == 3){return "Right In";}
  else{
    printf("Invalid graph index (legend)!\n");
    return "";
  }
}

void adjustMinMaxs(vector<Float_t> ymins, vector<Float_t> ymaxs, Float_t val, Int_t msmt){
  if(val < ymins[msmt]){ymins[msmt] = val;}
  if(val > ymaxs[msmt]){ymaxs[msmt] = val;}
}

Bool_t isCloseTo(Float_t num, Float_t ref){
  Bool_t val = TMath::Abs(num)>=0.99*TMath::Abs(ref) && TMath::Abs(num)<=1.01*TMath::Abs(ref);
  return val;
}

Bool_t isFlipLeft(Int_t num, Float_t HWienAngle, Float_t VWienAngle, Float_t PhiFG, Bool_t isSnail=true){
  if((num < 99 && isSnail) || (num < 4900 && !isSnail)){
    if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, -89.9008) && isCloseTo(PhiFG, 86.9014))
      return true;
    else if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 91.1902))
      return false;
  }
  else{
    if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, -90.5996) && isCloseTo(PhiFG, 88.0277))
      return true;
    else if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 89.9558))
      return false;
  }
  return false;
}

Int_t getGraphInd(Int_t snailNum, Float_t hWien, Float_t vWien, Float_t solWien, Float_t ihwp, Bool_t isSnail=true){
  if(isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp <= 0.5){return 0;}
  else if(isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp > 0.5){return 1;}
  else if(!isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp <= 0.5){return 2;}
  else if(!isFlipLeft(snailNum, hWien, vWien, solWien, isSnail) && ihwp > 0.5){return 3;}
  else{printf("Invalid graph index (graph)!\n"); return 0;}
}

void adjustGraphLimits(vector<TGraphErrors *> graphs, Float_t ymin, Float_t ymax, Float_t smallFac, Float_t largeFac, Float_t xmin, Float_t xmax){
  Float_t yminFac = smallFac; Float_t ymaxFac = largeFac;
  if(ymin < 0){yminFac = largeFac;}
  if(ymax < 0){ymaxFac = smallFac;}
  for(Int_t i = 0; i < 4; i++){
    graphs[i]->GetXaxis()->SetLimits(xmin, xmax);
    graphs[i]->GetYaxis()->SetRangeUser(yminFac*ymin, ymaxFac*ymax);
  }
}

void drawAndPrintGraphs(TString msmt, Int_t msmtNum, TCanvas *can, vector<TGraphErrors *> graphs){
  can->cd();
  can->SetGridx(); can->SetGridy();
  //TGraphErrors *g1 = new TGraphErrors(); TGraphErrors *g2 = new TGraphErrors();
  //g1->SetPoint(0,  0, 90.31); g1->SetPointError(0, 0.0, 0.72248);
  //g1->SetPoint(1, 44, 90.31); g1->SetPointError(1, 0.0, 0.72248);
  //g2->SetPoint(0,  0, 89.10); g2->SetPointError(0, 0.0, 0.71280);
  //g2->SetPoint(1, 44, 89.10); g2->SetPointError(1, 0.0, 0.71280);
  //g1->SetFillColorAlpha(kOrange + 7, 0.5); g2->SetFillColorAlpha(kGreen + 2, 0.5);
  //g1->SetLineColor(kOrange + 7); g2->SetLineColor(kGreen + 2);
  //g1->GetXaxis()->SetLimits(0, 45); g1->GetYaxis()->SetRangeUser(80, 100);
  //g2->GetXaxis()->SetLimits(0, 45); g2->GetYaxis()->SetRangeUser(80, 100);
  //g1->Draw("a3l"); 
  //g2->Draw("a3l");
  TLegend *leg = new TLegend(0.9, 0.75, 0.98, 0.9, "", "NDC");
  //TLegend *leg = new TLegend(0.9, 0.7, 0.98, 0.9);
  for(Int_t i = 0; i < 4; i++){
    //if(i == 1 || i == 3) continue;
    //if(i == 0 || i == 2) continue;
    leg->AddEntry(graphs[i], getLegendEntry(i).Data());
    graphs[i]->SetMarkerSize(1.5); graphs[i]->SetLineWidth(1);
    if(i == 0){graphs[i]->Draw("ap");}
    else{graphs[i]->Draw("p && same");}
    //graphs[i]->Draw("p && same");
  }
  //leg->AddEntry(g1, "Moller IHWP OUT");
  //leg->AddEntry(g2, "Moller IHWP IN");
  //leg->Draw("same");

  can->Print(Form("%s/grandOnline/plots/msmt%04i_%s.pdf", getenv("COMPMON_GRAND"), msmtNum, msmt.Data()), "pdf");
  //exit(0);
}

//||====================================||
//||                                    ||
//||  //====\\  |\   || ||              ||
//|| //      \\ |\\  || ||              ||
//|| \\         ||\\ || ||              ||
//||  \\====\\  || \\|| ||              ||
//||         \\ ||  \|| ||              ||
//||         || ||   \| ||              ||
//||         || ||   || ||              ||
//|| \\      // ||   || ||              ||
//||  \\====//  ||   || ||======        ||
//||                                    ||
//||====================================||

void plotPolSnl(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool signCorr=true, bool chi2=false, Float_t factor=1.0){
  FitPolVar var;
  vector<TGraphErrors *> graphs;
  vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  TFile *grand = new TFile(Form("%s/aggregates/%s", getenv("COMPMON_WEB"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("snl");
  TString chiID(""); TString signID("");
  if(chi2) chiID = "chi2";
  if(signCorr) signID = "sign";
  TCanvas *can = new TCanvas(Form("can%s%s_%s", signID.Data(), chiID.Data(), msmt.Data()), "Some Title", 1200, 400);
  TString titleAdd("");
  if(signCorr) titleAdd = " (Sign Corrected)";
  TString chiAdd("");
  if(chi2) chiAdd = " (Chi2 / NDF) ";

  Int_t snailNum, sign;
  Float_t vWien, hWien, solWien, ihwp;
  tree->SetBranchAddress("snailNum", &snailNum);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress(msmt.Data(), &var);
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s%s%s vs Snail", msmt.Data(), chiAdd.Data(), titleAdd.Data()));
    g->GetXaxis()->SetTitle("snailNum"); g->GetYaxis()->SetTitle(msmt.Data());
    if(chi2){g->GetYaxis()->SetTitle(Form("%s (Chi2 / NDF)", msmt.Data()));}
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(snailNum, hWien, vWien, solWien, ihwp, true);
    Float_t polVar = factor*var.mean; Float_t polErr = factor*var.meanErr;
    if(signCorr){polVar = factor*var.mean*sign;}
    if(chi2){polVar = var.Chi2*1.0/var.NDF; polErr = 0.0;}
    graphs[ind]->SetPoint(graphCounts[ind], snailNum, polVar);
    graphs[ind]->SetPointError(graphCounts[ind], 0.0, TMath::Abs(polErr));
    graphCounts[ind] += 1;
    if(polVar - polErr < ymin){ymin = polVar - polErr;}
    if(polVar + polErr > ymax){ymax = polVar + polErr;}
    if(snailNum - 1 < xmin){xmin = snailNum - 1;}
    if(snailNum + 1 > xmax){xmax = snailNum + 1;}
    //printf("<tr class=\"myRow\"><td class=\"myCell\">%i</td><td class=\"myCell\">%.4f +/ %.4f</td></tr>\n", snailNum, 1000*polVar, 1000*polErr);
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

void plotStandardSnl(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool isFloat=true){
  Int_t snailNum, sign, intVar;
  Float_t hWien, vWien, solWien, ihwp, fltVar;

  TFile *grand = new TFile(Form("%s/aggregates/%s", getenv("COMPMON_WEB"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("snl");
  TCanvas *can = new TCanvas(Form("can_%s", msmt.Data()), "", 1200, 400); 
  vector<TGraphErrors *> graphs; vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  tree->SetBranchAddress("snailNum", &snailNum);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("sign", &sign);
  if(isFloat)
    tree->SetBranchAddress(msmt.Data(), &fltVar);
  else
    tree->SetBranchAddress(msmt.Data(), &intVar);
  
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s vs Snail", msmt.Data()));
    g->GetXaxis()->SetTitle("snailNum"); g->GetYaxis()->SetTitle(msmt.Data());
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(snailNum, hWien, vWien, solWien, ihwp, true);
    if(isFloat){
      graphs[ind]->SetPoint(graphCounts[ind], snailNum, fltVar);
      if(fltVar < ymin){ymin = fltVar;}
      if(fltVar > ymax){ymax = fltVar;}
    }
    else{
      graphs[ind]->SetPoint(graphCounts[ind], snailNum, intVar);
      if(((Float_t)intVar) < ymin){ymin = (Float_t)intVar;}
      if(((Float_t)intVar) > ymax){ymax = (Float_t)intVar;}
    }
    if(snailNum - 1 < xmin){xmin = snailNum - 1;}
    if(snailNum + 1 > xmax){xmax = snailNum + 1;}
    graphCounts[ind] += 1;
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

//||====================================||
//||                                    ||
//|| ||====\\ ||     || |\   ||         ||
//|| ||    || ||     || |\\  ||         ||
//|| ||    || ||     || ||\\ ||         ||
//|| ||====// ||     || || \\||         ||
//|| ||\\     ||     || ||  \\|         ||
//|| || \\    ||     || ||   \|         ||
//|| ||  \\   ||     || ||   ||         ||
//|| ||   \\  \\     // ||   ||         ||
//|| ||    \\  \\===//  ||   ||         ||
//||                                    ||
//||====================================||

void plotPolRun(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool signCorr=true, bool chi2=false, Float_t factor=1.0){
  FitPolVar var;
  vector<TGraphErrors *> graphs;
  vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  TFile *grand = new TFile(Form("%s/aggregates/%s", getenv("COMPMON_WEB"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("run");
  TString chiID(""); TString signID("");
  if(chi2) chiID = "chi2";
  if(signCorr) signID = "sign";
  TString canName = Form("can%s%s_%s", signID.Data(), chiID.Data(), msmt.Data());
  TCanvas *can = new TCanvas(canName.Data(), "Some Title", 1200, 400);
  TString titleAdd(""); TString chiAdd("");
  if(signCorr) titleAdd = " (Sign Corrected)";
  if(chi2) chiAdd = " (Chi2 / NDF) ";

  Int_t runNum, sign;
  Float_t vWien, hWien, solWien, ihwp;
  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress(msmt.Data(), &var);
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s%s%s vs Run", msmt.Data(), chiAdd.Data(), titleAdd.Data()));
    g->GetXaxis()->SetTitle("runNum"); g->GetYaxis()->SetTitle(msmt.Data());
    if(chi2) g->GetYaxis()->SetTitle(Form("%s (Chi2 / NDF)", msmt.Data()));
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
    Int_t ind = getGraphInd(runNum, hWien, vWien, solWien, ihwp, false);
    Float_t polVar = factor*var.mean; Float_t polErr = factor*var.meanErr;
    if(sign == 0 || TMath::IsNaN(polVar) || var.NDF==0) continue;
    if(signCorr){polVar = factor*var.mean*sign;}
    if(chi2){polVar = var.Chi2*1.0/var.NDF; polErr = 0.0;}
    graphs[ind]->SetPoint(graphCounts[ind], runNum, polVar);
    graphs[ind]->SetPointError(graphCounts[ind], 0.0, TMath::Abs(polErr));
    graphCounts[ind] += 1;
    if(polVar - polErr < ymin){ymin = polVar - polErr;}
    if(polVar + polErr > ymax){ymax = polVar + polErr;}
    if(runNum - 5 < xmin){xmin = runNum - 5;}
    if(runNum + 5 > xmax){xmax = runNum + 5;}
    //printf("<tr class=\"myRow\"><td class=\"myCell\">%i</td><td class=\"myCell\">%.4f +/ %.4f</td></tr>\n", snailNum, 1000*polVar, 1000*polErr);
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

void plotStandardRun(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool isFloat=true){
  Int_t runNum, sign, intVar;
  Float_t hWien, vWien, solWien, ihwp, fltVar;

  TFile *grand = new TFile(Form("%s/aggregates/%s", getenv("COMPMON_WEB"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("run");
  TCanvas *can = new TCanvas(Form("can_%s", msmt.Data()), "", 1200, 400); 
  vector<TGraphErrors *> graphs; vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 1e10; Float_t xmax = -1e10;

  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("sign", &sign);
  if(isFloat)
    tree->SetBranchAddress(msmt.Data(), &fltVar);
  else
    tree->SetBranchAddress(msmt.Data(), &intVar);
  
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s vs Run", msmt.Data()));
    g->GetXaxis()->SetTitle("runNum"); g->GetYaxis()->SetTitle(msmt.Data());
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0) continue;
    Int_t ind = getGraphInd(runNum, hWien, vWien, solWien, ihwp, false);
    if(isFloat){
      graphs[ind]->SetPoint(graphCounts[ind], runNum, fltVar);
      if(fltVar < ymin){ymin = fltVar;}
      if(fltVar > ymax){ymax = fltVar;}
    }
    else{
      graphs[ind]->SetPoint(graphCounts[ind], runNum, intVar);
      if(((Float_t)intVar) < ymin){ymin = (Float_t)intVar;}
      if(((Float_t)intVar) > ymax){ymax = (Float_t)intVar;}
    }
    if(runNum - 5 < xmin){xmin = runNum - 5;}
    if(runNum + 5 > xmax){xmax = runNum + 5;}
    graphCounts[ind] += 1;
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

//||====================================||
//||                                    ||
//||  //====\\  \\    //  //====\\      ||
//|| //      \\  \\  //  //      \\     ||
//|| ||           \\//   ||             ||
//|| ||            \/    ||             ||
//|| ||            ||    ||             ||
//|| ||            ||    ||             ||
//|| ||            ||    ||             ||
//|| \\      //    ||    \\      //     ||
//||  \\====//     ||     \\====//      ||
//||                                    ||
//||====================================||

void plotPolCyc(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool signCorr=true, Float_t factor=1.0){
  PolVar var;
  vector<TGraphErrors *> graphs;
  vector<Int_t> graphCounts;

  TFile *grand = new TFile(Form("%s/aggregates/%s", getenv("COMPMON_WEB"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("cyc");
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 0; Float_t xmax = 1.05*tree->GetEntries("sign!=0 && CycleCut==0");
  TString signID("");
  if(signCorr) signID = "sign";
  TString canName = Form("can%s_%s", signID.Data(), msmt.Data());
  TCanvas *can = new TCanvas(canName.Data(), "Some Title", 1200, 400);
  TString titleAdd(""); TString chiAdd("");
  if(signCorr) titleAdd = " (Sign Corrected)";

  Int_t runNum, sign, cycCut;
  Float_t vWien, hWien, solWien, ihwp;
  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress("CycleCut", &cycCut);
  tree->SetBranchAddress(msmt.Data(), &var);
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s%s vs Cycle", msmt.Data(), titleAdd.Data()));
    g->GetXaxis()->SetTitle("cycNum"); g->GetYaxis()->SetTitle(msmt.Data());
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign == 0 || cycCut != 0) continue;
    Int_t ind = getGraphInd(runNum, hWien, vWien, solWien, ihwp, false);
    Float_t polVar = factor*var.mean; Float_t polErr = factor*var.meanErr;
    if(signCorr){polVar = factor*var.mean*sign;}
    if(TMath::IsNaN(polVar)) continue;
    graphs[ind]->SetPoint(graphCounts[ind], graphCounts[0] + graphCounts[1] + graphCounts[2] + graphCounts[3], polVar);
    graphs[ind]->SetPointError(graphCounts[ind], 0.0, polErr);
    graphCounts[ind] += 1;
    if(polVar - polErr < ymin){ymin = polVar - polErr;}
    if(polVar + polErr > ymax){ymax = polVar + polErr;}
    //if(runNum - 5 < xmin){xmin = runNum - 5;}
    //if(runNum + 5 > xmax){xmax = runNum + 5;}
    //printf("<tr class=\"myRow\"><td class=\"myCell\">%i</td><td class=\"myCell\">%.4f +/ %.4f</td></tr>\n", snailNum, 1000*polVar, 1000*polErr);
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

void plotStandardCyc(TString fname, TString msmt, Int_t msmtNum, Float_t smallFac, Float_t largeFac, bool signCorr=false){
  Int_t runNum, sign, cycCut;
  Float_t hWien, vWien, solWien, ihwp;
  DataVar data;

  TFile *grand = new TFile(Form("%s/aggregates/%s", getenv("COMPMON_WEB"), fname.Data()), "READ");
  TTree *tree = (TTree *)grand->Get("cyc");
  TCanvas *can = new TCanvas(Form("can_%s", msmt.Data()), "", 1200, 400); 
  vector<TGraphErrors *> graphs; vector<Int_t> graphCounts;
  Float_t ymin = 1e16; Float_t ymax = -1e16;
  Float_t xmin = 0; Float_t xmax = 1.05*tree->GetEntries("sign!=0 && CycleCut==0");

  tree->SetBranchAddress("runNum", &runNum);
  tree->SetBranchAddress("HWienAngle", &hWien);
  tree->SetBranchAddress("VWienAngle", &vWien);
  tree->SetBranchAddress("PhiFG", &solWien);
  tree->SetBranchAddress("ihwp", &ihwp);
  tree->SetBranchAddress("sign", &sign);
  tree->SetBranchAddress("CycleCut", &cycCut);
  tree->SetBranchAddress(msmt.Data(), &data);
  
  for(Int_t i = 0; i < 4; i++){
    TGraphErrors *g = new TGraphErrors();
    g->SetTitle(Form("%s vs Cycle", msmt.Data()));
    g->GetXaxis()->SetTitle("cycNum"); g->GetYaxis()->SetTitle(msmt.Data());
    g->SetMarkerStyle(getMarkerStyle(i)); g->SetMarkerColor(getColor(i));
    graphs.push_back(g); graphCounts.push_back(0);
  }

  for(Int_t i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if(sign==0 || cycCut==1) continue;
    Int_t ind = getGraphInd(runNum, hWien, vWien, solWien, ihwp, false);
    if(signCorr){
      graphs[ind]->SetPoint(graphCounts[ind], graphCounts[0] + graphCounts[1] + graphCounts[2] + graphCounts[3], data.mean*sign);
      graphs[ind]->SetPointError(graphCounts[ind], 0.0, data.meanErr);
      if(data.mean*sign - data.meanErr < ymin){ymin = data.mean*sign - data.meanErr;}
      if(data.mean*sign + data.meanErr > ymax){ymax  =data.mean*sign + data.meanErr;}
    }
    else{
      graphs[ind]->SetPoint(graphCounts[ind], graphCounts[0] + graphCounts[1] + graphCounts[2] + graphCounts[3], data.mean);
      graphs[ind]->SetPointError(graphCounts[ind], 0.0, data.meanErr);
      if(data.mean - data.meanErr < ymin){ymin = data.mean - data.meanErr;}
      if(data.mean + data.meanErr > ymax){ymax = data.mean + data.meanErr;}
    }
    //if(runNum - 5 < xmin){xmin = runNum - 5;}
    //if(runNum + 5 > xmax){xmax = runNum + 5;}
    graphCounts[ind] += 1;
  }

  adjustGraphLimits(graphs, ymin, ymax, smallFac, largeFac, xmin, xmax);
  drawAndPrintGraphs(msmt, msmtNum, can, graphs);
  grand->Close();
}

#endif
