#ifndef RUNWISEAGGREGATE_H
#define RUNWISEAGGREGATE_H

#include "../grandConstruction/vars.h"

Int_t snailNum, runNum, numCycles, numCyclesAcc;
Float_t runTime;

const Int_t nPols = 5;
TString polNames0[nPols] = {"Pol0", "Asym0", "Asym0NGC", "Asym0LasOn", "Asym0LasOff"};
Float_t pol0Mults[nPols] = {100.0, 1000.0, 1000.0, 1000.0, 1000.0};
vector<FitPolVar> pols0; vector<TF1 *> fits0;
vector<TH1F *> pol0Hists; vector<TH1F *> pol0Pulls;
vector<TH1F *> pol0Pull2Ds; vector<TPaveText *> pol0Texts;
vector<vector<TLine *>> pol0Lines;
vector<Float_t> pol0Min, pol0Max;
// TString polNames4[nPols] = {"Pol4", "Asym4", "Asym4NGC", "Asym4LasOn", "Asym4LasOff"};
// Float_t pol4Mults[nPols] = {100.0, 1000.0, 1000.0, 1000.0, 1000.0};
// vector<FitPolVar> pols4; vector<TF1 *> fits4;
// vector<TH1F *> pol4Hists; vector<TH1F *> pol4Pulls;
// vector<TH1F *> pol4Pull2Ds; vector<TPaveText *> pol4Texts;
// vector<vector<TLine *>> pol4Lines;
// vector<Float_t> pol4Min, pol4Max;

const Int_t nAccVars = 3;
TString varNames0[nAccVars] = {"Acc0LasOn", "Acc0LasOff", "Acc0BeamOff"};
vector<DataVar> vars0;
vector<TH1F *> var0Hists;
vector<vector<TLine *>> var0Lines;
vector<Float_t> var0Min, var0Max;
// TString varNames4[nAccVars] = {"Acc4LasOn", "Acc4LasOff", "Acc4BeamOff"};
// vector<DataVar> vars4;
// vector<TH1F *> var4Hists;
// vector<vector<TLine *>> var4Lines;
// vector<Float_t> var4Min, var4Max;

const Int_t nVars = 16;
TString varNames[nVars] = {"LaserPower", "BeamCurrent", "bpmAx", "bpmAy", "bpmBx", "bpmBy", 
                           "CentralRateLasOn", "CentralRateLasOff", "USbg1", "USbg2","DSbg1", "DSbg2", 
                           "HFingerLasOn", "HFingerLasOff", "VFingerLasOn", "VFingerLasOff"};
vector<DataVar> vars;
vector<TH1F *> varHists;
vector<vector<TLine *>> varLines;
vector<Float_t> varMin, varMax;

// const Int_t nBursts = 20;
// TString burstNames[nBursts] = {"BurstPosAcc0LasOn",   "BurstNegAcc0LasOn",   "BurstDiffAcc0LasOn",   "BurstSummAcc0LasOn",   "BurstAsym0LasOn",
//                                "BurstPosAcc0LasOff1", "BurstNegAcc0LasOff1", "BurstDiffAcc0LasOff1", "BurstSummAcc0LasOff1", "BurstAsym0LasOff1",
//                                "BurstPosAcc0LasOff2", "BurstNegAcc0LasOff2", "BurstDiffAcc0LasOff2", "BurstSummAcc0LasOff2", "BurstAsym0LasOff2",
//                                "BurstPosAcc0LasOff",  "BurstNegAcc0LasOff",  "BurstDiffAcc0LasOff",  "BurstSummAcc0LasOff",  "BurstAsym0LasOff"};
// vector<FitPolVar> burstVars;
// vector<TH1F *> burstHists;
// vector<vector<TLine *>> burstLines;
// vector<Float_t> burstMin, burstMax;

const Int_t nCombos = 5;
TString comboNames[nCombos] = {"BurstAsym0LasOn", "BurstAsym0LasOff", "BurstAsym0NGC", "BurstAsym0", "BurstPol0"};
Float_t comboMults[nCombos] = {1000.0, 1000.0, 1000.0, 1000.0, 100.0};
vector<FitPolVar> comboVars; vector<TF1 *> comboFits;
vector<TH1F *> comboHists; vector<TH1F *> comboPulls;
vector<TH1F *> comboPull2Ds; vector<TPaveText *> comboTexts;
vector<vector<TLine *>> comboLines;
vector<Float_t> comboMin, comboMax;

vector<Bool_t> isPol;

Int_t nRuns;

TString exptName(Int_t snailNum){
  if(snailNum < 100 || snailNum == 500) return "prex";
  else return "crex";
}

void initHists(Int_t snlNum){
  for(Int_t i = 0; i < nPols; i++){
    TH1F *h; TH1F *hPull2D;
    h = new TH1F(Form("h_%s", polNames0[i].Data()), Form("Snail %i: %s vs Run", snlNum, polNames0[i].Data()), nRuns, 0, nRuns);
    TH1F *hPull = new TH1F(Form("hPull_%s", polNames0[i].Data()), Form("%s Pull Plot", polNames0[i].Data()), 40, -8, 8);
    hPull2D = new TH1F(Form("hPull2D_%s", polNames0[i].Data()), "", nRuns, 0, nRuns);
    hPull2D->SetLineColor(kGreen); hPull2D->SetFillColor(kGreen);
    pol0Hists.push_back(h); pol0Pulls.push_back(hPull); pol0Pull2Ds.push_back(hPull2D);
    vector<TLine *> lines;
    pol0Lines.push_back(lines);
    pol0Min.push_back(1e16); pol0Max.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString hName = Form("Snail %i: %s vs Run", snlNum, varNames0[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail%i: %s RMS vs Run", snlNum, varNames0[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames0[(Int_t)i/2].Data()), hName.Data(), nRuns, 0, nRuns);
    var0Hists.push_back(h);
    vector<TLine *> lines;
    var0Lines.push_back(lines);
    var0Min.push_back(1e16); var0Max.push_back(-1e16);
  }
  /**
  for(Int_t i = 0; i < nPols; i++){
    TH1F *h; TH1F *hPull2D;
    h = new TH1F(Form("h_%s", polNames4[i].Data()), Form("Snail %i: %s vs Run", snlNum, polNames4[i].Data()), nRuns, 0, nRuns);
    TH1F *hPull = new TH1F(Form("hPull_%s", polNames4[i].Data()), Form("%s Pull Plot", polNames4[i].Data()), 40, -8, 8);
    hPull2D = new TH1F(Form("hPull2D_%s", polNames4[i].Data()), "", nRuns, 0, nRuns);
    hPull2D->SetLineColor(kGreen); hPull2D->SetFillColor(kGreen);
    pol4Hists.push_back(h); pol4Pulls.push_back(hPull); pol4Pull2Ds.push_back(hPull2D);
    vector<TLine *> lines;
    pol4Lines.push_back(lines);
    pol4Min.push_back(1e16); pol4Max.push_back(-1e16);
  }
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString hName = Form("Snail %i: %s vs Run", snlNum, varNames4[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s RMS vs Run", snlNum, varNames4[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames4[(Int_t)i/2].Data()), hName.Data(), nRuns, 0, nRuns);
    var4Hists.push_back(h);
    vector<TLine *> lines;
    var4Lines.push_back(lines);
    var4Min.push_back(1e16); var4Max.push_back(-1e16);
  }
  **/
  for(Int_t i = 0; i < 2*nVars; i++){
    TString hName = Form("Snail %i: %s vs Run", snlNum, varNames[(Int_t)i/2].Data());
    TString hNameAdd("");
    if(i % 2 == 1){hName = Form("Snail %i: %s RMS vs Run", snlNum, varNames[(Int_t)i/2].Data()); hNameAdd = "RMS";}
    TH1F *h = new TH1F(Form("h%s_%s", hNameAdd.Data(), varNames[(Int_t)i/2].Data()), hName.Data(), nRuns, 0, nRuns);
    varHists.push_back(h);
    vector<TLine *> lines;
    varLines.push_back(lines);
    varMin.push_back(1e16); varMax.push_back(-1e16);
  }
  for(Int_t i = 0; i < nCombos; i++){
    TH1F *h; TH1F *hPull2D;
    h = new TH1F(Form("h_%s", comboNames[i].Data()), Form("Snail %i: %s vs Runs", snlNum, comboNames[i].Data()), nRuns, 0, nRuns);
    TH1F *hPull = new TH1F(Form("hPull_%s", comboNames[i].Data()), Form("%s Pull Plot", comboNames[i].Data()), 40, -8, 8);
    hPull2D = new TH1F(Form("hPull2D_%s", comboNames[i].Data()), "", nRuns, 0, nRuns);
    hPull2D->SetLineColor(kGreen); hPull2D->SetFillColor(kGreen);
    comboHists.push_back(h); comboPulls.push_back(hPull); comboPull2Ds.push_back(hPull2D);
    vector<TLine *> lines;
    comboLines.push_back(lines);
    comboMin.push_back(1e16); comboMax.push_back(-1e16);
  }
}

void initFits(){
  for(Int_t i = 0; i < nPols; i++){
    TString fName0 = Form("f_%s", polNames0[i].Data());
    // TString fName4 = Form("f_%s", polNames4[i].Data());
    TF1 *f0 = new TF1(fName0.Data(), "pol0");
    // TF1 *f4 = new TF1(fName4.Data(), "pol0");
    fits0.push_back(f0);
    // fits4.push_back(f4);

    TPaveText *pt0 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
    // TPaveText *pt4 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
    pol0Texts.push_back(pt0);
    // pol4Texts.push_back(pt4);
  }
  for(Int_t i = 0; i < nCombos; i++){
    TString fName0 = Form("f_%s", comboNames[i].Data());
    TF1 *f0 = new TF1(fName0.Data(), "pol0");
    comboFits.push_back(f0);

    TPaveText *pt0 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
    comboTexts.push_back(pt0);
  }
}

void initBranches(TTree* run){
  run->SetBranchAddress("snailNum", &snailNum);
  run->SetBranchAddress("runNum", &runNum);
  run->SetBranchAddress("numCycles", &numCycles);
  run->SetBranchAddress("numCyclesAcc", &numCyclesAcc);
  run->SetBranchAddress("runTime", &runTime);

  for(Int_t i = 0; i < nPols; i++){FitPolVar pol; pols0.push_back(pol);}
  for(Int_t i = 0; i < nAccVars; i++){DataVar data; vars0.push_back(data);}
  // for(Int_t i = 0; i < nPols; i++){FitPolVar pol; pols4.push_back(pol);}
  // for(Int_t i = 0; i < nAccVars; i++){DataVar data; vars4.push_back(data);}
  for(Int_t i = 0; i < nVars; i++){DataVar data; vars.push_back(data);}
  for(Int_t i = 0; i < nCombos; i++){FitPolVar data; comboVars.push_back(data);}

  for(Int_t i = 0; i < nPols; i++){run->SetBranchAddress(polNames0[i].Data(), &pols0[i]);}
  for(Int_t i = 0; i < nAccVars; i++){run->SetBranchAddress(varNames0[i].Data(), &vars0[i]);}
  // for(Int_t i = 0; i < nPols; i++){run->SetBranchAddress(polNames4[i].Data(), &pols4[i]);}
  // for(Int_t i = 0; i < nAccVars; i++){run->SetBranchAddress(varNames4[i].Data(), &vars4[i]);}
  for(Int_t i = 0; i < nVars; i++){run->SetBranchAddress(varNames[i].Data(), &vars[i]);}
  for(Int_t i = 0; i < nCombos; i++){run->SetBranchAddress(comboNames[i].Data(), &comboVars[i]);}
}

void setPol0Strs(Int_t i){
  Float_t par = fits0[i]->GetParameter(0);
  Float_t parErr = fits0[i]->GetParError(0);
  Float_t chi2 = fits0[i]->GetChisquare();
  Int_t ndf = fits0[i]->GetNDF();

  pol0Texts[i]->AddText(Form("--------Fit Results--------"));
  pol0Texts[i]->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  pol0Texts[i]->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  pol0Texts[i]->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  pol0Texts[i]->SetBorderSize(1); pol0Texts[i]->SetFillColor(0);
}

/**
void setPol4Strs(Int_t i){
  Float_t par = fits4[i]->GetParameter(0);
  Float_t parErr = fits4[i]->GetParError(0);
  Float_t chi2 = fits4[i]->GetChisquare();
  Int_t ndf = fits4[i]->GetNDF();

  pol4Texts[i]->AddText(Form("--------Fit Results--------"));
  pol4Texts[i]->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  pol4Texts[i]->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  pol4Texts[i]->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  pol4Texts[i]->SetBorderSize(1); pol4Texts[i]->SetFillColor(0);
}
**/

void setComboStrs(Int_t i){
  Float_t par = comboFits[i]->GetParameter(0);
  Float_t parErr = comboFits[i]->GetParError(0);
  Float_t chi2 = comboFits[i]->GetChisquare();
  Int_t ndf = comboFits[i]->GetNDF();

  comboTexts[i]->AddText(Form("--------Fit Results--------"));
  comboTexts[i]->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  comboTexts[i]->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  comboTexts[i]->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  comboTexts[i]->SetBorderSize(1); comboTexts[i]->SetFillColor(0);
}

void makePol0Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < nPols; i++){
    TCanvas *c = new TCanvas(Form("cPol_%s", polNames0[i].Data()), "Pol Canvas", 1200, 600);
    TPad *pPol1 = new TPad(Form("pPol1_%s", polNames0[i].Data()), "Pol Avg", 0.0, 0.3, 0.7, 1.0);
    pPol1->SetGridx(1); pPol1->SetGridy(1);
    TPad *pPol2 = new TPad(Form("pPol2_%s", polNames0[i].Data()), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
    TPad *pPol3 = new TPad(Form("pPol3_%s", polNames0[i].Data()), "Pol Pull Graph", 0.0, 0.0, 0.7, 0.3);
    pPol3->SetGridx(1); pPol3->SetGridy(1);
    pPol1->Draw(); pPol2->Draw(); pPol3->Draw();

    pPol1->cd();
    pol0Hists[i]->Fit(fits0[i], "Q", "", 0, nRuns);
    Float_t par = fits0[i]->GetParameter(0);
    for(Int_t j = 0; j < nRuns; j++){
      Float_t nDiffs = (pol0Hists[i]->GetBinContent(j + 1) - par)/pol0Hists[i]->GetBinError(j + 1);
      pol0Pull2Ds[i]->SetBinContent(j + 1, nDiffs); pol0Pulls[i]->Fill(nDiffs);
    }
    
    setPol0Strs(i);
    pol0Hists[i]->SetStats(0);
    pol0Hists[i]->Draw("P");
    for(Int_t j = 0; j < pol0Lines[i].size(); j++){
      pol0Lines[i][j]->Draw("same");
    }
    pol0Texts[i]->Draw("same");

    pPol2->cd();
    pol0Pulls[i]->SetStats(220);
    pol0Pulls[i]->Draw();
    
    pPol3->cd();
    pol0Pull2Ds[i]->SetStats(0);
    pol0Pull2Ds[i]->Draw();

    c->Print(Form("%s/run_agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
  //return msmt;
}

void makePol0Lines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < pol0Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      Float_t xpos = bounds[j];
      TLine *line = new TLine(xpos, pol0Min[i], xpos, pol0Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      pol0Lines[i].push_back(line);
    }
  }
}

/**
void makePol4Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < nPols; i++){
    TCanvas *c = new TCanvas(Form("cPol_%s", polNames4[i].Data()), "Pol Canvas", 1200, 600);
    TPad *pPol1 = new TPad(Form("pPol1_%s", polNames4[i].Data()), "Pol Avg", 0.0, 0.3, 0.7, 1.0);
    pPol1->SetGridx(1); pPol1->SetGridy(1);
    TPad *pPol2 = new TPad(Form("pPol2_%s", polNames4[i].Data()), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
    TPad *pPol3 = new TPad(Form("pPol3_%s", polNames4[i].Data()), "Pol Pull Graph", 0.0, 0.0, 0.7, 0.3);
    pPol3->SetGridx(1); pPol3->SetGridy(1);
    pPol1->Draw(); pPol2->Draw(); pPol3->Draw();

    pPol1->cd();
    pol4Hists[i]->Fit(fits4[i], "Q", "", 0, nRuns);
    Float_t par = fits4[i]->GetParameter(0);
    for(Int_t j = 0; j < nRuns; j++){
      Float_t nDiffs = (pol4Hists[i]->GetBinContent(j + 1) - par)/pol4Hists[i]->GetBinError(j + 1);
      pol4Pull2Ds[i]->SetBinContent(j + 1, nDiffs); pol4Pulls[i]->Fill(nDiffs);
    }
    
    setComboStrs(i);
    pol4Hists[i]->SetStats(0);
    pol4Hists[i]->Draw("P");
    for(Int_t j = 0; j < pol4Lines[i].size(); j++){
      pol4Lines[i][j]->Draw("same");
    }
    pol4Texts[i]->Draw("same");

    pPol2->cd();
    pol4Pulls[i]->SetStats(220);
    pol4Pulls[i]->Draw();
    
    pPol3->cd();
    pol4Pull2Ds[i]->SetStats(0);
    pol4Pull2Ds[i]->Draw();

    c->Print(Form("%s/run_agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
  //return msmt;
}

void makePol4Lines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < pol4Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      Float_t xpos = bounds[j];
      TLine *line = new TLine(xpos, pol4Min[i], xpos, pol4Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      pol4Lines[i].push_back(line);
    }
  }
}
**/

void makeVar0Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString canName = Form("cAccVar_%s", varNames0[(Int_t)i/2].Data());
    if(i % 2 == 1)
      canName = Form("cAccVar_%s_RMS", varNames0[(Int_t)i/2].Data());
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Acc Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    var0Hists[i]->SetStats(0);
    var0Hists[i]->Draw("P");
    for(Int_t j = 0; j < var0Lines[i].size(); j++){
      var0Lines[i][j]->Draw("same");
    }
    c->Print(Form("%s/run_agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeVar0Lines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < var0Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], var0Min[i], bounds[j], var0Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      var0Lines[i].push_back(line);
    }
  }
}

/**
void makeVar4Plots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nAccVars; i++){
    TString canName = Form("cAccVar_%s", varNames4[(Int_t)i/2].Data());
    if(i % 2 == 1)
      canName = Form("cAccVar_%s_RMS", varNames4[(Int_t)i/2].Data());
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Acc Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    var4Hists[i]->SetStats(0);
    var4Hists[i]->Draw("P");
    for(Int_t j = 0; j < var4Lines[i].size(); j++){
      var4Lines[i][j]->Draw("same");
    }
    c->Print(Form("%s/run_agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeVar4Lines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < var4Hists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], var4Min[i], bounds[j], var4Max[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      var4Lines[i].push_back(line);
    }
  }
}
**/

void makeVarPlots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < 2*nVars; i++){
    TString canName = Form("cAccVar_%s", varNames[(Int_t)i/2].Data());
    if(i % 2 == 1)
      canName = Form("cAccVar_%s_RMS", varNames[(Int_t)i/2].Data());
    printf("MSMT #%i: %s\n", *msmt, canName.Data());
    TCanvas *c = new TCanvas(canName.Data(), "Var Canvas", 1200, 400);
    c->SetGridx(1); c->SetGridy(1);
    varHists[i]->SetStats(0);
    varHists[i]->Draw("P");
    for(Int_t j = 0; j < varLines[i].size(); j++){
      varLines[i][j]->Draw("same");
    }
    c->Print(Form("%s/run_agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
}

void makeVarLines(Int_t snlNum, vector<Float_t> bounds){
  for(Int_t i = 0; i < varHists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      TLine *line = new TLine(bounds[j], varMin[i], bounds[j], varMax[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      varLines[i].push_back(line);
    }
  }
}

void makeComboPlots(Int_t snlNum, Int_t *msmt){
  TString snailName = Form("snail%i", snlNum);
  TString outputDir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snailName.Data()));
  for(Int_t i = 0; i < nCombos; i++){
    TCanvas *c = new TCanvas(Form("cPol_%s", comboNames[(Int_t)i/2].Data()), "Pol Canvas", 1200, 600);
    TPad *pPol1 = new TPad(Form("pPol1_%s", comboNames[(Int_t)i/2].Data()), "Pol Avg", 0.0, 0.3, 0.7, 1.0);
    pPol1->SetGridx(1); pPol1->SetGridy(1);
    TPad *pPol2 = new TPad(Form("pPol2_%s", comboNames[(Int_t)i/2].Data()), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
    TPad *pPol3 = new TPad(Form("pPol3_%s", comboNames[(Int_t)i/2].Data()), "Pol Pull Graph", 0.0, 0.0, 0.7, 0.3);
    pPol3->SetGridx(1); pPol3->SetGridy(1);
    pPol1->Draw(); pPol2->Draw(); pPol3->Draw();

    pPol1->cd();
    comboHists[i]->Fit(comboFits[i], "Q", "", 0, nRuns);
    Float_t par = comboFits[i]->GetParameter(0);
    for(Int_t j = 0; j < nRuns; j++){
      Float_t nDiffs = (comboHists[i]->GetBinContent(j + 1) - par)/comboHists[i]->GetBinError(j + 1);
      comboPull2Ds[i]->SetBinContent(j + 1, nDiffs); comboPulls[i]->Fill(nDiffs);
    }
    
    setComboStrs(i);
    comboHists[i]->SetStats(0);
    comboHists[i]->Draw("P");
    for(Int_t j = 0; j < comboLines[i].size(); j++){
      comboLines[i][j]->Draw("same");
    }
    comboTexts[i]->Draw("same");

    pPol2->cd();
    comboPulls[i]->SetStats(220);
    comboPulls[i]->Draw();
    
    pPol3->cd();
    comboPull2Ds[i]->SetStats(0);
    comboPull2Ds[i]->Draw();

    c->Print(Form("%s/run_agg_plots_%04i.pdf", outputDir.Data(), (*msmt)++), "pdf");
  }
  //return msmt;
}

void makeComboLines(Int_t snlNum, vector<Float_t> bounds, vector<Float_t> cutBounds){
  for(Int_t i = 0; i < comboHists.size(); i++){
    for(Int_t j = 0; j < bounds.size(); j++){
      Float_t xpos = (i % 2 == 1) ? bounds[j] : cutBounds[j];
      TLine *line = new TLine(xpos, comboMin[i], xpos, comboMax[i]);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      comboLines[i].push_back(line);
    }
  }
}

#endif
