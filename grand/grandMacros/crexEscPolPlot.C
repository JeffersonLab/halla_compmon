#include "../grandOnline/makePlots.h"
#include "../grandConstruction/vars.h"


Int_t offset = 0;
Int_t nWiens = 0;
Int_t nIHWPs = 0;

TCanvas *can;
vector<vector<TGraphErrors *>> graphs;
TMultiGraph *mgraph;
vector<vector<Int_t>> counts;
vector<vector<TPaveText *>> texts;
TF1 *firstFit;
vector<vector<TF1 *>> escFits;
FitPolVar fitPolVar;
TString name = "Pol0";
TString title = "Polarization [pct]";

PolVar polVar;
Int_t cycNum, runNum, snailNum, cycCut;
Int_t snailTreeNum, nCycles, sign;
Float_t hWien, vWien, solWien, ihwp, laserPol;


vector<vector<Int_t>> getFitLimits(){
  vector<vector<Int_t>> limits;
  vector<Int_t> maxs;
  vector<Int_t> mins;
  vector<Int_t> wien;
  vector<Int_t> pols;

  mins.push_back(101); mins.push_back(138); mins.push_back(150); mins.push_back(190);
  maxs.push_back(137); maxs.push_back(150); maxs.push_back(185); maxs.push_back(223);
  wien.push_back(1); wien.push_back(0); wien.push_back(0); wien.push_back(1);
  pols.push_back(0); pols.push_back(1); pols.push_back(1); pols.push_back(0);

  limits.push_back(mins); limits.push_back(maxs);
  limits.push_back(wien); limits.push_back(pols);

  return limits;
}


vector<vector<vector<Float_t>>> getBoxLimits(){
  vector<vector<vector<Float_t>>> limits;
  vector<vector<Float_t>> xLo;
  vector<vector<Float_t>> xHi;
  vector<vector<Float_t>> yLo;
  vector<vector<Float_t>> yHi;
  vector<Float_t> xLoOut; vector<Float_t> xLoIn;
  vector<Float_t> xHiOut; vector<Float_t> xHiIn;
  vector<Float_t> yLoOut; vector<Float_t> yLoIn;
  vector<Float_t> yHiOut; vector<Float_t> yHiIn;

  xLoOut.push_back(0.10); xLoOut.push_back(0.30); xLoOut.push_back(0.50); xLoOut.push_back(0.70);
  xLoIn.push_back(0.10);  xLoIn.push_back(0.30);  xLoIn.push_back(0.50);  xLoIn.push_back(0.70);
  xHiOut.push_back(0.30); xHiOut.push_back(0.50); xHiOut.push_back(0.70); xHiOut.push_back(0.90);
  xHiIn.push_back(0.30);  xHiIn.push_back(0.50);  xHiIn.push_back(0.70);  xHiIn.push_back(0.90);
  yLoOut.push_back(0.80); yLoOut.push_back(0.80); yLoOut.push_back(0.80); yLoOut.push_back(0.80);
  yLoIn.push_back(0.70);  yLoIn.push_back(0.70);  yLoIn.push_back(0.70);  yLoIn.push_back(0.70);
  yHiOut.push_back(0.90); yHiOut.push_back(0.90); yHiOut.push_back(0.90); yHiOut.push_back(0.90);
  yHiIn.push_back(0.80);  yHiIn.push_back(0.80);  yHiIn.push_back(0.80);  yHiIn.push_back(0.80);

  xLo.push_back(xLoOut); xLo.push_back(xLoIn);
  xHi.push_back(xHiOut); xHi.push_back(xHiIn);
  yLo.push_back(yLoOut); yLo.push_back(yLoIn);
  yHi.push_back(yHiOut); yHi.push_back(yHiIn);

  limits.push_back(xLo); limits.push_back(xHi);
  limits.push_back(yLo); limits.push_back(yHi);

  return limits;
}


vector<vector<Int_t>> readEscargatoires(){
  Int_t startEsc = 301; Int_t endEsc = 395;

  vector<vector<Int_t>> allEscs;
  for(Int_t i = startEsc; i <= endEsc; i++){
    printf("Reading snails from file esc%i.list...\n", i);
    string snlNumStr;
    vector<Int_t> snailList;
    ifstream infile(Form("%s/escargatoires/esc%i.list", getenv("COMPMON_SNAILS"), i));
    if(infile.good()){
      while(getline(infile, snlNumStr)){
        Int_t snlNum = atoi(snlNumStr.c_str());
        if(snlNum == 0) continue;
        snailList.push_back(snlNum);
      }
      allEscs.push_back(snailList);
      snailList.clear();
    }
  }

  return allEscs;
}


vector<vector<Int_t>> getSlugNums(){
  Int_t startEsc = 301; Int_t endEsc = 395;

  vector<vector<Int_t>> allEscs;
  for(Int_t i = startEsc; i <= endEsc; i++){
    printf("Reading slugs from file escSlugs%i.list...\n", i);
    string slugNumStr;
    vector<Int_t> slugList;
    ifstream infile(Form("%s/escargatoires/escSlugs%i.list", getenv("COMPMON_SNAILS"), i));
    if(infile.good()){
      while(getline(infile, slugNumStr)){
        Int_t slugNum = atoi(slugNumStr.c_str());
        if(slugNum == 0) continue;
        slugList.push_back(slugNum);
      }
      allEscs.push_back(slugList);
      slugList.clear();
    }
  }

  return allEscs;
}


vector<Float_t> getAvgSlugNums(vector<vector<Int_t>> matchingSlugs){
  vector<Float_t> allCenters;
  for(Int_t i = 0; i < matchingSlugs.size(); i++){
    Int_t slugTotal = 0; Int_t nSlugs = 0;
    for(Int_t j = 0; j < matchingSlugs[i].size(); j++){
      Int_t slugNum = matchingSlugs[i][j];
      if(slugNum == 0) 
        continue;
      slugTotal += slugNum;
      nSlugs++;
    }

    Float_t avg = slugTotal*1.0/nSlugs;
    allCenters.push_back(avg);
  }

  return allCenters;
}


vector<Float_t> getPVMeanErrs(vector<vector<Int_t>> matchingSlugs){
  TFile *f = TFile::Open("/u/group/halla/parity/software/japan_offline/bmodAna/rootScripts/BeamMod/plotMacros/processed_respin2_data/results/CREX_All_IncludeBMOD_rcdb_slug_Avg_Outputs_main_det_corrections.root");
  TTree *maindet = (TTree *)f->Get("mini_eigen_lagr_allbpms_part_avg_det_asyms_det_weighted");

  vector<Float_t> allMeanErrs;
  for(Int_t i = 0; i < matchingSlugs.size(); i++){
    Float_t meanErrSum = 0.0;
    for(Int_t j = 0; j < matchingSlugs[i].size(); j++){
      printf("Making plots for slug %i...\n", matchingSlugs[i][j]);
      TString hName = Form("hSlug%i_meanErr", matchingSlugs[i][j]);
      maindet->Draw(Form("eigen_lagr_asym_manual_main_det_mean_err/1e-9>>%s", hName.Data()), Form("rcdb_slug==%i", matchingSlugs[i][j]), "goff");
      TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
      Float_t meanErr = (Float_t)h->GetMean();
      meanErrSum += 1.0/(meanErr*meanErr);
    }
    allMeanErrs.push_back(TMath::Sqrt(1.0/meanErrSum));
  }

  return allMeanErrs;
}


//==============================================================
//                                                            ||
//                      Main Functions                        ||
//                                                            ||
//==============================================================


void initOffset(){
  offset = 301;
  nWiens = 4;
  nIHWPs = 2;
}


void initBranches(TTree *snl, TTree *cyc){
  cyc->SetBranchAddress("snailNum", &snailNum);
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress(name.Data(), &polVar);
  cyc->SetBranchAddress("sign", &sign);
  cyc->SetBranchAddress("CycleCut", &cycCut);
  snl->SetBranchAddress("snailNum", &snailTreeNum);
  snl->SetBranchAddress("numCyclesAcc", &nCycles);
  snl->SetBranchAddress("HWienAngle", &hWien);
  snl->SetBranchAddress("VWienAngle", &vWien);
  snl->SetBranchAddress("PhiFG", &solWien);
  snl->SetBranchAddress("ihwp", &ihwp);
  snl->SetBranchAddress("LaserPolarization", &laserPol);
  snl->SetBranchAddress("sign", &sign);
}


void initFits(vector<vector<Int_t>> fitLimits, vector<vector<vector<Float_t>>> boxLimits){
  can = new TCanvas(Form("c%s", name.Data()), Form("%s Canvas", title.Data()), 1200, 800);
  can->SetGridx(); can->SetGridy();
  firstFit = new TF1(Form("firstFit_%s", name.Data()), "pol0");

  mgraph = new TMultiGraph();
  mgraph->SetTitle(Form("%s vs Escargatoire", title.Data()));
  mgraph->GetXaxis()->SetTitle("Avg Slug Num");
  mgraph->GetYaxis()->SetTitle(title.Data());

  for(Int_t i = 0; i < nWiens; i++){
    vector<TGraphErrors *> innerGraphs;
    vector<Int_t> innerCounts;
    vector<TF1 *> innerFits;
    vector<TPaveText *> innerTexts;
  
    Int_t colorIndBase = (fitLimits[2][i] == 0) ? 0 : 2;

    for(Int_t j = 0; j < nIHWPs; j++){
      Int_t colorInd = colorIndBase + j;
      TGraphErrors *thisGraph = new TGraphErrors();
      thisGraph->SetMarkerColor(getColor(colorInd));
      thisGraph->SetMarkerStyle(getMarkerStyle(colorInd));

      TF1 *thisEscFit = new TF1(Form("escFit_%s_wien%i_ihwp%i", name.Data(), i, j), Form("pol%i", fitLimits[3][i]));
      thisEscFit->SetLineColor(getColor(colorInd));

      TPaveText *thisText = new TPaveText(boxLimits[0][j][i], boxLimits[2][j][i], boxLimits[1][j][i], boxLimits[3][j][i], "blNDC");
      thisText->SetFillColor(0); thisText->SetBorderSize(1);
      
      innerGraphs.push_back(thisGraph);
      innerCounts.push_back(0);
      innerFits.push_back(thisEscFit);
      innerTexts.push_back(thisText);
    }

    graphs.push_back(innerGraphs);
    counts.push_back(innerCounts);
    escFits.push_back(innerFits);
    texts.push_back(innerTexts);
  }
}


Float_t uncertRatio(vector<Float_t> avgSlug, vector<Float_t> meanErrs, vector<Float_t> polErrs){
  Float_t numSum = 0.0;
  Float_t demSum = 0.0;
  for(Int_t i = 0; i < avgSlug.size(); i++){
    numSum += 1.0/(meanErrs[i]*meanErrs[i]*meanErrs[i]*meanErrs[i])*polErrs[i]*polErrs[i];
    demSum += 1.0/(meanErrs[i]*meanErrs[i]);
  }

  return TMath::Sqrt(numSum*1.0/demSum);
}


void crexEscPolPlot(){
  TString expt = experimentCode(2);

  vector<vector<Int_t>> allEscs = readEscargatoires();
  vector<vector<Int_t>> fitLimits = getFitLimits();
  vector<vector<vector<Float_t>>> boxLimits = getBoxLimits();
  vector<vector<Int_t>> allSlugs = getSlugNums();
  vector<Float_t> avgSlugs = getAvgSlugNums(allSlugs);
  vector<Float_t> meanErrs = getPVMeanErrs(allSlugs);
  vector<Float_t> polErrs;

  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");
  TTree *snl = (TTree *)f->Get("snl");

  initOffset();
  initBranches(snl, cyc);
  initFits(fitLimits, boxLimits);

  for(Int_t i = 0; i < allEscs.size(); i++){
    Float_t escHWien = 0.0; Float_t escVWien = 0.0; Float_t escSolWien = 0.0; 
    Float_t escIHWP = -1.0; Float_t escLaserPol = 0.0; 
    Int_t escCycles = 0; Int_t escSign = 0;
    for(Int_t j = 0; j < snl->GetEntries(); j++){
      snl->GetEntry(j);
      Bool_t escMatch = kFALSE;
      for(Int_t k = 0; k < allEscs[i].size(); k++){
        escMatch = escMatch || snailTreeNum == allEscs[i][k];
      }
      if(escMatch){
        escHWien = hWien; escVWien = vWien; escSolWien = solWien;
        escIHWP = ihwp; escLaserPol = laserPol; escCycles += nCycles;
        escSign = sign;
      }
    }

    TH1F *hist = new TH1F(Form("h%s_esc%i", name.Data(), i+offset), "", escCycles, 0, escCycles);;
    Int_t cycCount = 1;

    for(Int_t j = 0; j < cyc->GetEntries(); j++){
      cyc->GetEntry(j);
      Bool_t escMatch = kFALSE;
      for(Int_t k = 0; k < allEscs[i].size(); k++){
        escMatch = escMatch || snailNum == allEscs[i][k];
      }
      if(escMatch && cycCut == 0){
        hist->SetBinContent(cycCount, polVar.mean);
        hist->SetBinError(cycCount++, polVar.meanErr);
      }
    }
  
    hist->Fit(Form("firstFit_%s", name.Data()), "Q", "", 0, escCycles);

    Int_t wienInd = -1;
    Int_t ihwpInd = escIHWP;
    for(Int_t j = 0; j < nWiens; j++){
      if(avgSlugs[i] > 186.0 && avgSlugs[i] < 190.0){
        wienInd = 3;
        break;
      }
      if(avgSlugs[i] >= fitLimits[0][j] && avgSlugs[i] <= fitLimits[1][j]){
        //if(fitLimits[3][j] == 1) ihwpInd = 0;
        wienInd = j;
        break;
      }
    }
    Float_t yval = 100*escSign*firstFit->GetParameter(0);
    //if(wienInd == 1){yval -= -0.106*avgSlugs[i];}
    //if(wienInd == 2){yval -= -0.025*avgSlugs[i];}
    graphs[wienInd][ihwpInd]->SetPoint(counts[wienInd][ihwpInd], avgSlugs[i], yval);
    graphs[wienInd][ihwpInd]->SetPointError(counts[wienInd][ihwpInd]++, 0.0, 100*firstFit->GetParError(0));   

    //printf("Escargatoire %i gets Pol0 %.4f +/- %.4f\n", i+offset, 100*escSign*fits[j]->GetParameter(0), 100*fits[j]->GetParError(0));
    printf("%i,%.3f,%.4f,%.4f,%f\n", i+offset, avgSlugs[i], 100*escSign*firstFit->GetParameter(0), 100*firstFit->GetParError(0), meanErrs[i]);
    polErrs.push_back(100*firstFit->GetParError(0));
  }

  can->cd();
  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      mgraph->Add(graphs[i][j]);
    }
  }
  mgraph->Draw("ap");
  
  /**
  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      printf("GraphInd: (%i, %i) and fitting between %i-%i\n", i, j, fitLimits[0][i], fitLimits[1][i]);
      Int_t colorIndBase = (fitLimits[2][i] == 0) ? 0 : 2;
      graphs[i][j]->Fit(escFits[i][j], "Q", "", fitLimits[0][i], fitLimits[1][i]);
      TText *lineOut1 = texts[i][j]->AddText(Form("------%s------\n", getLegendEntry(colorIndBase + j).Data()));
      lineOut1->SetTextColor(getColor(colorIndBase + j));
      TText *lineOut2 = texts[i][j]->AddText(Form("Par0: %.3f +/- %.3f\n", escFits[i][j]->GetParameter(0), escFits[i][j]->GetParError(0)));
      lineOut2->SetTextColor(getColor(colorIndBase + j));
      if(fitLimits[3][i] == 1){
        TText *lineOut3 = texts[i][j]->AddText(Form("Par1: %.3f +/- %.3f\n", escFits[i][j]->GetParameter(1), escFits[i][j]->GetParError(1)));
        lineOut3->SetTextColor(getColor(colorIndBase + j));
      }
      //texts[i][j]->Draw("same");

      //escFits[i][j]->Draw("same");
    }
  }
  **/

  /**
  Float_t par0RelErr = escFits[2][0]->GetParError(0)*1.0/escFits[2][0]->GetParameter(0);
  Float_t par1RelErr = escFits[2][0]->GetParError(1)*1.0/escFits[2][0]->GetParameter(1);
  TF1 *loFit = new TF1("loFit", "pol1", fitLimits[0][2]-168.92, fitLimits[1][2]-168.92);
  loFit->SetParameter(0, escFits[2][0]->GetParameter(0)*(1 - par0RelErr));
  loFit->SetParameter(1, escFits[2][0]->GetParameter(1)*(1 - par1RelErr));
  loFit->SetLineColor(kGreen + 2);
  loFit->Draw("same");
  
  TF1 *hiFit = new TF1("hiFit", "pol1", fitLimits[0][2]-168.92, fitLimits[1][2]-168.92);
  hiFit->SetParameter(0, escFits[2][0]->GetParameter(0)*(1 + par0RelErr));
  hiFit->SetParameter(1, escFits[2][0]->GetParameter(1)*(1 + par1RelErr));
  hiFit->SetLineColor(kCyan + 2);
  hiFit->Draw("same");

  printf("Par 0 Limits: %f - %f - %f\n", escFits[2][0]->GetParameter(0)*(1-par0RelErr), escFits[2][0]->GetParameter(0), escFits[2][0]->GetParameter(0)*(1+par0RelErr));
  printf("Par 1 Limits: %f - %f - %f\n", escFits[2][0]->GetParameter(1)*(1-par1RelErr), escFits[2][0]->GetParameter(1), escFits[2][0]->GetParameter(1)*(1+par1RelErr));
  **/

  Float_t errRatio = uncertRatio(avgSlugs, meanErrs, polErrs);
  printf("Uncertainty Ratio is %f\n", errRatio);
}
