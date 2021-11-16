#include "../grandOnline/makePlots.h"

using namespace std;

const Int_t nWiens = 2;
const Int_t nIHWPs = 2;
Int_t escLo = 301;
Int_t escHi = 395;

TGraphErrors *gCompton[nWiens][nIHWPs];
TGraphErrors *gMoller[nWiens][nIHWPs];
Int_t nPtsCompton[nWiens][nIHWPs]; 
Int_t nPtsMoller[nWiens][nIHWPs];

Int_t markerStyles[nWiens][nIHWPs] = {{23, 22}, {20, 21}};
Int_t colors[nWiens][nIHWPs] = {{kViolet, kOrange}, {kRed, kBlue}};
TString legEntry[nWiens][nIHWPs] = {{"Left Out", "Left In"}, {"Right Out", "Right In"}};

Int_t mollerMarkerStyles[nWiens][nIHWPs] = {{47, 29}, {33, 34}};
Int_t mollerColors[nWiens][nIHWPs] = {{kTeal, kOrange+2}, {kGreen+1, kMagenta-4}};

vector<Float_t> avgSlugs;
vector<vector<Double_t>> polData;
vector<Int_t> signs;
vector<Int_t> ihwps;
vector<Int_t> wiens;

vector<vector<vector<Float_t>>> molSlugAvgs;
vector<vector<vector<Float_t>>> molPols;
vector<vector<vector<Float_t>>> molErrs;
vector<vector<vector<Int_t>>> molNums;

const Int_t nFuncs = 8;
Float_t par0[nFuncs] = {86.682, 86.735, 101.595, 101.552, 91.256, 91.109, 87.784, 87.413};
Float_t par1[nFuncs] = {0.0, 0.0, -0.106, -0.106, -0.025, -0.025, 0.0, 0.0};
Float_t fMins[nFuncs] = {100.0, 100.0, 138.0, 138.0, 150.0, 150.0, 190.0, 190.0};
Float_t fMaxs[nFuncs] = {137.0, 137.0, 150.0, 150.0, 185.0, 185.0, 222.0, 222.0};
Int_t fCols[nFuncs] = {kRed, kBlue, kViolet, kOrange, kViolet, kOrange, kRed, kBlue};
TF1 *funcs[nFuncs];


void initFunctions(){
  for(Int_t i = 0; i < nFuncs; i++){
    TF1 *f = new TF1(Form("f%i", i), "pol1", fMins[i], fMaxs[i]);
    f->SetParameter(0, par0[i]);
    f->SetParameter(1, par1[i]);
    f->SetLineColor(fCols[i]);

    funcs[i] = f;
  }
}


//=====================================================================================\\
//                                                                                     ||
//                                                                                     ||
//                                   Add Moller Data                                   ||
//                                                                                     ||
//                                                                                     ||
//=====================================================================================//


void initMollerVectors(){
  for(Int_t i = 0; i < nWiens; i++){
    vector<vector<Float_t>> molWienAvg;
    vector<vector<Float_t>> molWienPol;
    vector<vector<Float_t>> molWienErr;
    vector<vector<Int_t>> molWienNum;
    for(Int_t j = 0; j < nIHWPs; j++){
      vector<Float_t> molIHWPAvg;
      vector<Float_t> molIHWPPol;
      vector<Float_t> molIHWPErr;
      vector<Int_t> molIHWPNum;

      molWienAvg.push_back(molIHWPAvg);
      molWienPol.push_back(molIHWPPol);
      molWienErr.push_back(molIHWPErr);
      molWienNum.push_back(molIHWPNum);
    }

    molSlugAvgs.push_back(molWienAvg);
    molPols.push_back(molWienPol);
    molErrs.push_back(molWienErr);
    molNums.push_back(molWienNum);
  }
}


void fillMollerVectors(){
  initMollerVectors();

  const Int_t nLeftOut = 2;
  const Int_t nLeftIn = 1;
  const Int_t nRightOut = 7;
  const Int_t nRightIn = 7;

  Float_t avgLeftOut[nLeftOut] = {144.0, 175.0};
  Float_t polLeftOut[nLeftOut] = {86.79, 87.08};
  Float_t errLeftOut[nLeftOut] = {0.16, 0.23};
  Float_t numLeftOut[nLeftOut] = {3039, 3042};

  Float_t avgLeftIn[nLeftIn] = {144.0};
  Float_t polLeftIn[nLeftIn] = {86.23};
  Float_t errLeftIn[nLeftIn] = {0.16};
  Float_t numLeftIn[nLeftIn] = {3038};

  Float_t avgRightOut[nRightOut] = {100.0, 105.5, 120.0, 133.5, 194.5, 213.0, 222.0};
  Float_t polRightOut[nRightOut] = {87.02, 86.85, 87.12, 86.98, 87.68, 87.45, 87.62};
  Float_t errRightOut[nRightOut] = {0.20, 0.19, 0.51, 0.20, 0.20, 0.16, 0.18};
  Float_t numRightOut[nRightOut] = {3002, 3013, 3019, 3022, 3046, 3048, 3052};

  Float_t avgRightIn[nRightIn] = {100.0, 105.5, 120.0, 133.5, 194.5, 213.0, 222.0};
  Float_t polRightIn[nRightIn] = {86.65, 87.01, 86.39, 86.83, 87.43, 87.49, 87.11};
  Float_t errRightIn[nRightIn] = {0.19, 0.19, 0.17, 0.20, 0.20, 0.17, 0.18};
  Float_t numRightIn[nRightIn] = {3003, 3012, 3017, 3023, 3045, 3049, 3053};

  for(Int_t i = 0; i < nLeftOut; i++){
    molSlugAvgs[0][0].push_back(avgLeftOut[i]);
    molPols[0][0].push_back(polLeftOut[i]);
    molErrs[0][0].push_back(errLeftOut[i]);
    molNums[0][0].push_back(numLeftOut[i]);
  }

  for(Int_t i = 0; i < nLeftIn; i++){
    molSlugAvgs[0][1].push_back(avgLeftIn[i]);
    molPols[0][1].push_back(polLeftIn[i]);
    molErrs[0][1].push_back(errLeftIn[i]);
    molNums[0][1].push_back(numLeftIn[i]);
  }

  for(Int_t i = 0; i < nRightOut; i++){
    molSlugAvgs[1][0].push_back(avgRightOut[i]);
    molPols[1][0].push_back(polRightOut[i]);
    molErrs[1][0].push_back(errRightOut[i]);
    molNums[1][0].push_back(numRightOut[i]);
  }

  for(Int_t i = 0; i < nRightIn; i++){
    molSlugAvgs[1][1].push_back(avgRightIn[i]);
    molPols[1][1].push_back(polRightIn[i]);
    molErrs[1][1].push_back(errRightIn[i]);
    molNums[1][1].push_back(numRightIn[i]);
  }
}


void fillMollerGraphs(){
  fillMollerVectors();
 
  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      for(Int_t k = 0; k < molSlugAvgs[i][j].size(); k++){
        gMoller[i][j]->SetPoint(nPtsMoller[i][j], molSlugAvgs[i][j][k], molPols[i][j][k]);
        gMoller[i][j]->SetPointError(nPtsMoller[i][j], 0.0, molErrs[i][j][k]);
        nPtsMoller[i][j]++;
      }
    }
  }
}


void mollerResiduals(){
  TCanvas *cRes = new TCanvas("cRes", "Residuals Canvas", 1500, 500);
  TGraphErrors *gResiduals = new TGraphErrors();
  gResiduals->SetTitle("Moller Residuals with Compton Fits");
  gResiduals->GetXaxis()->SetTitle("Moller Msmt #");
  gResiduals->GetYaxis()->SetTitle("Residual");
  gResiduals->SetMarkerStyle(8);
  gResiduals->SetMarkerSize(1.5);

  TF1 *fResFit = new TF1("fResFit", "pol0");
  TPaveText *pt = new TPaveText(0.75, 0.75, 0.90, 0.90, "blNDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(1);

  Int_t nPtsRes = 0;
  
  for(Int_t i = 0; i < nWiens; i++){
    for(Int_t j = 0; j < nIHWPs; j++){
      for(Int_t k = 0; k < molSlugAvgs[i][j].size(); k++){
        Int_t fitInd = -1;
        for(Int_t l = 0; l < nFuncs; l++){
          if(molSlugAvgs[i][j][k] >= fMins[l] && molSlugAvgs[i][j][k] <= fMaxs[l]){
            fitInd = l+j;
            break;
          }
        }
        Float_t res = molPols[i][j][k] - (par1[fitInd]*molSlugAvgs[i][j][k] + par0[fitInd]);
        printf("For Moller Msmt %i Avg Slug is %.1f, pol is %.2f, and fit parameters are %.3f and %.3f\n", molNums[i][j][k], molSlugAvgs[i][j][k], molPols[i][j][k], par0[fitInd], par1[fitInd]);
        gResiduals->SetPoint(nPtsRes, molNums[i][j][k], res);
        gResiduals->SetPointError(nPtsRes++, 0.0, molErrs[i][j][k]);
      }
    }
  }

  gResiduals->Fit(fResFit, "Q");

  pt->AddText("======Fit Parameters======");
  pt->AddText(Form("Par 0: %.4f +/- %.4f", fResFit->GetParameter(0), fResFit->GetParError(0)));
  pt->AddText(Form("Chi2 / NDF: %.3f / %i", fResFit->GetChisquare(), fResFit->GetNDF()));

  cRes->cd();
  gResiduals->Draw("ap");
  pt->Draw("same");
}


//=====================================================================================\\
//                                                                                     ||
//                                                                                     ||
//                             Calculate Avg Slugs for Esc                             ||
//                                                                                     ||
//                                                                                     ||
//=====================================================================================//


void getAvgSlugNums(){
  for(Int_t i = escLo; i <= escHi; i++){
    printf("Reading slugs from file escSlugs%i.list...\n", i);
    string slugNumStr;
    vector<Int_t> slugList;
    Int_t slugTotal = 0.0;
    ifstream infile(Form("%s/escargatoires/escSlugs%i.list", getenv("COMPMON_SNAILS"), i));
    if(infile.good()){
      while(getline(infile, slugNumStr)){
        Int_t slugNum = atoi(slugNumStr.c_str());
        if(slugNum == 0) 
          continue;
        slugList.push_back(slugNum);
        slugTotal += slugNum;
      }
      avgSlugs.push_back(slugTotal*1.0/((Float_t)slugList.size()));
      slugList.clear();
    }
  }
}


//=====================================================================================\\
//                                                                                     ||
//                                                                                     ||
//                         Translate snails to escargatoires                           ||
//                                                                                     ||
//                                                                                     ||
//=====================================================================================//


vector<Int_t> readEscSnails(Int_t esc){
  vector<Int_t> escSnls;
  printf("Reading snails from file esc%i.list...\n", esc);
  string snlNumStr;
  ifstream infile(Form("%s/escargatoires/esc%i.list", getenv("COMPMON_SNAILS"), esc));
  if(infile.good()){
    while(getline(infile, snlNumStr)){
      Int_t snlNum = atoi(snlNumStr.c_str());
      if(snlNum == 0) 
        continue;
      escSnls.push_back(snlNum);
    }
  }

  return escSnls;
}


vector<vector<Double_t>> getEscSnlData(Int_t esc, Int_t nMeas, vector<Int_t> escSnls, TString *names){
  TFile *f = TFile::Open(Form("%s/aggregates/crexGrandCompton.root", getenv("COMPMON_WEB")));
  TTree *snl = (TTree *)f->Get("snl");

  vector<vector<Double_t>> snlData;
  for(Int_t j = 0; j < escSnls.size(); j++){
    vector<Double_t> snlMeas;
    snlData.push_back(snlMeas);

    for(Int_t k = 0; k < nMeas; k++){
      TString hName = Form("esc%i_snl%i_%s", esc, escSnls[j], names[k].Data());
      snl->Draw(Form("%s>>%s", names[k].Data(), hName.Data()), Form("snailNum==%i", escSnls[j]), "goff");
      TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
      if(k < 2){
        snlData[j].push_back((Double_t)100.0*h->GetMean());
      }
      else{
        snlData[j].push_back((Double_t)1.0*h->GetMean());
      }
    }
  }

  return snlData;
}


void getEscPol0(){
  const Int_t nMeas = 5;
  vector<Double_t> polMeans;
  vector<Double_t> polErrs;
  TString names[nMeas] = {"Pol0.mean", "Pol0.meanErr", "sign", "ihwp", "VWienAngle"};

  for(Int_t i = escLo; i <= escHi; i++){
    vector<Int_t> escSnls = readEscSnails(i);
    vector<vector<Double_t>> snlData = getEscSnlData(i, nMeas, escSnls, names);
    
    Double_t wgtSum = 0.0;
    for(Int_t j = 0; j < snlData.size(); j++){
      wgtSum += 1.0/(snlData[j][1]*snlData[j][1]);
    }

    Double_t polAvg = 0.0;
    Double_t errAvg = 0.0;
    for(Int_t j = 0; j < escSnls.size(); j++){
      Double_t relWgt = (1.0/(snlData[j][1]*snlData[j][1]))*1.0/wgtSum;
      polAvg += relWgt*snlData[j][0];
      errAvg += 1.0/(snlData[j][1]*snlData[j][1]);
    }    

    polMeans.push_back(polAvg);
    polErrs.push_back(TMath::Sqrt(1.0/errAvg));
  }

  polData.push_back(polMeans);
  polData.push_back(polErrs);
}


void getIndData(){
  const Int_t nMeas = 5;
  TString names[nMeas] = {"Pol0.mean", "Pol0.meanErr", "sign", "ihwp", "VWienAngle"};

  for(Int_t i = escLo; i <= escHi; i++){
    vector<Int_t> escSnls = readEscSnails(i);
    vector<vector<Double_t>> snlData = getEscSnlData(i, nMeas, escSnls, names);

    Int_t escSign = (snlData[0][2] > 0.0) ? 1 : -1;
    Int_t escIHWP = (snlData[0][3] < 0.5) ? 0 : 1;
    Int_t escWien = (snlData[0][4] < 0.0) ? 0 : 1;

    signs.push_back(escSign);
    ihwps.push_back(escIHWP);
    wiens.push_back(escWien); 
  }
}


//=====================================================================================\\
//                                                                                     ||
//                                                                                     ||
//                                     Init Plots                                      ||
//                                                                                     ||
//                                                                                     ||
//=====================================================================================//


void initPlots(){
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

  initFunctions();
}


void crexEscwisePolWithMoller(){
  getAvgSlugNums();
  getEscPol0();
  getIndData();
  initPlots();

  for(Int_t i = 0; i < polData[0].size(); i++){
    if(signs[i] == 0) 
      continue;
    Int_t wien = wiens[i];
    Int_t ihwpInd = ihwps[i];
    gCompton[wien][ihwpInd]->SetPoint(nPtsCompton[wien][ihwpInd], avgSlugs[i], signs[i]*polData[0][i]);
    gCompton[wien][ihwpInd]->SetPointError(nPtsCompton[wien][ihwpInd], 0.0, polData[1][i]);
    nPtsCompton[wien][ihwpInd]++;
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
  mg->GetXaxis()->SetTitle("Average Slug Number");
  mg->GetYaxis()->SetTitle("Beam Polarization [pct]");
  // mg->GetXaxis()->SetLimits(xmin, xmax); mg->GetXaxis()->SetRangeUser(xmin, xmax);
  // mg->GetYaxis()->SetLimits(84, 90); mg->GetYaxis()->SetRangeUser(84, 90);
  mg->Draw("a");
  leg->Draw("same");

  for(Int_t i = 0; i < nFuncs; i++){
    funcs[i]->Draw("same");
  }

  //mg->Draw("a && same");

  mollerResiduals();
}
