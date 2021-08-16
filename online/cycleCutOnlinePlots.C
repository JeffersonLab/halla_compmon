#include "utils.h"
#include "../grand/grandOnline/makePlots.h"

using namespace std;

Int_t msmtNum = 0;

TString getCycleCut(Int_t prexOrCrex, Int_t runNum, TString cutVar){
  TString cut = "";
  TString expt = experimentCode(prexOrCrex);
  TFile *grand = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree* cycRaw = (TTree *)grand->Get("cyc");
  TTree* cyc = cycRaw->CopyTree(Form("runNum==%i && CycleCut==0", runNum));
  Int_t cycCut, start, stop;
  cyc->SetBranchAddress("CycleCut", &cycCut);
  cyc->SetBranchAddress("firstOffStartMPS", &start);
  cyc->SetBranchAddress("lastOffEndMPS", &stop);

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(i > 0){
      cut += " || ";
    }
    cut += Form("(%s>=%i && %s<=%i)", cutVar.Data(), start, cutVar.Data(), stop);
  }

  return ("(" + cut + ")");
}

vector<vector<Int_t>> getCycleListAlt(Int_t prexOrCrex, Int_t runNum){
  vector<vector<Int_t>> runCycles;
  TString expt = experimentCode(prexOrCrex);
  TFile *grand = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree* cycRaw = (TTree *)grand->Get("cyc");
  TTree* cyc = cycRaw->CopyTree(Form("runNum==%i && CycleCut==0", runNum));
  Int_t start, stop;
  cyc->SetBranchAddress("firstOffStartMPS", &start);
  cyc->SetBranchAddress("lastOffEndMPS", &stop);
  
  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    vector<Int_t> cycLimits;
    cycLimits.push_back(start);
    cycLimits.push_back(stop);
    runCycles.push_back(cycLimits);
  }

  return runCycles;
}

Float_t avg(vector<Float_t> data){
  Float_t tot = 0.0;
  for(Float_t datum : data){tot += datum;}
  return tot/((int)data.size());
}

Bool_t mpsInCycleAlt(vector<vector<int>> runCycles, Int_t mpsNum){
  for(Int_t i = 0; i < runCycles.size(); i++){
    if(mpsNum >= runCycles[i][0] && mpsNum <= runCycles[i][1]){
      return kTRUE;
    }
  }
  return kFALSE;
}

void makeHisto(TTree *somethingwise, TH1F *h, TString name, TString title, TString msmt, TString cuts){
  h = (TH1F *)gDirectory->Get(name.Data()); if(h){delete h;}
  somethingwise->Draw(Form("%s>>%s", msmt.Data(), name.Data()), cuts.Data(), "goff");
  h = (TH1F *)gDirectory->Get(name.Data()); h->SetTitle(title.Data());
  h->GetXaxis()->SetTitle(msmt.Data());
}

void drawSubCanvas(TCanvas *c, Int_t subCanvas, TH1F *h){
  c->cd(subCanvas);
  h->Draw();
  vector<TString> s1 = hist_stats(h);
  TPaveText *pt1 = new TPaveText(0.75, 0.7, 0.95, 0.9, "blNDC"); 
  pt1->SetBorderSize(1); pt1->SetFillColor(0);
  for(int i = 0; i < s1.size(); i++){
    pt1->AddText(s1[i].Data())->SetTextColor(kBlue);
  } 
  pt1->Draw("same");
}

void drawSubCanvasTwo(TCanvas *c, Int_t subCanvas, TH1F *hON, TH1F *hOFF, TString name, TString title){
  c->cd(subCanvas);
  THStack *hs = new THStack(name.Data(), title.Data());
  hON->SetLineColor(kGreen + 2); hOFF->SetLineColor(kRed);
  hs->Add(hON); hs->Add(hOFF);
  hs->Draw("nostack");

  vector<TString> s1 = hist_stats(hON);
  TPaveText *pt1 = new TPaveText(0.75, 0.7, 0.95, 0.9, "blNDC");
  pt1->AddText("---- Laser ON ----")->SetTextColor(kGreen + 2);
  pt1->SetBorderSize(1); pt1->SetFillColor(0);
  for(int i = 0; i < s1.size(); i++){
    pt1->AddText(s1[i].Data())->SetTextColor(kGreen + 2);
  } 
  pt1->Draw("same");

  vector<TString> s2 = hist_stats(hOFF);
  TPaveText *pt2 = new TPaveText(0.75, 0.5, 0.95, 0.7, "blNDC"); 
  pt2->AddText("---- Laser OFF ----")->SetTextColor(kRed);
  pt2->SetBorderSize(1); pt2->SetFillColor(0);
  for(int i = 0; i < s2.size(); i++){
    pt2->AddText(s2[i].Data())->SetTextColor(kRed);
  } 
  pt2->Draw("same");
}

void makeGraph(TCanvas *c, Int_t canSub, TString name, TString title, vector<Float_t> xData, vector<Float_t> yData, Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, 
          TString msmtX, TString msmtY){
  Double_t xLoFac = 0.98; Double_t xHiFac = 1.02; Double_t yLoFac = 0.95; Double_t yHiFac = 1.05;
  if(xmin < 0){xLoFac = 1.02;} if(xmax < 0){xHiFac = 0.98;}
  if(ymin < 0){yLoFac = 1.05;} if(ymax < 0){yHiFac = 0.95;}
  TGraph *g = (TGraph *)gDirectory->Get(name.Data());
  if(g){delete g;}
  g = new TGraph(xData.size(), xData.data(), yData.data());
  g->SetName(name.Data()); g->SetTitle(title.Data());
  g->GetXaxis()->SetTitle(msmtX.Data()); g->GetYaxis()->SetTitle(msmtY.Data());
  g->GetXaxis()->SetLimits(xmin*xLoFac, xmax*xHiFac); g->GetYaxis()->SetRangeUser(ymin*yLoFac, ymax*yHiFac);

  c->cd(canSub);
  g->Draw("ap");
}

void makeMultiGraph(TCanvas *c, Int_t canSub, TMultiGraph* mg, Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, TString msmtX, TString msmtY){
  Double_t xLoFac = 0.98; Double_t xHiFac = 1.02; Double_t yLoFac = 0.95; Double_t yHiFac = 1.05;
  if(xmin < 0){xLoFac = 1.02;} if(xmax < 0){xHiFac = 0.98;}
  if(ymin < 0){yLoFac = 1.05;} if(ymax < 0){yHiFac = 0.95;}
  mg->GetXaxis()->SetTitle(msmtX.Data()); mg->GetYaxis()->SetTitle(msmtY.Data());
  mg->GetXaxis()->SetLimits(xmin*xLoFac, xmax*xHiFac); mg->GetYaxis()->SetRangeUser(ymin*yLoFac, ymax*yHiFac);

  c->cd(canSub);
  mg->Draw("ap");
}


void breakdownPlots(Int_t prexOrCrex, TChain* somethingwise, int run_num, TString msmt, TString tree_id, TString msmt_id, int max_event){
  printf("Plotting %s...\n", msmt_id.Data());
  
  TString states[2] = {"ON", "OFF"};
  TString cycCuts = getCycleCut(prexOrCrex, run_num, "mpsCoda");
  TString laserON("(laserState==0 || laserState==1)"); TString laserOFF("(laserState==2 || laserState==3)");
  TString beamON = Form("beamState==1"); TString beamOFF("beamState==0");
  TString laserCuts[2] = {laserON, laserOFF}; TString beamCuts[2] = {beamON, beamOFF};
  TString extraCuts(Form("mpsCoda<%i && %s", max_event, cycCuts.Data()));

  TString name = Form("%s_%s", tree_id.Data(), msmt_id.Data());
  TString title = Form("Run %i %s, %s: All", run_num, tree_id.Data(), msmt_id.Data());
  TH1F *h; makeHisto(somethingwise, h, name, title, msmt, extraCuts.Data());
  for(int i = 0; i < 2; i++){
    TString cuts1 = Form("%s && %s", beamCuts[i].Data(), extraCuts.Data());
    TString name1 = Form("%s_beam%i_%s", tree_id.Data(), (Int_t)!i, msmt_id.Data());
    TString title1 = Form("Run %i %s, %s: Beam %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data());
    TH1F *h1; makeHisto(somethingwise, h1, name1, title1, msmt, cuts1);
    for(int j = 0; j < 2; j++){
      TString cuts2 = Form("%s && %s && %s", beamCuts[i].Data(), laserCuts[j].Data(), extraCuts.Data());
      TString name2 = Form("%s_beam%i_las%i_%s", tree_id.Data(), (Int_t)!i, (Int_t)!j, msmt_id.Data());
      TString title2 = Form("Run %i, %s, %s: Beam %s, Laser %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data(), states[j].Data());
      TH1F *h2; makeHisto(somethingwise, h2, name2, title2, msmt, cuts2);
      TString cuts2A = Form("%s && %s", laserCuts[j].Data(), extraCuts.Data());
      TString name2A = Form("%s_las%i_%s", tree_id.Data(), (Int_t)!j, msmt_id.Data());
      TString title2A = Form("Run %i, %s, %s: Laser %s", run_num, tree_id.Data(), msmt_id.Data(), states[j].Data());
      TH1F *h2A; makeHisto(somethingwise, h2A, name2A, title2A, msmt, cuts2A);
    }
  }
  
  TCanvas *c = new TCanvas(Form("c%s", msmt.Data()), Form("%s Canvas", msmt.Data()), 1200, 800);
  c->Divide(2, 2);
  drawSubCanvas(c, 1, (TH1F *)gDirectory->Get(Form("%s_%s", tree_id.Data(), msmt_id.Data())));
  drawSubCanvas(c, 2, (TH1F *)gDirectory->Get(Form("%s_beam0_%s", tree_id.Data(), msmt_id.Data())));
  drawSubCanvas(c, 3, (TH1F *)gDirectory->Get(Form("%s_beam1_las1_%s", tree_id.Data(), msmt_id.Data())));
  drawSubCanvas(c, 4, (TH1F *)gDirectory->Get(Form("%s_beam1_las0_%s", tree_id.Data(), msmt_id.Data())));

  c->SaveAs(Form("%s/runs/Run%i/mps_%04i.png", getenv("COMPMON_WEB"), run_num, msmtNum++), "png");
}

void mpsGraphs(Int_t prexOrCrex, TChain* somethingwise, int run_num, TString tree_id, int max_event){
  cout<<"Plotting mps graphs..."<<endl;

  //Initialize strings for plots
  const Int_t nPlots = 4;
  const Int_t nVars = 7;
  TString names[nPlots] = {"beam ON, laser ON", "beam ON, laser OFF", "beam OFF, laser ON", "beam OFF, laser OFF"};
  TString msmt_ids[nVars] = {"acc0_time", "bcm_time", "cavPower_time", "bpmAx_time", "bpmAy_time", "bpmBx_time", "bpmBy_time"};
  TString ytitles[nVars] = {"Acc0/NAcc0", "BCM", "Cav Power", "BPM Ax", "BPM Ay", "BPM Bx", "BPM By"};

  //Initialize Acc0, BPM, BCM plots y_vals[msmt_id][config]
  vector<vector<int>> runCycles = getCycleListAlt(prexOrCrex, run_num);
  vector<vector<Float_t>> mpsCoda_vals; vector<vector<vector<Float_t>>> y_vals;
  vector<TCanvas *> cans;
  for(int i = 0; i < nVars; i++){
    vector<vector<Float_t>> real_vals;
    for(int j = 0; j < nPlots; j++){
      if(i == 0){vector<Float_t> mpsCoda_val; mpsCoda_vals.push_back(mpsCoda_val);}
      vector<Float_t> real_val; real_vals.push_back(real_val);
    }
    y_vals.push_back(real_vals);
    TCanvas *c = new TCanvas(Form("cMPS%s", msmt_ids[i].Data()), Form("%s vs Time Canvas", ytitles[i].Data()), 1200, 800);
    c->Divide(2, 2);
    cans.push_back(c);
  }

  Int_t mpsCoda, laserState, beamState;
  Double_t acc0, nacc0;
  Float_t bpmAx, bpmAy, bpmBx, bpmBy, bcm, cavPower;

  TString cycCuts = getCycleCut(prexOrCrex, run_num, "mpsCoda");
  TString treeCuts = Form("mpsCoda<%i && beamState<2 && laserState<4 && %s", max_event, cycCuts.Data());

  somethingwise->SetBranchAddress("mpsCoda", &mpsCoda);
  somethingwise->SetBranchAddress("laserState", &laserState); somethingwise->SetBranchAddress("beamState", &beamState);
  somethingwise->SetBranchAddress("Acc0", &acc0); somethingwise->SetBranchAddress("NAcc0", &nacc0);
  somethingwise->SetBranchAddress("bcm", &bcm); somethingwise->SetBranchAddress("cavPowerCalibrated", &cavPower);
  somethingwise->SetBranchAddress("bpmAx", &bpmAx); somethingwise->SetBranchAddress("bpmAy", &bpmAy);
  somethingwise->SetBranchAddress("bpmBx", &bpmBx); somethingwise->SetBranchAddress("bpmBy", &bpmBy);

  Float_t xmin = 1e16; Float_t xmax = -1e16;
  Float_t ymins[nVars] = { 1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16}; 
  Float_t ymaxs[nVars] = {-1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16};

  for(int i = 0; i < somethingwise->GetEntries(); i++){
    somethingwise->GetEntry(i);
    if(laserState > 3 || beamState > 1 || mpsCoda > max_event || (!mpsInCycleAlt(runCycles, mpsCoda))){continue;}
    Int_t laserInd = (laserState == 0 || laserState == 1) ? 0 : 1;
    Int_t beamInd = beamState == 1 ? 0 : 2;
    mpsCoda_vals[beamInd + laserInd].push_back(i);
    if(mpsCoda < xmin){xmin = i;}
    if(mpsCoda > xmax){xmax = i;}
    Double_t all_vals[nVars] = {acc0/nacc0, bcm, cavPower, bpmAx, bpmAy, bpmBx, bpmBy};
    for(int j = 0; j < y_vals.size(); j++){
      y_vals[j][beamInd + laserInd].push_back(all_vals[j]);
      if(all_vals[j] < ymins[j]){ymins[j] = all_vals[j];}
      if(all_vals[j] > ymaxs[j]){ymaxs[j] = all_vals[j];}
    }
  }

  for(int i = 0; i < nVars; i++){
    for(int j = 0; j < nPlots; j++){
      Int_t beam = (Int_t)(j < 2);
      Int_t laser = (Int_t)(j % 2 == 0);
	    TString name = Form("%s_beam%i_las%i_%s", tree_id.Data(), beam, laser, msmt_ids[i].Data()); 
      TString title = Form("Run %i %s, %s", run_num, tree_id.Data(), names[j].Data());
      makeGraph(cans[i], j+1, name, title, mpsCoda_vals[j], y_vals[i][j], xmin, xmax, ymins[i], ymaxs[i], "mpsCoda", ytitles[i]);
    }
    cans[i]->SaveAs(Form("%s/runs/Run%i/mps_%04i.png", getenv("COMPMON_WEB"), run_num, msmtNum++), "png");
  }
}

void quartetPlots(Int_t prexOrCrex, TChain* somethingwise, int run_num, int max_event){
  cout<<"Plotting quartets..."<<endl;
  
  TString states[2] = {"ON", "OFF"};
  TString cycCuts = getCycleCut(prexOrCrex, run_num, "firstMPSnumber");
  TString laserON("(laserState==0 || laserState==1)"); TString laserOFF("(laserState==2 || laserState==3)");
  TString beamON = Form("beamState==1"); TString beamOFF("beamState==0");
  TString laserCuts[2] = {laserON, laserOFF}; TString beamCuts[2] = {beamON, beamOFF};
  TString extraCuts(Form("firstMPSnumber < %i && %s", max_event, cycCuts.Data()));

  TString posH0("PosHelAcc0/PosHelNSamples0");
  TString negH0("NegHelAcc0/NegHelNSamples0");
  TString diff0 = Form("%s - %s", posH0.Data(), negH0.Data());
  TString summ0 = Form("%s + %s", posH0.Data(), negH0.Data());
  TString posH4("PosHelAcc4");
  TString negH4("NegHelAcc4");
  TString diff4 = Form("%s - %s", posH4.Data(), negH4.Data());
  TString summ4 = Form("%s + %s", posH4.Data(), negH4.Data());
  TString asymUSbg1("(PosHelScalerRun4  - NegHelScalerRun4) /(PosHelScalerRun4  + NegHelScalerRun4)");
  TString asymUSbg2("(PosHelScalerRun13 - NegHelScalerRun13)/(PosHelScalerRun13 + NegHelScalerRun13)");
  TString asymDSbg1("(PosHelScalerRun0  - NegHelScalerRun0) /(PosHelScalerRun0 + NegHelScalerRun0)");
  TString asymDSbg2("(PosHelScalerRun1  - NegHelScalerRun1) /(PosHelScalerRun1 + NegHelScalerRun1)");
  TString asymHoriz("(PosHelScalerRun10 - NegHelScalerRun10)/(PosHelScalerRun10 + NegHelScalerRun10)");
  TString asymVert ("(PosHelScalerRun11 - NegHelScalerRun11)/(PosHelScalerRun11 + NegHelScalerRun11)");
  
  const Int_t nVars = 18;
  TString meas[nVars] = {posH0, negH0, diff0, summ0, posH4, negH4, diff4, summ4, asymUSbg1, asymUSbg2, asymDSbg1, asymDSbg2, asymHoriz, asymVert, "", "", "", ""};
  TString meas_id[nVars] = {"PosHelAcc0/NSamples0", "NegHelAcc0/NSamples0", "Diffs0", "Sums0", "PosHelAcc4", "NegHelAcc4", "Diffs4", "Sums4", "Asyms USbg1", "Asyms USbg2", "Asyms DSbg1", "Asyms DSbg2", 
                            "Asyms Horiz Finger", "Asyms Vert Finger", "Asyms0", "Asyms4", "Asyms0", "Asyms4"};
  TString meas_small[nVars] = {"posH0", "negH0", "diff0", "summ0", "posH4", "negH4", "diff4", "summ4", "USbg1", "USbg2", "DSbg1", "DSbg2", 
                               "horizFing", "vertFing", "asym0subOn",  "asym4subOn", "asym0subOff", "asym4subOff"};
  TString tree_id("quartetwise");

  const Int_t nCans = 3;
  const Int_t nGraphs = 6;
  TCanvas *cans[nCans];
  TString titles[nCans][nGraphs];
  for(Int_t i = 0; i < nCans; i++){
    TCanvas *c = new TCanvas(Form("cQuartet%i", i+1), Form("Quartet Canvas %i", i+1), 1200, 800);
    c->Divide(2, 3);
    cans[i] = c;
    for(Int_t j = 0; j < nGraphs; j++){
      Int_t base = 4*((Int_t)(i==1)) + 8*((Int_t)(i==2));
      Int_t ind = ((Int_t)((i<2 && j<4) || i==2))*(base + j) + ((Int_t)(i<2 && j==4))*(14 + i);
      TString title = Form("Run %i %s, %s", run_num, tree_id.Data(), meas_id[ind].Data());
      titles[i][j] = title;
    }
  }

  Float_t summOn0 = 0.0; Float_t summOff0 = 0.0;
  Float_t summOn4 = 0.0; Float_t summOff4 = 0.0;
  for(int i = 0; i < nVars; i++){
    if(i==(nVars - 4)){
      meas[i]   = Form("1000*(%s)/(%s - %f)", diff0.Data(), summ0.Data(), summOff0);
      meas[i+1] = Form("1000*(%s)/(%s - %f)", diff4.Data(), summ4.Data(), summOff4);
      meas[i+2] = Form("1000*(%s)/(%f - %f)", diff0.Data(), summOn0, summOff0);
      meas[i+3] = Form("1000*(%s)/(%f - %f)", diff4.Data(), summOn4, summOff4);
    }
    TString hname = Form("%s_%s", tree_id.Data(), meas_small[i].Data());
    TString title = Form("Run %i %s, %s: All", run_num, tree_id.Data(), meas_id[i].Data());
    TH1F *h; makeHisto(somethingwise, h, hname, title, meas[i], extraCuts.Data());
    for(int j = 0; j < 2; j++){
      TString cuts1 = Form("%s && %s", beamCuts[j].Data(), extraCuts.Data());
      TString hname1 = Form("%s_beam%i_%s", tree_id.Data(), (Int_t)!j, meas_small[i].Data());
      TString title1 = Form("Run %i %s, %s: Beam %s", run_num, tree_id.Data(), meas_id[i].Data(), states[j].Data());
      TH1F *h1; makeHisto(somethingwise, h1, hname1, title1, meas[i], cuts1);
      for(int k = 0; k < 2; k++){
        TString cuts2 = Form("%s && %s", laserCuts[k].Data(), extraCuts.Data());
        TString hname2 = Form("%s_las%i_%s", tree_id.Data(), (Int_t)!k, meas_small[i].Data());
        TString title2 = Form("Run %i %s, %s: Laser %s", run_num, tree_id.Data(), meas_id[i].Data(), states[k].Data());
        TH1F *h2; makeHisto(somethingwise, h2, hname2, title2, meas[i], cuts2);
        TString cuts3 = Form("%s && %s && %s", beamCuts[j].Data(), laserCuts[k].Data(), extraCuts.Data());
        TString hname3 = Form("%s_beam%i_las%i_%s", tree_id.Data(), (Int_t)!j, (Int_t)!k, meas_small[i].Data());
        TString title3 = Form("Run %i %s, %s: Beam %s, Laser %s", run_num, tree_id.Data(), meas_id[i].Data(), states[j].Data(), states[k].Data());
        TH1F *h3; makeHisto(somethingwise, h3, hname3, title3, meas[i], cuts3); h3 = (TH1F *)gDirectory->Get(hname3.Data());
        if(i==3 && j==0 && k==0){summOn0 = h3->GetMean();}
        if(i==7 && j==0 && k==0){summOn4 = h3->GetMean();}
        if(i==3 && j==0 && k==1){summOff0 = h3->GetMean();}
        if(i==7 && j==0 && k==1){summOff4 = h3->GetMean();}
      }
    }
    if(i < 14){
      Int_t canInd = ((Int_t)(i>=4 && i<=7)) + 2*((Int_t)(i>=8 && i<=13));
      Int_t gInd = ((Int_t)(i<=3))*i + ((Int_t)(i>=4 && i<=7))*(i - 4) + ((Int_t)(i>=8 && i<=13))*(i - 8);
      
      TH1F *hON = (TH1F *)gDirectory->Get(Form("%s_beam1_las1_%s", tree_id.Data(), meas_small[i].Data()));
      TH1F *hOFF = (TH1F *)gDirectory->Get(Form("%s_beam1_las0_%s", tree_id.Data(), meas_small[i].Data()));

      drawSubCanvasTwo(cans[canInd], gInd+1, hON, hOFF, Form("%s_stack_%s", tree_id.Data(), meas_small[i].Data()), titles[canInd][gInd]);
    }
    else if(i >= 16){
      TH1F *hON =  (TH1F *)gDirectory->Get(Form("%s_beam1_las1_%s", tree_id.Data(), meas_small[i].Data()));
      TH1F *hOFF = (TH1F *)gDirectory->Get(Form("%s_beam1_las0_%s", tree_id.Data(), meas_small[i].Data()));
      
      drawSubCanvasTwo(cans[i - 16], 5, hON, hOFF, Form("%s_stack_%s", tree_id.Data(), meas_small[i].Data()), titles[i - 16][4]);
    }
  }

  for(Int_t i = 0; i < nCans; i++){
    cans[i]->SaveAs(Form("%s/runs/Run%i/qrt_%04i.png", getenv("COMPMON_WEB"), run_num, msmtNum++), "png");
  }
  
}

void quartetGraphs(Int_t prexOrCrex, TChain* somethingwise, int run_num, TString tree_id, int max_event){
  cout<<"Plotting quartets vs time..."<<endl;

  TString states[2] = {"OFF", "ON"};
  TString cycCuts = getCycleCut(prexOrCrex, run_num, "firstMPSnumber");
  const Int_t nMsmts = 18;
  TString msmts[nMsmts] = {"posH0", "negH0", "diff0", "summ0", "bcm", "posH4", "negH4", "diff4", "summ4",
                           "USbg1", "USbg2", "DSbg1", "DSbg2", "vertFing", "horizFing", "cavPow", "asym0sub", "asym4sub"};

  vector<vector<int>> runCycles = getCycleListAlt(prexOrCrex, run_num);
  vector<vector<vector<Float_t>>> x_vals; vector<vector<vector<Float_t>>> y_vals;
  
  for(int i = 0; i < nMsmts; i++){
    vector<vector<Float_t>> x; vector<vector<Float_t>> y;
    vector<Float_t> x_l0_b0; vector<Float_t> x_l0_b1; vector<Float_t> x_l1_b0; vector<Float_t> x_l1_b1;
    vector<Float_t> y_l0_b0; vector<Float_t> y_l0_b1; vector<Float_t> y_l1_b0; vector<Float_t> y_l1_b1;
    x.push_back(x_l0_b0); x.push_back(x_l0_b1); x.push_back(x_l1_b0); x.push_back(x_l1_b1);
    y.push_back(y_l0_b0); y.push_back(y_l0_b1); y.push_back(y_l1_b0); y.push_back(y_l1_b1);
    x_vals.push_back(x); y_vals.push_back(y);
  }

  const Int_t nCans = 3;
  const Int_t nGraphs = 6;
  TCanvas *cans[nCans];
  TMultiGraph *mgs[nCans][nGraphs];
  for(Int_t i = 0; i < nCans; i++){
    TCanvas *c = new TCanvas(Form("cQuartetGraphs%i", i+1), Form("Quartet Graph Canvas %i", i+1), 1200, 800);
    c->Divide(2, 3);
    cans[i] = c;
    for(Int_t j = 0; j < nGraphs; j++){
      Int_t base = 5*((Int_t)(i==1)) + 9*((Int_t)(i==2));
      Int_t ind = ((Int_t)((i<2 && j<4) || i==2))*(base + j) + ((Int_t)(i<2 && j==4))*(16 + i) + 5*((Int_t)(i<2 && j==5));
      TString name = Form("%s_%s_time", tree_id.Data(), msmts[ind].Data());
      TString title = Form("Run %i %s, %s vs time", run_num, tree_id.Data(), msmts[ind].Data());
      mgs[i][j] = new TMultiGraph(name.Data(), title.Data());
    }
  }

  Int_t laserState, beamState, firstMPSnumber;
  Double_t posAcc, negAcc, posSamp, negSamp, posAcc4, negAcc4, posSamp4, negSamp4;
  Int_t posUSbg1, negUSbg1, posUSbg2, negUSbg2, posDSbg1, negDSbg1, posDSbg2, negDSbg2,
        posHoriz, negHoriz, posVert, negVert;
  Float_t posBCM, negBCM, posCavPow, negCavPow;

  somethingwise->SetBranchAddress("PosHelBCM", &posBCM); somethingwise->SetBranchAddress("NegHelBCM", &negBCM);
  somethingwise->SetBranchAddress("PosHelCavPowerCalibrated", &posCavPow); somethingwise->SetBranchAddress("NegHelCavPowerCalibrated", &negCavPow);
  somethingwise->SetBranchAddress("laserState", &laserState); somethingwise->SetBranchAddress("beamState", &beamState);
  somethingwise->SetBranchAddress("PosHelAcc0", &posAcc); somethingwise->SetBranchAddress("PosHelNSamples0", &posSamp);
  somethingwise->SetBranchAddress("NegHelAcc0", &negAcc); somethingwise->SetBranchAddress("NegHelNSamples0", &negSamp);
  somethingwise->SetBranchAddress("PosHelAcc4", &posAcc4); somethingwise->SetBranchAddress("PosHelNSamples4", &posSamp4);
  somethingwise->SetBranchAddress("NegHelAcc4", &negAcc4); somethingwise->SetBranchAddress("NegHelNSamples4", &negSamp4);
  somethingwise->SetBranchAddress("PosHelScalerRun4",  &posUSbg1); somethingwise->SetBranchAddress("NegHelScalerRun4",  &negUSbg1);
  somethingwise->SetBranchAddress("PosHelScalerRun13", &posUSbg2); somethingwise->SetBranchAddress("NegHelScalerRun13", &negUSbg2);
  somethingwise->SetBranchAddress("PosHelScalerRun0",  &posDSbg1); somethingwise->SetBranchAddress("NegHelScalerRun0",  &negDSbg1);
  somethingwise->SetBranchAddress("PosHelScalerRun1",  &posDSbg2); somethingwise->SetBranchAddress("NegHelScalerRun1",  &negDSbg2);
  somethingwise->SetBranchAddress("PosHelScalerRun10", &posHoriz); somethingwise->SetBranchAddress("NegHelScalerRun10", &negHoriz);
  somethingwise->SetBranchAddress("PosHelScalerRun11", &posVert);  somethingwise->SetBranchAddress("NegHelScalerRun11", &negVert);
  somethingwise->SetBranchAddress("firstMPSnumber", &firstMPSnumber);

  Double_t ymins[nMsmts] = { 1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  
                             1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16};
  Double_t ymaxs[nMsmts] = {-1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, 
                            -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16};
  for(int i = 0; i < somethingwise->GetEntries(); i++){
    somethingwise->GetEntry(i);
    if(laserState > 3 || beamState > 1 || firstMPSnumber > max_event || (!mpsInCycleAlt(runCycles, firstMPSnumber))){continue;}
    Int_t laserInd = (Int_t)(laserState==0 || laserState==1);
    Int_t beamInd  = (Int_t)(beamState==1);
    Double_t msmt[nMsmts - 2] = {posAcc/posSamp, negAcc/negSamp, posAcc/posSamp - negAcc/negSamp, posAcc/posSamp + negAcc/negSamp, 
                                 (posBCM + negBCM)/2.0, posAcc4, negAcc4, posAcc4 - negAcc4, posAcc4 + negAcc4,
                                 (posUSbg1 + negUSbg1)*30.0, (posUSbg2 + negUSbg2)*30.0, (posDSbg1 + negDSbg1)*30.0, (posDSbg2 + negDSbg2)*30.0,
                                 (posHoriz + negHoriz)*30.0, (posVert + negVert)*30.0, (posCavPow + negCavPow)/2.0};
    for(int j = 0; j < nMsmts - 2; j++){
      if(msmt[j] < ymins[j]){ymins[j] = msmt[j];}
      if(msmt[j] > ymaxs[j]){ymaxs[j] = msmt[j];}
      x_vals[j][2*laserInd + beamInd].push_back(i);
      y_vals[j][2*laserInd + beamInd].push_back(msmt[j]);
    }
  }

  Float_t summOn0  = avg(y_vals[3][3]);
  Float_t summOff0 = avg(y_vals[3][1]);
  Float_t summOn4  = avg(y_vals[8][3]);
  Float_t summOff4 = avg(y_vals[8][1]);
  for(int i = 0; i < y_vals[4][3].size(); i++){
    Float_t asym0 = 1000*y_vals[2][3][i]/(summOn0 - summOff0);
    Float_t asym4 = 1000*y_vals[7][3][i]/(summOn4 - summOff4);
    if(asym0 < ymaxs[16]){ymins[16] = asym0;} if(asym0 > ymaxs[16]){ymaxs[16] = asym0;}
    if(asym4 < ymins[17]){ymins[17] = asym4;} if(asym4 > ymaxs[17]){ymaxs[17] = asym4;}
    x_vals[16][3].push_back(x_vals[2][3][i]);
    x_vals[17][3].push_back(x_vals[7][3][i]);
    y_vals[16][3].push_back(asym0);
    y_vals[17][3].push_back(asym4);
  }
  for(int i = 0; i < y_vals[4][1].size(); i++){
    Float_t asym0 = 1000*y_vals[2][1][i]/(summOn0 - summOff0);
    Float_t asym4 = 1000*y_vals[7][1][i]/(summOn4 - summOff4);
    if(asym0 < ymins[16]){ymins[16] = asym0;} if(asym0 > ymaxs[16]){ymaxs[16] = asym0;}
    if(asym4 < ymins[17]){ymins[17] = asym4;} if(asym4 > ymaxs[17]){ymaxs[17] = asym4;}
    x_vals[16][1].push_back(x_vals[2][1][i]);
    x_vals[17][1].push_back(x_vals[7][1][1]);
    y_vals[16][1].push_back(asym0);
    y_vals[17][1].push_back(asym4);
  }

  for(int i = 0; i < x_vals.size(); i++){
    TGraph *g1 = new TGraph(x_vals[i][0].size(), x_vals[i][0].data(), y_vals[i][0].data());
    TGraph *g2 = new TGraph(x_vals[i][1].size(), x_vals[i][1].data(), y_vals[i][1].data());
    g2->SetMarkerColor(kRed);
    TGraph *g3 = new TGraph(x_vals[i][2].size(), x_vals[i][2].data(), y_vals[i][2].data());
    TGraph *g4 = new TGraph(x_vals[i][3].size(), x_vals[i][3].data(), y_vals[i][3].data());
    g4->SetMarkerColor(kGreen + 2);

    if(i != 4 && i != 15){
      Int_t canInd = ((Int_t)((i>=5 && i<=8) || i==17)) + 2*((Int_t)(i>=9 && i<=14));
      Int_t gInd = ((Int_t)(i<=3))*i + ((Int_t)(i>=5 && i<=8))*(i - 5) + ((Int_t)(i>=9 && i<=14))*(i - 9) + ((Int_t)(i>=16))*4;

      mgs[canInd][gInd]->Add(g2);
      mgs[canInd][gInd]->Add(g4);
      mgs[canInd][gInd]->Add(g1);
      mgs[canInd][gInd]->Add(g3);

      makeMultiGraph(cans[canInd], gInd+1, mgs[canInd][gInd], 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);
    }
    else if(i == 4){
      mgs[0][5]->Add(g1); mgs[0][5]->Add(g2); mgs[0][5]->Add(g3); mgs[0][5]->Add(g4);
      mgs[1][5]->Add(g1); mgs[1][5]->Add(g2); mgs[1][5]->Add(g3); mgs[1][5]->Add(g4);

      makeMultiGraph(cans[0], 6, mgs[0][5], 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);
      makeMultiGraph(cans[1], 6, mgs[1][5], 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);
    }
  
    /**
    TString name1 = Form("%s_all_%s_time", tree_id.Data(), msmts[i].Data());
    TString title1 = Form("Run %i %s, %s vs time: all", run_num, tree_id.Data(), msmts[i].Data());
    TMultiGraph *mg1 = new TMultiGraph(name1.Data(), title1.Data());
    mg1->Add(g1); mg1->Add(g2); mg1->Add(g3); mg1->Add(g4);
    makeMultiGraph(cans[i], 1, mg1, 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);

    TString name2 = Form("%s_beam0_%s_time", tree_id.Data(), msmts[i].Data());
    TString title2 = Form("Run %i %s, %s vs time: beam OFF", run_num, tree_id.Data(), msmts[i].Data());
    TMultiGraph *mg2 = new TMultiGraph(name2.Data(), title2.Data());
    mg2->Add(g1); mg2->Add(g3);
    makeMultiGraph(cans[i], 2, mg2, 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);

    TString name3 = Form("%s_beam1_las1_%s_time", tree_id.Data(), msmts[i].Data());
    TString title3 = Form("Run %i %s, %s vs time: beam ON, laser ON", run_num, tree_id.Data(), msmts[i].Data());
    TMultiGraph *mg3 = new TMultiGraph(name3.Data(), title3.Data());
    mg3->Add(g4);
    makeMultiGraph(cans[i], 3, mg3, 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);

    TString name4 = Form("%s_beam1_las0_%s_time", tree_id.Data(), msmts[i].Data());
    TString title4 = Form("Run %i %s, %s vs time: beam ON, laser OFF", run_num, tree_id.Data(), msmts[i].Data());
    TMultiGraph *mg4 = new TMultiGraph(name4.Data(), title4.Data());
    mg4->Add(g3);
    makeMultiGraph(cans[i], 4, mg4, 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);
    **/
  }

  for(Int_t i = 0; i < nCans; i++){
    cans[i]->SaveAs(Form("%s/runs/Run%i/qrt_%04i.png", getenv("COMPMON_WEB"), run_num, msmtNum++), "png");
  }
}


void cycleCutOnlinePlots(Int_t runNum, Int_t maxEvt=1e9){
  Int_t prexOrCrex = (runNum < 4800) ? 1 : 2;

  Int_t mpswise = 0; Int_t quartetwise = 1;
  Int_t pulserwise = 2; Int_t triggerwise = 3;
  Int_t epicswise = 4; Int_t runwise = 5;
  Int_t snapshots = 6;
  vector<TChain *> runChains = loadChain(runNum);

  gStyle->SetOptStat(0);
  const Int_t nVars = 7;
  TString vars[nVars] = {"Acc0/NAcc0", "bcm", "cavPowerCalibrated", "bpmAx", "bpmAy", "bpmBx", "bpmBy"};
  TString names[nVars] = {"acc0", "BCM", "lasPow", "bpmAx", "bpmAy", "bpmBx", "bpmBy"};
  for(Int_t i = 0; i < nVars; i++){
    breakdownPlots(prexOrCrex, runChains[mpswise], runNum, vars[i], "mpswise", names[i], maxEvt);
  }

  mpsGraphs(prexOrCrex, runChains[mpswise], runNum, "mpswise", maxEvt);
  quartetPlots(prexOrCrex, runChains[quartetwise], runNum, maxEvt);
  quartetGraphs(prexOrCrex, runChains[quartetwise], runNum, "quartetwise", maxEvt);

  gSystem->Exec(Form("convert %s/runs/Run%i/mps*.png %s/runs/Run%i/mpsCycCut.pdf", getenv("COMPMON_WEB"), runNum, getenv("COMPMON_WEB"), runNum));
  gSystem->Exec(Form("convert %s/runs/Run%i/qrt*.png %s/runs/Run%i/qrtCycCut.pdf", getenv("COMPMON_WEB"), runNum, getenv("COMPMON_WEB"), runNum));
  gSystem->Exec(Form("rm -f %s/runs/Run%i/mps*.png", getenv("COMPMON_WEB"), runNum));
  gSystem->Exec(Form("rm -f %s/runs/Run%i/qrt*.png", getenv("COMPMON_WEB"), runNum));
}
