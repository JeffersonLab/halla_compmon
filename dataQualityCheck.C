#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TChain.h>
#include <TColor.h>
#include <TString.h>

#include <vector>
#include <fstream>

#include "utils.h"

using namespace std;

/**vector<bool> accept(float snapshot[], int length, int minY, int maxY, int minX, int min_ped, int max_ped, float pedestal, int snapClock){
  float startSum = 0; float endSum = 0; int avgSamples = 20;
  for(int i = 40; i < 40 + avgSamples; i++){startSum += snapshot[i];}
  for(int i = length - 40 - 1; i > (length - 40 - 1 - avgSamples); i--){endSum += snapshot[i];}
  float startAvg = startSum*1.0/avgSamples; float endAvg = endSum*1.0/avgSamples;
  //bool pctDiff = abs(maxY - startAvg)/maxY < 0.03 && abs(maxY - endAvg)/maxY < 0.03;
  float background = (startAvg + endAvg)/2.0;
  bool peakHeight = abs(minY - background) > abs(maxY - background);
  bool saturated = (minY == 0);
  bool correct_pedestal = abs(pedestal - background)/pedestal < 0.03 && pedestal < 3900;
  bool narrow_pedestal = TMath::Abs(max_ped - min_ped) < 40;
  bool correct_time = snapClock>10e3 && snapClock<803e3;

  vector<bool> cuts; cuts.push_back(peakHeight); cuts.push_back(not saturated); 
  cuts.push_back(correct_pedestal); cuts.push_back(narrow_pedestal); cuts.push_back(correct_time);
  return cuts;
}**/

/**vector<int> get_stats(float snapshot[], int length){
  int minY = 1e6; int maxY = 0; int minX = 0; int sum = 0;
  int min_ped = 1e6; int max_ped  = 0; double pedestal = 0;
  for(int i = 40; i < length - 40; i++){
    if(snapshot[i] < minY){minY = snapshot[i]; minX = i;}
    if(snapshot[i] > maxY){maxY = snapshot[i];}
    if(i < 80){
      pedestal += snapshot[i];
      if(snapshot[i] < min_ped){min_ped = snapshot[i];}
      if(snapshot[i] > max_ped){max_ped = snapshot[i];}
    }
  }
  for(int i = 40; i < length - 40; i++){sum += pedestal/40.0 - snapshot[i];}
  vector<int> statsVec; statsVec.push_back(minY); statsVec.push_back(maxY);
  statsVec.push_back(minX); statsVec.push_back(sum); statsVec.push_back(min_ped);
  statsVec.push_back(max_ped); statsVec.push_back(pedestal*1000.0/40.0);
  return statsVec;
}**/

vector<bool> accept(float snapshot[], int length, int minY, int maxY, int minX, int min_ped, int max_ped, float pre_ped, float post_ped, int snapClock){
  double pedestal = (pre_ped + post_ped)/2.0;
  bool peakHeight = (pedestal - minY) > (maxY - pedestal);
  bool saturated = (minY == 0);
  bool correct_pedestal = TMath::Abs(pre_ped - post_ped)/pedestal < 0.03 && pedestal < 3900 && min_ped > minY;
  bool narrow_pedestal = TMath::Abs(max_ped - min_ped) < 50;
  bool correct_time = snapClock>10e3 && snapClock<803e3;

  //cout<<"Pedestal: "<<pedestal<<"; Pre Ped: "<<pre_ped<<"; Post Ped: "<<post_ped<<"; minY: "<<minY<<"; maxY: "<<maxY<<"; Min Ped: "<<
  //min_ped<<"; Max Ped: "<<max_ped<<endl;
  vector<bool> cuts; cuts.push_back(peakHeight); cuts.push_back(not saturated); 
  cuts.push_back(correct_pedestal); cuts.push_back(narrow_pedestal); cuts.push_back(correct_time);
  return cuts;
}

vector<int> get_stats(float snapshot[], int length){
  int minY = 1e6; int maxY = 0; int minX = 0; int sum = 0;
  int min_ped = 1e6; int max_ped  = 0; 
  int pedSamp = 40; double pre_ped = 0; double post_ped = 0;
  for(int i = pedSamp; i < 2*pedSamp; i++){
    pre_ped += snapshot[i];
    if(snapshot[i] < min_ped){min_ped = snapshot[i];}
    if(snapshot[i] > max_ped){max_ped = snapshot[i];}
  }
  for(int i = length - pedSamp - 1; i > length - 2*pedSamp - 1; i--){
    post_ped += snapshot[i];
    if(snapshot[i] < min_ped){min_ped = snapshot[i];}
    if(snapshot[i] > max_ped){max_ped = snapshot[i];}
  }
  for(int i = 2*pedSamp; i < length - 2*pedSamp - 1; i++){
    sum += (pre_ped + post_ped)*1.0/(2*pedSamp) - snapshot[i];
    if(snapshot[i] < minY){minY = snapshot[i]; minX = i;}
    if(snapshot[i] > maxY){maxY = snapshot[i];}
  }
  vector<int> statsVec; statsVec.push_back(minY); statsVec.push_back(maxY);
  statsVec.push_back(minX); statsVec.push_back(sum); statsVec.push_back(min_ped);
  statsVec.push_back(max_ped); statsVec.push_back(pre_ped*1000.0/pedSamp);
  statsVec.push_back(post_ped*1000.0/pedSamp);
  return statsVec;
}

inline bool file_exists(TString name){
  ifstream f(name.Data());
  return f.good();
}

vector<Double_t> count_acc0_limits(TTree *somethingwise, TString xvar){
  Double_t xmin = 1e10; Double_t xmax = -1e10; 
  Double_t ymin = 1e10; Double_t ymax = -1e10;
  Int_t x;
  Double_t acc0, nacc0;
  somethingwise->SetBranchAddress(xvar.Data(), &x);
  somethingwise->SetBranchAddress("Acc0", &acc0);
  somethingwise->SetBranchAddress("NAcc0", &nacc0);
  for(int i = 0; i < somethingwise->GetEntries(); i++){
    somethingwise->GetEntry(i);
    Float_t y = acc0/nacc0;
    if(x < xmin){xmin = x;}
    if(x > xmax){xmax = x;}
    if(y < ymin){ymin = y;}
    if(y > ymax){ymax = y;}
  }

  vector<Double_t> results; results.push_back(xmin); results.push_back(xmax);
  results.push_back(ymin); results.push_back(ymax);
  return results;
}

/**void make_histo(TString name, TString title, Float_t nbins, Float_t xmin, Float_t xmax, TString msmt){
  TH1F *h = (TH1F *)gDirectory->Get(name.Data());
  if(h){delete h;}
  h = new TH1F(name, title, nbins, xmin, xmax);
  h->GetXaxis()->SetTitle(msmt.Data());
}**/

void make_histo(TTree *somethingwise, TH1F *h, TString name, TString title, TString msmt, TString cuts){
  h = (TH1F *)gDirectory->Get(name.Data()); if(h){delete h;}
  somethingwise->Draw(Form("%s>>%s", msmt.Data(), name.Data()), cuts.Data(), "goff");
  h = (TH1F *)gDirectory->Get(name.Data()); h->SetTitle(title.Data());
  h->GetXaxis()->SetTitle(msmt.Data());
}

void make_graph(TString name){
  TGraph *g = (TGraph *)gDirectory->Get(name.Data());
  if(g){delete g;}
  g = new TGraph();
  g->SetName(name.Data());
}

void make_and_write_graph(TString name, TString title, vector<Float_t> xData, vector<Float_t> yData, TFile *outfile, Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, 
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
  g->Write(); delete g;
}

void write_histo(TString name, TFile *outfile){
  TH1F *h = (TH1F *)gDirectory->Get(name.Data());
  if(h){outfile->cd(); h->Write(); delete h;}
  else{printf("Could not write histogram %s", name.Data());}
}

void write_graph(TString name, TFile *outfile, TString title, Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax, TString msmtX, TString msmtY){
  TH2F *g = (TH2F *)gDirectory->Get(name.Data());
  Double_t xLoFac = 0.98; Double_t xHiFac = 1.02; Double_t yLoFac = 0.95; Double_t yHiFac = 1.05;
  if(xmin < 0){xLoFac = 1.02;} if(xmax < 0){xHiFac = 0.98;}
  if(ymin < 0){yLoFac = 1.05;} if(ymax < 0){yHiFac = 0.95;}
  if(g){
		outfile->cd(); 
    g->SetTitle(title.Data());
    g->GetXaxis()->SetTitle(msmtX.Data()); g->GetYaxis()->SetTitle(msmtY.Data());
    g->GetXaxis()->SetRangeUser(xmin*xLoFac, xmax*xHiFac); g->GetYaxis()->SetRangeUser(ymin*yLoFac, ymax*yHiFac);
    g->Write(); delete g;
  }
	else{printf("Could not write graph %s", name.Data());}
}

Float_t avg(vector<Float_t> data){
  Float_t tot = 0.0;
  for(Float_t datum : data){tot += datum;}
  return tot/((int)data.size());
}

int translate_day(TString day){
  if(day.EqualTo("Mon")){return 1;}
  else if(day.EqualTo("Tue")){return 2;}
  else if(day.EqualTo("Wed")){return 3;}
  else if(day.EqualTo("Thu")){return 4;}
  else if(day.EqualTo("Fri")){return 5;}
  else if(day.EqualTo("Sat")){return 6;}
  else if(day.EqualTo("Sun")){return 7;}
  else{return 0;}
}

int translate_month(TString month){
  if(month.EqualTo("Jan")){return 1;}
  else if(month.EqualTo("Feb")){return 2;}
  else if(month.EqualTo("Mar")){return 3;}
  else if(month.EqualTo("Apr")){return 4;}
  else if(month.EqualTo("May")){return 5;}
  else if(month.EqualTo("Jun")){return 6;}
  else if(month.EqualTo("Jul")){return 7;}
  else if(month.EqualTo("Aug")){return 8;}
  else if(month.EqualTo("Sep")){return 9;}
  else if(month.EqualTo("Oct")){return 10;}
  else if(month.EqualTo("Nov")){return 11;}
  else if(month.EqualTo("Dec")){return 12;}
  else{return 0;}
}

void epicsPlots(TChain *epicswise, int run_num, TFile *outfile){
  cout<<"Plotting epics times..."<<endl;

  string *date_ptr = 0;
  epicswise->SetBranchAddress("epics_datestring", &date_ptr);
  epicswise->GetEntry(0);
  
  TString date_str(*date_ptr); TObjArray *date_split = date_str.Tokenize(" ");
  TString day(        ( (TObjString *)(date_split->At(0)) )->String());
  TString month(      ( (TObjString *)(date_split->At(1)) )->String());
  TString date_num_str( ( (TObjString *)(date_split->At(2)) )->String());
  TString time_str    ( ( (TObjString *)(date_split->At(3)) )->String());
  TString year_str    ( ( (TObjString *)(date_split->At(5)) )->String());

  TObjArray *time_ptr = time_str.Tokenize(":");
  TString hour_str(   ( (TObjString *)(time_ptr->At(0)) )->String());
  TString minute_str( ( (TObjString *)(time_ptr->At(1)) )->String());
  TString second_str( ( (TObjString *)(time_ptr->At(2)) )->String());

  Int_t date_num = atoi(date_num_str.Data());
  Int_t year     = atoi(year_str.Data());
  Int_t hour     = atoi(hour_str.Data());
  Int_t minute   = atoi(minute_str.Data());
  Int_t second   = atoi(second_str.Data());

  TH1F *h_day      = new TH1F("h_day",      "Day of Week",    7,    1,    8); h_day->Fill(translate_day(day));
  TH1F *h_month    = new TH1F("h_month",    "Month of Year", 12,    1,   13); h_month->Fill(translate_month(month));
  TH1F *h_date_num = new TH1F("h_date_num", "Day Number",    31,    1,   32); h_date_num->Fill(date_num);
  TH1F *h_year     = new TH1F("h_year",     "Year",          50, 2001, 2051); h_year->Fill(year);
  TH1F *h_hour     = new TH1F("h_hour",     "Hour",          24,    0,   24); h_hour->Fill(hour);
  TH1F *h_minute   = new TH1F("h_minute",   "Minute",        60,    0,   60); h_minute->Fill(minute);
  TH1F *h_second   = new TH1F("h_second",   "Second",        60,    0,   60); h_second->Fill(second);

  write_histo("h_day",    outfile); write_histo("h_month", outfile); write_histo("h_date_num", outfile);
  write_histo("h_year",   outfile); write_histo("h_hour",  outfile); write_histo("h_minute",   outfile);
  write_histo("h_second", outfile);

  TString qw1_name("epics_qw1"); TString hw1_name("epics_hw1");
  epicswise->Project(qw1_name.Data(), qw1_name.Data(), "");
  epicswise->Project(hw1_name.Data(), hw1_name.Data(), "");
  TH1F *h_qw1 = (TH1F *)gDirectory->Get(qw1_name.Data());
  TH1F *h_hw1 = (TH1F *)gDirectory->Get(hw1_name.Data());
  write_histo(qw1_name, outfile);
  write_histo(hw1_name, outfile);
  TString beamCut("epics_bcm_average > 5");
  TString Ax_name("epics_bpmAx"); TString Ay_name("epics_bpmAy"); TString Bx_name("epics_bpmBx"); TString By_name("epics_bpmBy");
  TString vWien("epics_VWienAngle"); TString hWien("epics_HWienAngle"); TString solWien("epics_PhiFG");
  TString tableX("epics_tablePosX"); TString tableY("epics_tablePosY");
  TString targetPos("epics_targetPos");
  epicswise->Project(Ax_name.Data(), Ax_name.Data(), beamCut.Data()); epicswise->Project(Ay_name.Data(), Ay_name.Data(), beamCut.Data());
  epicswise->Project(Bx_name.Data(), Bx_name.Data(), beamCut.Data()); epicswise->Project(By_name.Data(), By_name.Data(), beamCut.Data());
  epicswise->Project(vWien.Data(), vWien.Data(), "epics_VWienAngle != 0"); epicswise->Project(hWien.Data(), hWien.Data(), "epics_HWienAngle != 0"); 
  epicswise->Project(solWien.Data(), solWien.Data(), "epics_PhiFG != 0");
  epicswise->Project(tableX.Data(), tableX.Data(), ""); epicswise->Project(tableY.Data(), tableY.Data(), "");
  epicswise->Project(targetPos.Data(), targetPos.Data(), "");
  write_histo(Ax_name, outfile); write_histo(Ay_name, outfile); write_histo(Bx_name, outfile); write_histo(By_name, outfile);
  if((TH1F *)gDirectory->Get(vWien.Data())) write_histo(vWien, outfile); 
  if((TH1F *)gDirectory->Get(hWien.Data())) write_histo(hWien, outfile); 
  if((TH1F *)gDirectory->Get(solWien.Data())) write_histo(solWien, outfile);
  write_histo(tableX, outfile); write_histo(tableY, outfile);
  if((TH1F *)gDirectory->Get(targetPos.Data())) write_histo(targetPos, outfile);
}

void epicsPlotsLast(TChain *epicswise, int run_num, TFile *outfile){
  cout<<"Plotting epics (last evt) times..."<<endl;

  string *date_ptr = 0;
  epicswise->SetBranchAddress("epics_datestring", &date_ptr);
  epicswise->GetEntry(epicswise->GetEntries()-1);
  
  TString date_str(*date_ptr); TObjArray *date_split = date_str.Tokenize(" ");
  TString day(        ( (TObjString *)(date_split->At(0)) )->String());
  TString month(      ( (TObjString *)(date_split->At(1)) )->String());
  TString date_num_str( ( (TObjString *)(date_split->At(2)) )->String());
  TString time_str    ( ( (TObjString *)(date_split->At(3)) )->String());
  TString year_str    ( ( (TObjString *)(date_split->At(5)) )->String());

  TObjArray *time_ptr = time_str.Tokenize(":");
  TString hour_str(   ( (TObjString *)(time_ptr->At(0)) )->String());
  TString minute_str( ( (TObjString *)(time_ptr->At(1)) )->String());
  TString second_str( ( (TObjString *)(time_ptr->At(2)) )->String());

  Int_t date_num = atoi(date_num_str.Data());
  Int_t year     = atoi(year_str.Data());
  Int_t hour     = atoi(hour_str.Data());
  Int_t minute   = atoi(minute_str.Data());
  Int_t second   = atoi(second_str.Data());

  TH1F *h_day      = new TH1F("h_day_last",      "Day of Week",    7,    1,    8); h_day->Fill(translate_day(day));
  TH1F *h_month    = new TH1F("h_month_last",    "Month of Year", 12,    1,   13); h_month->Fill(translate_month(month));
  TH1F *h_date_num = new TH1F("h_date_num_last", "Day Number",    31,    1,   32); h_date_num->Fill(date_num);
  TH1F *h_year     = new TH1F("h_year_last",     "Year",          50, 2001, 2051); h_year->Fill(year);
  TH1F *h_hour     = new TH1F("h_hour_last",     "Hour",          24,    0,   24); h_hour->Fill(hour);
  TH1F *h_minute   = new TH1F("h_minute_last",   "Minute",        60,    0,   60); h_minute->Fill(minute);
  TH1F *h_second   = new TH1F("h_second_last",   "Second",        60,    0,   60); h_second->Fill(second);

  write_histo("h_day_last",    outfile); write_histo("h_month_last", outfile); write_histo("h_date_num_last", outfile);
  write_histo("h_year_last",   outfile); write_histo("h_hour_last",  outfile); write_histo("h_minute_last",   outfile);
  write_histo("h_second_last", outfile);
}

void runwisePlots(TChain* runwise, int run_num, TFile* outfile){
  TString runwise_vars[7] = {"FADC_ped_value", "cavityPowerOnMin", "cavityPowerOffMax", 
                             "BCMOnMin", "BCMOffMax", "FADC_ithrnear", "FADC_ithrfar"};
  for(TString var : runwise_vars){
    runwise->Project(Form("runwise_%s", var.Data()), var.Data(), "");
    write_histo(Form("runwise_%s", var.Data()), outfile);
  }
}

void snapshotPlots(TChain* snapshots, int run_num, TFile* outfile, int max_event){
  cout<<"Plotting snapshots..."<<endl;

  float snapshot[300];
  //float bcm;
  int randomTime, mpsCoda, numSamples, snapClock;
  int rej_cuts[6] = {0, 0, 0, 0, 0, 0};
  int randoms = 0;

  snapshots->SetBranchAddress("randomTime", &randomTime);
  snapshots->SetBranchAddress("snap", &snapshot);
  snapshots->SetBranchAddress("mpsCoda", &mpsCoda);
  snapshots->SetBranchAddress("numSamples", &numSamples);
  snapshots->SetBranchAddress("snapClock", &snapClock);

  TString g_name("snapshots_sumPeakHeight"); TString title = Form("Run %i, snapshots: Peak Height vs. Sum", run_num);
  TString g_peds("snapshots_pedestals"); TString pedTitle = Form("Run %i, snapshots: Pedestal vs. Time", run_num);
  TString h_peds("snapshots_pedestals_histo");
  TString acc_name("snapshots_accepted"); TString acc_title = Form("Run %i, snapshots, Accepted", run_num);
  TString rej_names[5] = {"snapshots_rej_peak", "snapshots_rej_sat", "snapshots_rej_corrPed", 
                          "snapshots_rej_narrPed", "snapshots_rej_time"};
  TGraph *g_heightSum = (TGraph *)gDirectory->Get(g_name);
  TGraph *g_pedTime = (TGraph *)gDirectory->Get(g_peds);
  TH2F *h_pedTime = (TH2F *)gDirectory->Get(h_peds);
  if(g_heightSum){delete g_heightSum;}
  if(g_pedTime){delete g_pedTime;}
  if(h_pedTime){delete h_pedTime;}
  g_heightSum = new TGraph(); g_heightSum->SetName(g_name);
  g_pedTime = new TGraph(); g_pedTime->SetName(g_peds);
  h_pedTime = new TH2F(h_peds.Data(), pedTitle.Data(), 400, 0, (Int_t)snapshots->GetEntries(), 600, 3750, 3830);
  TGraph *g_acc = (TGraph *)gDirectory->Get(acc_name);
  if(g_acc){delete g_acc;}
  g_acc = new TGraph(); g_acc->SetName(acc_name.Data());
  TGraph *g_rejs[5];
  for(int i = 0; i < 5; i++){
    g_rejs[i] = (TGraph *)gDirectory->Get(rej_names[i]);
    if(g_rejs[i]){delete g_rejs[i];}
    g_rejs[i] = new TGraph(); g_rejs[i]->SetName(rej_names[i].Data());
  }

  int accepted = 0;
  for(int i = 0; i < snapshots->GetEntries(); i++){
    snapshots->GetEntry(i);
    vector<int> stats; stats = get_stats(snapshot, numSamples);
    int minY = stats[0]; int maxY = stats[1]; int minX = stats[2]; int sum = stats[3];
    int min_ped = stats[4]; int max_ped = stats[5]; double pre_ped = stats[6]*1.0/1000.0; double post_ped = stats[7]*1.0/1000.0;
    vector<bool> cuts = accept(snapshot, numSamples, minY, maxY, minX, min_ped, max_ped, pre_ped, post_ped, snapClock);
    bool allCuts = true;
    if(randomTime == 1){randoms++; continue;}
    for(int j = 0; j < 5; j++){
      //if(j >= 0){continue;}
      allCuts = allCuts && cuts[j];
      if(not cuts[j]){rej_cuts[j]++;}
    }
    if(mpsCoda > max_event || not allCuts){continue;}
    //if((pre_ped + post_ped)/2.0 > 3815){printf("Found snapshot with high pedestal: %.2f. It's at MPS %i\n", (pre_ped + post_ped)/2.0, mpsCoda);}
    g_pedTime->SetPoint(accepted, accepted, (pre_ped + post_ped)/2.0);
    h_pedTime->Fill(i, (pre_ped + post_ped)/2.0);
    g_heightSum->SetPoint(accepted, sum, (pre_ped + post_ped)/2.0 - minY);
    accepted++;
  }

  g_acc->SetPoint(0, 1, accepted*1.0/(snapshots->GetEntries() - randoms));

  outfile->cd();
  g_heightSum->SetTitle(title);
  g_heightSum->GetXaxis()->SetTitle("Sum"); g_heightSum->GetYaxis()->SetTitle("Peak Height");
  g_heightSum->Write();
  g_pedTime->SetTitle(pedTitle);
  g_pedTime->GetXaxis()->SetTitle("Entry"); g_pedTime->GetYaxis()->SetTitle("Pedestal Average");
  g_pedTime->Write();
  h_pedTime->GetXaxis()->SetTitle("Entry"); h_pedTime->GetYaxis()->SetTitle("Pedestal Average");
  h_pedTime->Write();
  g_acc->SetTitle(acc_title);
  g_acc->GetXaxis()->SetTitle("mpsCoda"); g_acc->GetYaxis()->SetTitle("Accepted?");
  g_acc->Write();

  for(int i = 0; i < 5; i++){
    g_rejs[i]->SetPoint(0, 1, rej_cuts[i]*1.0/(snapshots->GetEntries() - randoms));
    g_rejs[i]->SetTitle(rej_names[i]);
    g_rejs[i]->Write();
  }
}

void trig_sum_plots(TChain* somethingwise, int run_num, TFile* outfile, TString msmt, TString tree_id, TString msmt_id, int max_event){
  cout<<"Plotting "<<msmt_id.Data()<<"..."<<endl;
  
  TString states[2] = {"ON", "OFF"}; TString hels[2] = {"Neg", "Pos"};
  TString laserON("(laserState==0 || laserState == 1)"); TString laserOFF("(laserState==2 || laserState==3)");
  TString beamON("beamState==1"); TString beamOFF("beamState==0");
  TString hel0("helicityStateReported==0"); TString hel1("helicityStateReported==1");
  TString laserCuts[2] = {laserON, laserOFF}; TString beamCuts[2] = {beamON, beamOFF}; TString helCuts[2] = {hel0, hel1};
  TString extraCuts(Form("sumIsRandom==0 && abs(sumPre - 3790)<10 && sum>0 && sum<50e3 && mpsCoda<%i", max_event));

  TString name = Form("%s_%s", tree_id.Data(), msmt_id.Data());
  TString title = Form("Run %i %s, %s: All", run_num, tree_id.Data(), msmt_id.Data());
  TH1F *h; make_histo(somethingwise, h, name, title, msmt, extraCuts); write_histo(name, outfile);
  for(int i = 0; i < 2; i++){
    TString cuts1 = Form("%s && %s", beamCuts[i].Data(), extraCuts.Data());
    TString name1 = Form("%s_beam%i_%s", tree_id.Data(), (Int_t)!i, msmt_id.Data());
    TString title1 = Form("Run %i %s, %s: Beam %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data());
    TH1F *h1; make_histo(somethingwise, h1, name1, title1, msmt, cuts1); write_histo(name1, outfile);
    for(int j = 0; j < 2; j++){
      TString cuts2 = Form("%s && %s && %s", beamCuts[i].Data(), laserCuts[j].Data(), extraCuts.Data());
      TString name2 = Form("%s_beam%i_las%i_%s", tree_id.Data(), (Int_t)!i, (Int_t)!j, msmt_id.Data());
      TString title2 = Form("Run %i, %s, %s: Beam %s, Laser %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data(), states[j].Data());
      TH1F *h2; make_histo(somethingwise, h2, name2, title2, msmt, cuts2); write_histo(name2, outfile);
      TString cuts2A = Form("%s && %s", laserCuts[j].Data(), extraCuts.Data());
      TString name2A = Form("%s_las%i_%s", tree_id.Data(), (Int_t)!j, msmt_id.Data());
      TString title2A = Form("Run %i, %s, %s: Laser %s", run_num, tree_id.Data(), msmt_id.Data(), states[j].Data());
      TH1F *h2A; make_histo(somethingwise, h2A, name2A, title2A, msmt, cuts2A); write_histo(name2A, outfile);
		  for(int k = 0; k < 2; k++){
        TString cuts3 = Form("%s && %s && %s && %s", beamCuts[i].Data(), laserCuts[j].Data(), helCuts[k].Data(), extraCuts.Data());
			  TString name3 = Form("%s_beam%i_las%i_hel%i_%s", tree_id.Data(), (Int_t)!i, (Int_t)!j, k, msmt_id.Data()); 
        TString title3 = Form("Run %i %s, %s: Beam %s, Laser %s, Hel %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data(), states[j].Data(), hels[k].Data());
       	TH1F *h3; make_histo(somethingwise, h3, name3, title3, msmt, cuts3); write_histo(name3, outfile);
      }
    }
  }
}

void breakdown_plots(TChain* somethingwise, int run_num, TFile* outfile, TString msmt, TString tree_id, TString msmt_id, int max_event){
  cout<<"Plotting "<<msmt_id.Data()<<"..."<<endl;
  
  TString states[2] = {"ON", "OFF"}; TString hels[2] = {"Neg", "Pos"};
  TString laserON("(laserState==0 || laserState == 1)"); TString laserOFF("(laserState==2 || laserState==3)");
  TString beamON("beamState==1"); TString beamOFF("beamState==0");
  TString hel0("helicityStateReported==0"); TString hel1("helicityStateReported==1");
  TString laserCuts[2] = {laserON, laserOFF}; TString beamCuts[2] = {beamON, beamOFF}; TString helCuts[2] = {hel0, hel1};
  TString extraCuts(Form("mpsCoda<%i", max_event));

  TString name = Form("%s_%s", tree_id.Data(), msmt_id.Data());
  TString title = Form("Run %i %s, %s: All", run_num, tree_id.Data(), msmt_id.Data());
  TH1F *h; make_histo(somethingwise, h, name, title, msmt, extraCuts.Data()); write_histo(name, outfile);
  for(int i = 0; i < 2; i++){
    TString cuts1 = Form("%s && %s", beamCuts[i].Data(), extraCuts.Data());
    TString name1 = Form("%s_beam%i_%s", tree_id.Data(), (Int_t)!i, msmt_id.Data());
    TString title1 = Form("Run %i %s, %s: Beam %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data());
    TH1F *h1; make_histo(somethingwise, h1, name1, title1, msmt, cuts1); write_histo(name1, outfile);
    for(int j = 0; j < 2; j++){
      TString cuts2 = Form("%s && %s && %s", beamCuts[i].Data(), laserCuts[j].Data(), extraCuts.Data());
      TString name2 = Form("%s_beam%i_las%i_%s", tree_id.Data(), (Int_t)!i, (Int_t)!j, msmt_id.Data());
      TString title2 = Form("Run %i, %s, %s: Beam %s, Laser %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data(), states[j].Data());
      TH1F *h2; make_histo(somethingwise, h2, name2, title2, msmt, cuts2); write_histo(name2, outfile);
      TString cuts2A = Form("%s && %s", laserCuts[j].Data(), extraCuts.Data());
      TString name2A = Form("%s_las%i_%s", tree_id.Data(), (Int_t)!j, msmt_id.Data());
      TString title2A = Form("Run %i, %s, %s: Laser %s", run_num, tree_id.Data(), msmt_id.Data(), states[j].Data());
      TH1F *h2A; make_histo(somethingwise, h2A, name2A, title2A, msmt, cuts2A); write_histo(name2A, outfile);
		  for(int k = 0; k < 2; k++){
        TString cuts3 = Form("%s && %s && %s && %s", beamCuts[i].Data(), laserCuts[j].Data(), helCuts[k].Data(), extraCuts.Data());
			  TString name3 = Form("%s_beam%i_las%i_hel%i_%s", tree_id.Data(), (Int_t)!i, (Int_t)!j, k, msmt_id.Data()); 
        TString title3 = Form("Run %i %s, %s: Beam %s, Laser %s, Hel %s", run_num, tree_id.Data(), msmt_id.Data(), states[i].Data(), states[j].Data(), hels[k].Data());
        TH1F *h3; make_histo(somethingwise, h3, name3, title3, msmt, cuts3); write_histo(name3, outfile);
      }
    }
  }
}

void mps_graphs(TChain* somethingwise, int run_num, TFile* outfile, TString tree_id, int max_event){
  cout<<"Plotting mps graphs..."<<endl;

  //Initialize strings for plots
  TString names[8] = {"beam ON, laser ON, hel 0",  "beam ON, laser ON, hel 1",  "beam ON, laser OFF, hel 0",  "beam ON, laser OFF, hel 1",
                      "beam OFF, laser ON, hel 0", "beam OFF, laser ON, hel 1", "beam OFF, laser OFF, hel 0", "beam OFF, laser OFF, hel 1"};
  TString msmt_ids[7] = {"acc0_time", "bcm_time", "cavPower_time", "bpmAx_time", "bpmAy_time", "bpmBx_time", "bpmBy_time"};
  TString ytitles[7] = {"Acc0/NAcc0", "BCM", "Cav Power", "BPM Ax", "BPM Ay", "BPM Bx", "BPM By"};

  //Initialize Acc0, BPM, BCM plots y_vals[msmt_id][config]
  vector<vector<Float_t>> mpsCoda_vals; vector<vector<vector<Float_t>>> y_vals;
  for(int i = 0; i < 7; i++){
    vector<vector<Float_t>> real_vals;
    for(int j = 0; j < 8; j++){
      if(i == 0){vector<Float_t> mpsCoda_val; mpsCoda_vals.push_back(mpsCoda_val);}
      vector<Float_t> real_val; real_vals.push_back(real_val);
    }
    y_vals.push_back(real_vals);
  }

  Int_t mpsCoda, helicityStateReported, laserState, beamState;
  Double_t acc0, nacc0, bpmAx, bpmAy, bpmBx, bpmBy;
  Float_t bcm, cavPower;

  somethingwise->SetBranchAddress("mpsCoda", &mpsCoda); somethingwise->SetBranchAddress("helicityStateReported", &helicityStateReported);
  somethingwise->SetBranchAddress("laserState", &laserState); somethingwise->SetBranchAddress("beamState", &beamState);
  somethingwise->SetBranchAddress("Acc0", &acc0); somethingwise->SetBranchAddress("NAcc0", &nacc0);
  somethingwise->SetBranchAddress("bcm", &bcm); somethingwise->SetBranchAddress("cavPowerCalibrated", &cavPower);

  Float_t xmin = 1e16; Float_t xmax = -1e16;
  Float_t ymins[7] = { 1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16}; 
  Float_t ymaxs[7] = {-1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16};

  for(int i = 0; i < somethingwise->GetEntries(); i++){
    somethingwise->GetEntry(i);
    if(mpsCoda > max_event){continue;}
    Int_t laserInd = 0;
    Int_t beamInd = ((Int_t)!beamState)*4;
    //if(beamState == 2) beamInd = 4;
    if(beamState == 2) continue;
    if(laserState==2 || laserState==3){laserInd = 2;}
    mpsCoda_vals[beamInd + laserInd + helicityStateReported].push_back(i);
    if(mpsCoda < xmin){xmin = i;}
    if(mpsCoda > xmax){xmax = i;}
    Double_t all_vals[3] = {acc0/nacc0, bcm, cavPower};
    for(int j = 0; j < y_vals.size(); j++){
      y_vals[j][beamInd + laserInd + helicityStateReported].push_back(all_vals[j]);
      if(all_vals[j] < ymins[j]){ymins[j] = all_vals[j];}
      if(all_vals[j] > ymaxs[j]){ymaxs[j] = all_vals[j];}
    }
    if(beamState == 1){

    }
  }

  for(int i = 0; i < y_vals.size(); i++){
    for(int j = 0; j < mpsCoda_vals.size(); j++){
      Int_t beam = (Int_t)(j < 4);
      Int_t laser = (Int_t)(j == 0 or j == 1 or j == 4 or j == 5);
      Int_t hel = j % 2;
	    TString name = Form("%s_beam%i_las%i_hel%i_%s", tree_id.Data(), beam, laser, hel, msmt_ids[i].Data()); 
      TString title = Form("Run %i %s, %s", run_num, tree_id.Data(), names[j].Data());
      make_and_write_graph(name, title, mpsCoda_vals[j], y_vals[i][j], outfile, xmin, xmax, ymins[i], ymaxs[i], "mpsCoda", ytitles[i]);
    }
  }
}

void quartet_plots(TChain* somethingwise, int run_num, TFile* outfile, int max_event){
  cout<<"Plotting quartets..."<<endl;
  
  TString states[2] = {"ON", "OFF"};
  TString laserON("(laserState==0 || laserState == 1)"); TString laserOFF("(laserState==2 || laserState==3)");
  TString beamON("beamState==1"); TString beamOFF("beamState==0");
  TString laserCuts[2] = {laserON, laserOFF}; TString beamCuts[2] = {beamON, beamOFF};
  TString extraCuts(Form("firstMPSnumber < %i", max_event));

  TString posH0("PosHelAcc0/PosHelNSamples0");
  TString negH0("NegHelAcc0/NegHelNSamples0");
  TString diff0 = Form("%s - %s", posH0.Data(), negH0.Data());
  TString summ0 = Form("%s + %s", posH0.Data(), negH0.Data());
  TString asym0 = Form("(%s)/(%s)", diff0.Data(), summ0.Data());
  TString posH4("PosHelAcc4");
  TString negH4("NegHelAcc4");
  TString diff4 = Form("%s - %s", posH4.Data(), negH4.Data());
  TString summ4 = Form("%s + %s", posH4.Data(), negH4.Data());
  TString asym4 = Form("(%s)/(%s)", diff4.Data(), summ4.Data());
  TString asymUSbg1("(PosHelScalerRun4  - NegHelScalerRun4) /(PosHelScalerRun4  + NegHelScalerRun4)");
  TString asymUSbg2("(PosHelScalerRun13 - NegHelScalerRun13)/(PosHelScalerRun13 + NegHelScalerRun13)");
  TString asymDSbg1("(PosHelScalerRun0  - NegHelScalerRun0) /(PosHelScalerRun0 + NegHelScalerRun0)");
  TString asymDSbg2("(PosHelScalerRun1  - NegHelScalerRun1) /(PosHelScalerRun1 + NegHelScalerRun1)");
  TString asymHoriz("(PosHelScalerRun10 - NegHelScalerRun10)/(PosHelScalerRun10 + NegHelScalerRun10)");
  TString asymVert ("(PosHelScalerRun11 - NegHelScalerRun11)/(PosHelScalerRun11 + NegHelScalerRun11)");

  TString meas[20] = {posH0, negH0, diff0, summ0, asym0, posH4, negH4, diff4, summ4, asym4, asymUSbg1, asymUSbg2, asymDSbg1, asymDSbg2, asymHoriz, asymVert, "", "", "", ""};
  TString meas_id[20] = {"PosHelAcc0/NSamples0", "NegHelAcc0/NSamples0", "Diffs0", "Sums0", "Asyms0 (Raw)", "PosHelAcc4", "NegHelAcc4", "Diffs4", "Sums4", "Asyms4 (Raw)",
                         "Asyms USbg1", "Asyms USbg2", "Asyms DSbg1", "Asyms DSbg2", "Asyms Horiz Finger", "Asyms Vert Finger", "Asyms0", "Asyms4", "Asyms0", "Asyms4"};
  TString meas_small[20] = {"posH0", "negH0", "diff0", "summ0", "asym0", "posH4", "negH4", "diff4", "summ4", "asym4", 
                            "USbg1", "USbg2", "DSbg1", "DSbg2", "horizFing", "vertFing", "asym0subOn",  "asym4subOn", "asym0subOff", "asym4subOff"};
  TString tree_id("quartetwise");

  Float_t summOn0 = 0.0; Float_t summOff0 = 0.0;
  Float_t summOn4 = 0.0; Float_t summOff4 = 0.0;
  for(int i = 0; i < 20; i++){
    if(i==16){
      meas[i]   = Form("1000*(%s)/(%s - %f)", diff0.Data(), summ0.Data(), summOff0);
      meas[i+1] = Form("1000*(%s)/(%s - %f)", diff4.Data(), summ4.Data(), summOff4);
      meas[i+2] = Form("1000*(%s)/(%f - %f)", diff0.Data(), summOn0, summOff0);
      meas[i+3] = Form("1000*(%s)/(%f - %f)", diff4.Data(), summOn4, summOff4);
    }
    TString hname = Form("%s_%s", tree_id.Data(), meas_small[i].Data());
    TString title = Form("Run %i %s, %s: All", run_num, tree_id.Data(), meas_id[i].Data());
    TH1F *h; make_histo(somethingwise, h, hname, title, meas[i], extraCuts.Data()); write_histo(hname, outfile);
    for(int j = 0; j < 2; j++){
      TString cuts1 = Form("%s && %s", beamCuts[j].Data(), extraCuts.Data());
      TString hname1 = Form("%s_beam%i_%s", tree_id.Data(), (Int_t)!j, meas_small[i].Data());
      TString title1 = Form("Run %i %s, %s: Beam %s", run_num, tree_id.Data(), meas_id[i].Data(), states[j].Data());
      TH1F *h1; make_histo(somethingwise, h1, hname1, title1, meas[i], cuts1); write_histo(hname1, outfile);
      for(int k = 0; k < 2; k++){
        TString cuts2 = Form("%s && %s", laserCuts[k].Data(), extraCuts.Data());
        TString hname2 = Form("%s_las%i_%s", tree_id.Data(), (Int_t)!k, meas_small[i].Data());
        TString title2 = Form("Run %i %s, %s: Laser %s", run_num, tree_id.Data(), meas_id[i].Data(), states[k].Data());
        TH1F *h2; make_histo(somethingwise, h2, hname2, title2, meas[i], cuts2); write_histo(hname2, outfile);
        TString cuts3 = Form("%s && %s && %s", beamCuts[j].Data(), laserCuts[k].Data(), extraCuts.Data());
        TString hname3 = Form("%s_beam%i_las%i_%s", tree_id.Data(), (Int_t)!j, (Int_t)!k, meas_small[i].Data());
        TString title3 = Form("Run %i %s, %s: Beam %s, Laser %s", run_num, tree_id.Data(), meas_id[i].Data(), states[j].Data(), states[k].Data());
        TH1F *h3; make_histo(somethingwise, h3, hname3, title3, meas[i], cuts3); h3 = (TH1F *)gDirectory->Get(hname3.Data());
        if(i==3 && j==0 && k==0){summOn0 = h3->GetMean(); meas[i+1] = Form("1000*(%s)/(%s - %f)", diff0.Data(), summ0.Data(), summOn0);}
        if(i==8 && j==0 && k==0){summOn4 = h3->GetMean(); meas[i+1] = Form("1000*(%s)/(%s - %f)", diff4.Data(), summ4.Data(), summOn4);}
        if(i==3 && j==0 && k==1){summOff0 = h3->GetMean();}
        if(i==8 && j==0 && k==1){summOff4 = h3->GetMean();}
        write_histo(hname3, outfile);
      }
    }
  }
}

void quartet_graphs(TChain* somethingwise, int run_num, TFile* outfile, TString tree_id, int max_event){
  cout<<"Plotting quartets vs time..."<<endl;

  TString states[2] = {"OFF", "ON"};
  Int_t n_msmts = 20;
  TString msmts[20] = {"posH0", "negH0", "diff0", "summ0", "asym0", "bcm", "posH4", "negH4", "diff4", "summ4", "asym4",
                       "USbg1", "USbg2", "DSbg1", "DSbg2", "vertFing", "horizFing", "cavPow", "asym0sub", "asym4sub"};

  vector<vector<vector<Float_t>>> x_vals; vector<vector<vector<Float_t>>> y_vals;
  for(int i = 0; i < n_msmts; i++){
    vector<vector<Float_t>> x; vector<vector<Float_t>> y;
    vector<Float_t> x_l0_b0; vector<Float_t> x_l0_b1; vector<Float_t> x_l1_b0; vector<Float_t> x_l1_b1;
    vector<Float_t> y_l0_b0; vector<Float_t> y_l0_b1; vector<Float_t> y_l1_b0; vector<Float_t> y_l1_b1;
    x.push_back(x_l0_b0); x.push_back(x_l0_b1); x.push_back(x_l1_b0); x.push_back(x_l1_b1);
    y.push_back(y_l0_b0); y.push_back(y_l0_b1); y.push_back(y_l1_b0); y.push_back(y_l1_b1);
    x_vals.push_back(x); y_vals.push_back(y);
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

  //Double_t posH_ymin =  5.5; Double_t posH_ymax = 12.0;
  //Double_t negH_ymin =  5.5; Double_t negH_ymax = 12.0;
  //Double_t diff_ymin = -2.0; Double_t diff_ymax = 2.0;
  //Double_t summ_ymin = 11.0; Double_t summ_ymax = 22.0;
  //Double_t asym_ymin = -0.2; Double_t asym_ymax = 0.2;
  //Double_t bcm_ymin = 0.0;   Double_t bcm_ymax = 120.0;

  //Double_t ymins[6] = {posH_ymin, negH_ymin, diff_ymin, summ_ymin, asym_ymin, bcm_ymin};
  //Double_t ymaxs[6] = {posH_ymax, negH_ymax, diff_ymax, summ_ymax, asym_ymax, bcm_ymax};
  Double_t ymins[18] = { 1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  
                         1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16,  1e16};
  Double_t ymaxs[18] = {-1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, 
                        -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16, -1e16};
  for(int i = 0; i < somethingwise->GetEntries(); i++){
    somethingwise->GetEntry(i);
    if(firstMPSnumber > max_event){continue;}
    Int_t laserInd = (Int_t)(laserState==0 || laserState==1);
    Int_t beamInd  = (Int_t)(beamState==1);
    Double_t msmt[18] = {posAcc/posSamp, negAcc/negSamp, posAcc/posSamp - negAcc/negSamp, posAcc/posSamp + negAcc/negSamp, 
                        1000*(posAcc/posSamp - negAcc/negSamp)/(posAcc/posSamp + negAcc/negSamp), (posBCM + negBCM)/2.0, 
                         posAcc4, negAcc4, posAcc4 - negAcc4, posAcc4 + negAcc4, 1000*(posAcc4 - negAcc4)/(posAcc4 + negAcc4),
                        (posUSbg1 + negUSbg1)*30.0, (posUSbg2 + negUSbg2)*30.0, (posDSbg1 + negDSbg1)*30.0, (posDSbg2 + negDSbg2)*30.0,
                        (posHoriz + negHoriz)*30.0, (posVert + negVert)*30.0, (posCavPow + negCavPow)/2.0};
    for(int j = 0; j < n_msmts; j++){
      if(msmt[j] < ymins[j]){ymins[j] = msmt[j];}
      if(msmt[j] > ymaxs[j]){ymaxs[j] = msmt[j];}
      x_vals[j][2*laserInd + beamInd].push_back(i);
      y_vals[j][2*laserInd + beamInd].push_back(msmt[j]);
    }
  }

  Float_t summOn0  = avg(y_vals[3][3]); Float_t yOnMin0  = 1e16; Float_t yOnMax0  = -1e16;
  Float_t summOff0 = avg(y_vals[3][1]); Float_t yOffMin0 = 1e16; Float_t yOffMax0 = -1e16;
  Float_t summOn4  = avg(y_vals[9][3]); Float_t yOnMin4  = 1e16; Float_t yOnMax4  = -1e16;
  Float_t summOff4 = avg(y_vals[9][1]); Float_t yOffMin4 = 1e16; Float_t yOffMax4 = -1e16;
  vector<Float_t> asym0sub; vector<Float_t> asym0subOff;
  vector<Float_t> asym4sub; vector<Float_t> asym4subOff; 
  for(int i = 0; i < y_vals[4][3].size(); i++){
    y_vals[4][3][i]  = 1000*y_vals[2][3][i]/(y_vals[3][3][i] - summOff0);
    y_vals[10][3][i] = 1000*y_vals[8][3][i]/(y_vals[9][3][i] - summOff4);
    Float_t asym0 = 1000*y_vals[2][3][i]/(summOn0 - summOff0);
    Float_t asym4 = 1000*y_vals[8][3][i]/(summOn4 - summOff4);
    if(asym0 < yOnMin0){yOnMin0 = asym0;} if(asym0 > yOnMax0){yOnMax0 = asym0;}
    if(asym4 < yOnMin4){yOnMin4 = asym4;} if(asym4 > yOnMax4){yOnMax4 = asym4;}
    asym0sub.push_back(asym0);
    asym4sub.push_back(asym4);
    
  }
  for(int i = 0; i < y_vals[4][1].size(); i++){
    Float_t asym0 = 1000*y_vals[2][1][i]/(summOn0 - summOff0);
    Float_t asym4 = 1000*y_vals[8][1][i]/(summOn4 - summOff4);
    if(asym0 < yOffMin0){yOffMin0 = asym0;} if(asym0 > yOffMax0){yOffMax0 = asym0;}
    if(asym4 < yOffMax4){yOffMax4 = asym4;} if(asym4 > yOffMax4){yOffMax4 = asym4;}
    asym0subOff.push_back(asym0);
    asym4subOff.push_back(asym4);
  }

  for(int i = 0; i < x_vals.size(); i++){
    for(int las = 0; las < 2; las++){
      for(int beam = 0; beam < 2; beam++){
        TString name = Form("%s_beam%i_las%i_%s_time", tree_id.Data(), beam, las, msmts[i].Data());
        TString title = Form("Run %i %s, %s vs time: beam %s, laser %s", run_num, tree_id.Data(), msmts[i].Data(), states[beam].Data(), states[las].Data());
        make_and_write_graph(name, title, x_vals[i][2*las+beam], y_vals[i][2*las+beam], outfile, 0, somethingwise->GetEntries(), ymins[i], ymaxs[i], "Entry$", msmts[i]);
      }
    }
  }
  TString nameOn0 = Form("%s_beam1_las1_asym0subOn_time", tree_id.Data());
  TString titleOn0 = Form("Run %i %s, Asym0 vs time: beam ON, laser ON", run_num, tree_id.Data());
  make_and_write_graph(nameOn0, titleOn0, x_vals[4][3], asym0sub, outfile, 0, somethingwise->GetEntries(), yOnMin0, yOnMax0, "Entry$", "asym (Acc0) / ppt");
  TString nameOff0 = Form("%s_beam1_las0_asym0subOff_time", tree_id.Data());
  TString titleOff0 = Form("Run %i %s, Asym0 vs time: beam ON, laser OFF", run_num, tree_id.Data());
  make_and_write_graph(nameOff0, titleOff0, x_vals[4][1], asym0subOff, outfile, 0, somethingwise->GetEntries(), yOffMin0, yOffMax0, "Entry$", "asym (Acc0) / ppt");
  TString nameOn4 = Form("%s_beam1_las1_asym4subOn_time", tree_id.Data());
  TString titleOn4 = Form("Run %i %s, Asym4 vs time: beam ON, laser ON", run_num, tree_id.Data());
  make_and_write_graph(nameOn4, titleOn4, x_vals[4][3], asym4sub, outfile, 0, somethingwise->GetEntries(), yOnMin4, yOnMax4, "Entry$", "asym (Acc4) / ppt");
  TString nameOff4 = Form("%s_beam1_las0_asym4subOff_time", tree_id.Data());
  TString titleOff4 = Form("Run %i %s, Asym4 vs time: beam ON, laser OFF", run_num, tree_id.Data());
  make_and_write_graph(nameOff4, titleOff4, x_vals[4][1], asym4subOff, outfile, 0, somethingwise->GetEntries(), yOffMin4, yOffMax4, "Entry$", "asym (Acc4) / ppt");
}



void dataQualityCheck(int run_num, int max_evt=1e9){
  TString compmon_out_path(getenv("COMPMON_PLOTFILES"));
  //TString fname = Form("%s/compmon_%i.root", rootfiles_path.Data(), run_num);
  TString outname = Form("%s/compton_online_run_%i.root", compmon_out_path.Data(), run_num);
  //TChain *snapshots = new TChain("snapshots");
  //TChain *triggerwise = new TChain("triggerwise");
  //TChain *mpswise = new TChain("mpswise");
  //TChain *quartetwise = new TChain("quartetwise");
  //TChain *epicswise = new TChain("epicswise");
  //TChain *runwise = new TChain("runwise");
  //snapshots->Add(fname.Data());
  //triggerwise->Add(fname.Data());
  //mpswise->Add(fname.Data());
  //quartetwise->Add(fname.Data());
  //epicswise->Add(fname.Data());
  //runwise->Add(fname.Data());
  Int_t mpswise = 0; Int_t quartetwise = 1;
  Int_t pulserwise = 2; Int_t triggerwise = 3;
  Int_t epicswise = 4; Int_t runwise = 5;
  Int_t snapshots = 6;
  vector<TChain *> runChains = loadChain(run_num);
  TFile *fout = new TFile(outname.Data(), "RECREATE");
  epicsPlots(runChains[epicswise], run_num, fout);
  epicsPlotsLast(runChains[epicswise], run_num, fout);
  runwisePlots(runChains[runwise], run_num, fout);
  snapshotPlots(runChains[snapshots], run_num, fout, max_evt);
  trig_sum_plots(runChains[triggerwise], run_num, fout, "sum", "triggerwise", "sums", max_evt);
  breakdown_plots(runChains[mpswise], run_num, fout, "Acc0/NAcc0", "mpswise", "acc0", max_evt);
  breakdown_plots(runChains[mpswise], run_num, fout, "numTriggers", "mpswise", "trigs", max_evt);
  mps_graphs(runChains[mpswise], run_num, fout, "mpswise", max_evt);
  quartet_plots(runChains[quartetwise], run_num, fout, max_evt);
  quartet_graphs(runChains[quartetwise], run_num, fout, "quartetwise", max_evt);

  TString hname("mpswise_ihwp"); TString htitle = Form("Run %i mpswise, IHWP State, all evts", run_num);
  TH1F *h_ihwp = new TH1F(hname.Data(), htitle.Data(), 2, 0, 2);
  runChains[mpswise]->Project(hname.Data(), "epics_ihwp_in", "mpsCoda > 300");
  fout->cd(); h_ihwp->Write(); delete h_ihwp;
  
  fout->Close();
}
