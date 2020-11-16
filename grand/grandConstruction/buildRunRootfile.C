#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "buildGrandRootfile.h"
#include "plot.h"

using namespace std;

/**
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

void epicsPlots(vector<vector<TString>> runEpcVars, Int_t runNum, TTree *epicswise, TFile *outfile){
  printf("  Plotting run epics vars...\n");
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

  TH1F *h_month    = new TH1F(Form("h%i_month", runNum),    "Month of Year", 12,    1,   13); h_month->Fill(translate_month(month));
  TH1F *h_date_num = new TH1F(Form("h%i_date_num", runNum), "Day Number",    31,    1,   32); h_date_num->Fill(date_num);
  TH1F *h_year     = new TH1F(Form("h%i_year", runNum),     "Year",          50, 2001, 2051); h_year->Fill(year);
  TH1F *h_hour     = new TH1F(Form("h%i_hour", runNum),     "Hour",          24,    0,   24); h_hour->Fill(hour);
  TH1F *h_minute   = new TH1F(Form("h%i_minute", runNum),   "Minute",        60,    0,   60); h_minute->Fill(minute);
  TH1F *h_second   = new TH1F(Form("h%i_second", runNum),   "Second",        60,    0,   60); h_second->Fill(second);

  outfile->cd(); h_year->Write(); h_month->Write(); h_date_num->Write();
  h_hour->Write(); h_minute->Write(); h_second->Write();

  for(Int_t i = 0; i < runEpcVars.size(); i++){
    TString hName = Form("h%i_%s", runNum, runEpcVars[i][0].Data());
    epicswise->Project(hName.Data(), runEpcVars[i][1].Data(), runEpcVars[i][2].Data());
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    outfile->cd(); h->Write();
  }
}
**/

void buildRunRootfile(Int_t runNum){
  TFile *runOut = new TFile(Form("%s/Run%i_Plots.root", getenv("COMPMON_RUNPLOTS"), runNum), "RECREATE");
  //TFile *infile = new TFile(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum), "READ");
  vector<TChain *> runChains = loadChain(runNum);
  TChain *mpswise = runChains[0]; TChain *quartetwise = runChains[1];
  TChain *epicswise = runChains[4]; TChain *triggerwise = runChains[3];

  //vector<vector<TString>> cycMPSVars = readVarsFile("cyc", "mpswise");
  //vector<vector<TString>> cycQrtVars = readVarsFile("cyc", "quartetwise");
  //vector<vector<TString>> runMPSVars = readVarsFile("run", "mpswise");
  //vector<vector<TString>> runEpcVars = readVarsFile("run", "epicswise");

  //TTree *mpswise = (TTree *)gDirectory->Get("mpswise");
  //TTree *quartetwise = (TTree *)gDirectory->Get("quartetwise");
  //TTree *epicswise = (TTree *)gDirectory->Get("epicswise");
  //TTree *triggerwise = (TTree *)gDirectory->Get("triggerwise");
  vector<vector<int>> cycles = findCycles(runNum);
  for(Int_t c = 0; c < cycles.size(); c++){
    printf("Plotting cycle %i/%i...\n", c + 1, (Int_t)cycles.size());
    //Int_t firstOffStartMPS = cycles[c][0]; Int_t firstOffEndMPS = cycles[c][1];
    //Int_t onStartMPS = cycles[c][2]; Int_t onEndMPS = cycles[c][3];
    //Int_t lastOffStartMPS = cycles[c][4]; Int_t lastOffEndMPS = cycles[c][5];
    //TString cycCut = Form("((mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i))",
    //                    firstOffStartMPS, firstOffEndMPS, onStartMPS, onEndMPS, lastOffStartMPS, lastOffEndMPS);
    printf("  Plotting mpswise vars...\n");
    cycMPSPlots(mpswise, runNum, c+1, cycles[c], runOut);
    /**
    for(Int_t i = 0; i < cycMPSVars.size(); i++){
      //printf("  Ploting %s...\n", cycMPSVars[i][0].Data());
      TString hName = Form("h%i.%i_%s", runNum, c + 1, cycMPSVars[i][0].Data());
      mpswise->Draw(Form("%s>>%s", cycMPSVars[i][1].Data(), hName.Data()), Form("(%s) && (%s)", cycMPSVars[i][2].Data(), cycCut.Data()), "goff");
      TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
      runOut->cd(); h->Write();
    }
    **/
    printf("  Plotting quartetwise vars...\n");
    cycQrtPlots(quartetwise, runNum, c+1, cycles[c], runOut);
    /**
    for(Int_t i = 0; i < cycQrtVars.size(); i++){
      //printf("  Ploting %s...\n", cycQrtVars[i][0].Data());
      TString hName = Form("h%i.%i_%s", runNum, c + 1, cycQrtVars[i][0].Data());
      quartetwise->Draw(Form("%s>>%s", cycQrtVars[i][1].Data(), hName.Data()), Form("(%s) && (%s)", cycQrtVars[i][2].Data(), cycCut.Data()), "goff");
      TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
      runOut->cd(); h->Write();
    }
    **/
    printf("  Plotting triggerwise vars...\n");
    cycTrgPlots(triggerwise, runNum, c+1, cycles[c], runOut);
  }
  printf("  Plotting run mps vars...\n");
  runMPSPlots(mpswise, runNum, runOut);
  runEpicsPlots(epicswise, mpswise, runNum, runOut);
  /**
  for(Int_t i = 0; i < runMPSVars.size(); i++){
    TString hName = Form("h%i_%s", runNum, runMPSVars[i][0].Data());
    mpswise->Draw(Form("%s>>%s", runMPSVars[i][1].Data(), hName.Data()), runMPSVars[i][2].Data(), "goff");
    TH1F *h = (TH1F *)gDirectory->Get(hName.Data());
    runOut->cd(); h->Write();
  }
  epicsPlots(runEpcVars, runNum, epicswise, runOut);
  **/

  //runOut->Write(); 
  runOut->Close();
  printf("...Done!\n");
}
