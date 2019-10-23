#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

std::vector<TString> hist_stats(TH1F* h){
  std::vector<TString> stats;
  
  stats.push_back(Form("Entries: %i", (Int_t)h->GetEntries()));
  stats.push_back(Form("Mean: %.4f +/- %.4f", h->GetMean(), h->GetMeanError()));
  stats.push_back(Form("Std Dev: %.4f +/- %.4f", h->GetRMS(), h->GetRMSError()));

  return stats;
}

void getCanvasDims(int n_cycles, int *x){
  if(     n_cycles ==  1){x[0] = 1; x[1] = 1;}
  else if(n_cycles ==  2){x[0] = 2; x[1] = 1;}
  else if(n_cycles ==  3){x[0] = 1; x[1] = 3;}
  else if(n_cycles ==  4){x[0] = 2; x[1] = 2;}
  else if(n_cycles ==  5){x[0] = 2; x[1] = 3;}
  else if(n_cycles ==  6){x[0] = 2; x[1] = 3;}
  else if(n_cycles ==  7){x[0] = 2; x[1] = 4;}
  else if(n_cycles ==  8){x[0] = 2; x[1] = 4;}
  else if(n_cycles ==  9){x[0] = 3; x[1] = 3;}
  else if(n_cycles == 10){x[0] = 2; x[1] = 5;}
  else if(n_cycles == 11){x[0] = 3; x[1] = 4;}
  else if(n_cycles == 12){x[0] = 3; x[1] = 4;}
  else if(n_cycles == 13){x[0] = 4; x[1] = 4;}
  else if(n_cycles == 14){x[0] = 4; x[1] = 4;}
  else if(n_cycles == 15){x[0] = 4; x[1] = 4;}
  else if(n_cycles == 16){x[0] = 4; x[1] = 4;}
  else if(n_cycles == 17){x[0] = 3; x[1] = 6;}
  else if(n_cycles == 18){x[0] = 3; x[1] = 6;}
  else if(n_cycles == 19){x[0] = 4; x[1] = 5;}
  else if(n_cycles == 20){x[0] = 4; x[1] = 5;}
  else if(n_cycles == 21){x[0] = 4; x[1] = 6;}
  else if(n_cycles == 22){x[0] = 4; x[1] = 6;}
  else if(n_cycles == 23){x[0] = 4; x[1] = 6;}
  else if(n_cycles == 24){x[0] = 4; x[1] = 6;}
  else if(n_cycles >= 25){x[0] = 5; x[1] = 5;}
  else                   {x[0] = 0; x[1] = 0;}
}

vector<vector<int>> findDivisions(ifstream &infile, int run_num){
  string read_str;
  int miniruns_in_this_run = 0;
  vector<vector<int>> run;
  while(getline(infile, read_str)){
    vector<int> limits; stringstream ss(read_str);
    for(int i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    miniruns_in_this_run++;
    vector<int> minirun; minirun.push_back(run_num); minirun.push_back(miniruns_in_this_run);
    minirun.push_back(limits[0]); minirun.push_back(limits[1]);
    minirun.push_back(limits[2]); minirun.push_back(limits[3]);
    minirun.push_back(limits[4]); minirun.push_back(limits[5]);
    run.push_back(minirun);
  }
  return run;
}

void cyclePlotByRun(int runNum){
  gStyle->SetOptTitle(0);
  string run_num_str;
  ifstream infile(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), runNum));
  vector<vector<int>> cycles = findDivisions(infile, runNum);
  
  //TCanvas *c = new TCanvas(Form("c_%i", runNum), Form("Run %i Cycle Canvas", runNum), 1200, 1200);
  int n_cycles = cycles.size();
  int divs[2] = {-1, -1}; getCanvasDims(n_cycles, divs);
  //c->Divide(divs[0], divs[1]);
  TString pos("PosHelAcc0/PosHelNSamples0");
  TString neg("NegHelAcc0/NegHelNSamples0");
  TString msmts[5] = {pos, neg, Form("%s - %s", pos.Data(), neg.Data()), Form("%s + %s", pos.Data(), neg.Data()),
                      Form("(%s - %s)/(%s + %s)", pos.Data(), neg.Data(), pos.Data(), neg.Data())};
  TString ids[5] = {"PosHel", "NegHel", "Diffs", "Sums", "Asyms"};

  vector<TCanvas *> cans;
  for(int msmt = 0; msmt < 5; msmt++){
    TCanvas *c = new TCanvas(Form("c_%i_%i", runNum, msmt), Form("Run %i Cycle Canvas", runNum), 1200, 1200);
    c->Divide(divs[0], divs[1]);
    cans.push_back(c);
  }
  if(n_cycles > 25){n_cycles = 25;}
  TFile *rootfile = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *quartetwise = (TTree *)rootfile->Get("quartetwise");
  for(int msmt = 0; msmt < 5; msmt++){
    for(int div = 1; div <= n_cycles; div++){
      quartetwise->Draw(Form("%s>>h_%s_div%i_offStart", msmts[msmt].Data(), ids[msmt].Data(), div), 
                        Form("(laserState==2 || laserState==3) && beamState==1 && firstMPSnumber>=%i && firstMPSnumber<%i", cycles[div - 1][2], cycles[div - 1][3]), "goff");
      quartetwise->Draw(Form("%s>>h_%s_div%i_on", msmts[msmt].Data(), ids[msmt].Data(), div), 
                        Form("(laserState==0 || laserState==1) && beamState==1 && firstMPSnumber>=%i && firstMPSnumber<%i", cycles[div - 1][4], cycles[div - 1][5]), "goff");
      quartetwise->Draw(Form("%s>>h_%s_div%i_offEnd", msmts[msmt].Data(), ids[msmt].Data(), div), 
                        Form("(laserState==2 || laserState==3) && beamState==1 && firstMPSnumber>=%i && firstMPSnumber<%i", cycles[div - 1][6], cycles[div - 1][7]), "goff");
      cans[msmt]->cd(div);
      TPad *pName  = new TPad(Form("p_%s_div%i_name", msmts[msmt].Data(), div),  "Hist Title", 0.0, 0.9, 1.0, 1.0);
      TPad *pStat1 = new TPad(Form("p_%s_div%i_stat1", msmts[msmt].Data(), div), "Stats 1", 0.0, 0.7, 0.33, 0.9);
      TPad *pStat2 = new TPad(Form("p_%s_div%i_stat2", msmts[msmt].Data(), div), "Stats 2", 0.33, 0.7, 0.67, 0.9);
      TPad *pStat3 = new TPad(Form("p_%s_div%i_stat3", msmts[msmt].Data(), div), "Stats 3", 0.67, 0.7, 1.0, 0.9);
      TPad *pHist  = new TPad(Form("p_%s_div%i_hist", msmts[msmt].Data(), div),  "Histogram", 0.0, 0.0, 1.0, 0.7);
      pName->Draw(); pStat1->Draw(); pStat2->Draw(); pStat3->Draw(); pHist->Draw();

      TH1F *hOff1 = (TH1F *)gDirectory->Get(Form("h_%s_div%i_offStart", ids[msmt].Data(), div));
      TH1F *hOn   = (TH1F *)gDirectory->Get(Form("h_%s_div%i_on", ids[msmt].Data(), div));
      TH1F *hOff2 = (TH1F *)gDirectory->Get(Form("h_%s_div%i_offEnd", ids[msmt].Data(), div));

      TString title = Form("Quartetwise: %s, Cycle %i", ids[msmt].Data(), div);
      hOff1->SetTitle(title.Data()); hOn->SetTitle(title.Data()); hOff2->SetTitle(title.Data());
      hOff1->SetLineColor(kRed); hOn->SetLineColor(kGreen + 2); hOff2->SetLineColor(kRed + 2);

      THStack *hs = new THStack(Form("%s_cycle%i_stack", ids[msmt].Data(), div), 
                                Form("Run %i quartetwise, %s", runNum, ids[msmt].Data()));
      hs->Add(hOn); hs->Add(hOff1); hs->Add(hOff2);
      pHist->cd(); hs->Draw("nostack");

      std::vector<TString> sOff1 = hist_stats(hOff1); std::vector<TString> sOn = hist_stats(hOn); std::vector<TString> sOff2 = hist_stats(hOff2);
      TPaveText *ptOff1 = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptOff1->SetBorderSize(1); ptOff1->SetFillColor(0);
      TPaveText *ptOn   = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptOn->SetBorderSize(1);   ptOn->SetFillColor(0);
      TPaveText *ptOff2 = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptOff2->SetBorderSize(1); ptOff2->SetFillColor(0);
      for(int i = 0; i < sOff1.size(); i++){ptOff1->AddText(sOff1[i].Data())->SetTextColor(kRed);}
      for(int i = 0; i < sOn.size(); i++){ptOn->AddText(sOn[i].Data())->SetTextColor(kGreen + 2);}
      for(int i = 0; i < sOff2.size(); i++){ptOff2->AddText(sOff2[i].Data())->SetTextColor(kRed + 2);}
      pStat1->cd(); ptOff1->Draw(); 
      pStat2->cd(); ptOn->Draw(); 
      pStat3->cd(); ptOff2->Draw();
      
      TPaveText *ptName = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC"); ptName->SetBorderSize(0); ptName->SetFillColor(16);
      ptName->AddText(Form("Cycle %i, %s", div, ids[msmt].Data()));
      pName->cd(); ptName->Draw();
    }
    if(msmt == 0){cans[msmt]->Print(Form("%s/runs/Run%i/cycle_qVars.pdf(", getenv("COMPMON_WEB"), runNum), "pdf");}
    else if(msmt == 4){cans[msmt]->Print(Form("%s/runs/Run%i/cycle_qVars.pdf)", getenv("COMPMON_WEB"), runNum), "pdf");}
    else{cans[msmt]->Print(Form("%s/runs/Run%i/cycle_qVars.pdf", getenv("COMPMON_WEB"), runNum), "pdf");}
  }
}

void cyclePlot(TString group_fname){
  string run_num_str;
  ifstream infile(Form("%s/%s.list", getenv("COMPMON_SNAILS"), group_fname.Data()));
  while(getline(infile, run_num_str)){
    int run_num = atoi(run_num_str.c_str());
    cyclePlotByRun(run_num);
  }
  infile.close();
}
