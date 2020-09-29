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

using namespace std;

void makeRunPedestalPlots(Int_t runNum, vector<vector<Int_t>> cycles){
  printf("Plotting pedestals for run %i...\n", runNum);
  TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *snapshots = (TTree *)f->Get("snapshots");

  float snapshot[300];
  int randomTime, mpsCoda, numSamples, snapClock, beamState, laserState;
  int rej_cuts[6] = {0, 0, 0, 0, 0, 0};
  int randoms = 0;
  TH1F *hON = new TH1F(Form("hON_run%i", runNum), Form("Run %i Pedestal Mean", runNum), 240, 3770, 3800);
  hON->SetLineColor(kGreen + 2); hON->SetStats(0); 
  hON->GetXaxis()->SetTitle("Pedestal Mean [RAU]");
  TH1F *hOFF = new TH1F(Form("hOFF_run%i", runNum), Form("Run %i Pedestal Mean", runNum), 240, 3770, 3800);
  hOFF->SetLineColor(kRed); hOFF->SetStats(0);
  hOFF->GetXaxis()->SetTitle("Pedestal Mean [RAU]");

  snapshots->SetBranchAddress("randomTime", &randomTime);
  snapshots->SetBranchAddress("snap", &snapshot);
  snapshots->SetBranchAddress("mpsCoda", &mpsCoda);
  snapshots->SetBranchAddress("numSamples", &numSamples);
  snapshots->SetBranchAddress("snapClock", &snapClock);
  snapshots->SetBranchAddress("beamState", &beamState);
  snapshots->SetBranchAddress("laserState", &laserState);
  
  for(Int_t i = 0; i < snapshots->GetEntries(); i++){
    snapshots->GetEntry(i);
    bool inRange = false;
    bool limitExceeded = false;
    Float_t sum = 0;
    for(Int_t j = 0; j < 8; j++){
      if(mpsCoda >= cycles[j][0] && mpsCoda <= cycles[j][1]){
        inRange = true;
      }
    }
    for(Int_t j = 0; j < 40; j++){
      sum += snapshot[j];
      if(snapshot[j] < 3770) limitExceeded = true;
    }
    if(!inRange || limitExceeded || beamState != 1) continue;
    sum /= 40.0;
    if(laserState == 0 || laserState == 1) hON->Fill(sum);
    else if(laserState == 2 || laserState == 3) hOFF->Fill(sum);
  }

  THStack *hs = new THStack(Form("hsRun%i", runNum), Form("Run %i Pedestal Mean", runNum));
  hs->Add(hON); hs->Add(hOFF);

  TCanvas *c = new TCanvas(Form("cRun%i", runNum), Form("Run %i Pedestals", runNum), 1200, 800);
  c->cd();
  //hON->Draw();
  //hOFF->Draw("same");
  hs->Draw("nostack");

  TPaveText *ptON = new TPaveText(0.78, 0.82, 0.9, 0.9, "blNDC");
  ptON->SetFillColor(0); ptON->SetBorderSize(1);
  TPaveText *ptOFF = new TPaveText(0.78, 0.74, 0.9, 0.82, "blNDC");
  ptOFF->SetFillColor(0); ptOFF->SetBorderSize(1);
  ptON->AddText("----Laser ON----")->SetTextColor(kGreen + 2);
  ptON->AddText(Form("Mean: %.2f +/- %.2f", hON->GetMean(), hON->GetMeanError()))->SetTextColor(kGreen + 2);
  ptOFF->AddText("----Laser OFF----")->SetTextColor(kRed);
  ptOFF->AddText(Form("Mean: %.2f +/- %.2f", hOFF->GetMean(), hOFF->GetMeanError()))->SetTextColor(kRed);
  ptON->Draw("same");
  ptOFF->Draw("same");
}

void snapshotPedestal(){
  Int_t runNum1 = 4364; //Specifically chosen for stability
  Int_t runNum2 = 4614;
  Int_t runNum3 = 4554;
  Int_t runNum4 = 4430;
  Int_t runNum5 = 4345;

  vector<vector<Int_t>> run1;
  vector<Int_t> cycle1_1; cycle1_1.push_back(150833); cycle1_1.push_back(208232); run1.push_back(cycle1_1);
  vector<Int_t> cycle1_2; cycle1_2.push_back(193601); cycle1_2.push_back(250160); run1.push_back(cycle1_2);
  vector<Int_t> cycle1_3; cycle1_3.push_back(237841); cycle1_3.push_back(291808); run1.push_back(cycle1_3);
  vector<Int_t> cycle1_4; cycle1_4.push_back(279649); cycle1_4.push_back(334080); run1.push_back(cycle1_4);
  vector<Int_t> cycle1_5; cycle1_5.push_back(150833); cycle1_5.push_back(208232); run1.push_back(cycle1_5);
  vector<Int_t> cycle1_6; cycle1_6.push_back(193601); cycle1_6.push_back(250160); run1.push_back(cycle1_6);
  vector<Int_t> cycle1_7; cycle1_7.push_back(237841); cycle1_7.push_back(291808); run1.push_back(cycle1_7);
  vector<Int_t> cycle1_8; cycle1_8.push_back(279649); cycle1_8.push_back(334080); run1.push_back(cycle1_8);

  vector<vector<Int_t>> run2;
  vector<Int_t> cycle2_1; cycle2_1.push_back(99066);  cycle2_1.push_back(156817); run2.push_back(cycle2_1);
  vector<Int_t> cycle2_2; cycle2_2.push_back(194442); cycle2_2.push_back(266057); run2.push_back(cycle2_2);
  vector<Int_t> cycle2_3; cycle2_3.push_back(248802); cycle2_3.push_back(294633); run2.push_back(cycle2_3);
  vector<Int_t> cycle2_4; cycle2_4.push_back(287602); cycle2_4.push_back(320289); run2.push_back(cycle2_4);
  vector<Int_t> cycle2_5; cycle2_5.push_back(297130); cycle2_5.push_back(366673); run2.push_back(cycle2_5);
  vector<Int_t> cycle2_6; cycle2_6.push_back(357530); cycle2_6.push_back(447113); run2.push_back(cycle2_6);
  vector<Int_t> cycle2_7; cycle2_7.push_back(403650); cycle2_7.push_back(500729); run2.push_back(cycle2_7);
  vector<Int_t> cycle2_8; cycle2_8.push_back(694394); cycle2_8.push_back(775201); run2.push_back(cycle2_8);

  vector<vector<Int_t>> run3;
  vector<Int_t> cycle3_1; cycle3_1.push_back(62370);  cycle3_1.push_back(162921); run3.push_back(cycle3_1);
  vector<Int_t> cycle3_2; cycle3_2.push_back(200058); cycle3_2.push_back(299993); run3.push_back(cycle3_2);
  vector<Int_t> cycle3_3; cycle3_3.push_back(264074); cycle3_3.push_back(354457); run3.push_back(cycle3_3);
  vector<Int_t> cycle3_4; cycle3_4.push_back(337330); cycle3_4.push_back(418561); run3.push_back(cycle3_4);
  vector<Int_t> cycle3_5; cycle3_5.push_back(390970); cycle3_5.push_back(472673); run3.push_back(cycle3_5);
  vector<Int_t> cycle3_6; cycle3_6.push_back(523938); cycle3_6.push_back(667153); run3.push_back(cycle3_6); 
  vector<Int_t> cycle3_7; cycle3_7.push_back(596978); cycle3_7.push_back(720817); run3.push_back(cycle3_7);
  vector<Int_t> cycle3_8; cycle3_8.push_back(881058); cycle3_8.push_back(926377); run3.push_back(cycle3_8);

  vector<vector<Int_t>> run4;
  vector<Int_t> cycle4_1; cycle4_1.push_back(5752);   cycle4_1.push_back(58495);  run4.push_back(cycle4_1);
  vector<Int_t> cycle4_2; cycle4_2.push_back(175008); cycle4_2.push_back(230319); run4.push_back(cycle4_2);
  vector<Int_t> cycle4_3; cycle4_3.push_back(212952); cycle4_3.push_back(278543); run4.push_back(cycle4_3);
  vector<Int_t> cycle4_4; cycle4_4.push_back(392448); cycle4_4.push_back(447407); run4.push_back(cycle4_4);
  vector<Int_t> cycle4_5; cycle4_5.push_back(430408); cycle4_5.push_back(494007); run4.push_back(cycle4_5);
  vector<Int_t> cycle4_6; cycle4_6.push_back(477112); cycle4_6.push_back(538567); run4.push_back(cycle4_6);
  vector<Int_t> cycle4_7; cycle4_7.push_back(523832); cycle4_7.push_back(582775); run4.push_back(cycle4_7);
  vector<Int_t> cycle4_8; cycle4_8.push_back(671768); cycle4_8.push_back(726255); run4.push_back(cycle4_8);

  vector<vector<Int_t>> run5;
  vector<Int_t> cycle5_1; cycle5_1.push_back(23702);  cycle5_1.push_back(107661); run5.push_back(cycle5_1);
  vector<Int_t> cycle5_2; cycle5_2.push_back(137894); cycle5_2.push_back(220981); run5.push_back(cycle5_2);
  vector<Int_t> cycle5_3; cycle5_3.push_back(361870); cycle5_3.push_back(426565); run5.push_back(cycle5_3);
  vector<Int_t> cycle5_4; cycle5_4.push_back(417262); cycle5_4.push_back(466693); run5.push_back(cycle5_4);
  vector<Int_t> cycle5_5; cycle5_5.push_back(496534); cycle5_5.push_back(545013); run5.push_back(cycle5_5);
  vector<Int_t> cycle5_6; cycle5_6.push_back(574854); cycle5_6.push_back(638573); run5.push_back(cycle5_6);
  vector<Int_t> cycle5_7; cycle5_7.push_back(668150); cycle5_7.push_back(722557); run5.push_back(cycle5_7);
  vector<Int_t> cycle5_8; cycle5_8.push_back(847270); cycle5_8.push_back(898509); run5.push_back(cycle5_8);

  makeRunPedestalPlots(runNum1, run1);
  makeRunPedestalPlots(runNum2, run2);
  makeRunPedestalPlots(runNum3, run3);
  makeRunPedestalPlots(runNum4, run4);
  makeRunPedestalPlots(runNum5, run5);
}
