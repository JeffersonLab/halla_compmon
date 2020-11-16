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

TChain *T;
TChain *snapshots;

void makeRunPedestalPlots(){

    T = new  TChain("ComptonG4");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4345.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4346.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4353.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4361.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4364.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4436.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4592.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4601.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4615.root");
    T->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4511.root");
    
    snapshots = new TChain("snapshots");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4345.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4346.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4353.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4361.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4364.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4436.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4592.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4601.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4615.root");
    snapshots->Add("/compton/data/cmuwork/rootfiles/prex/compmon_4511.root");


    // TTree *snapshots = (TTree *)chain->Get("snapshots");

  //  float snapshot[];
  int randomTime, mpsCoda, numSamples, snapClock, beamState, laserState;
  int randoms = 0;
  TH1F *hON = new TH1F("hON", "hON", 400,3700, 3900);
  hON->SetLineColor(kGreen + 2); hON->SetStats(1); 
  hON->GetXaxis()->SetTitle("ped [RAU]");
  TH1F *hOFF = new TH1F("hOFF", "hOFF", 400, 3700, 3900);
  hOFF->SetLineColor(kRed); hOFF->SetStats(1);
  hOFF->GetXaxis()->SetTitle("ped [RAU]");  
  TH1F *hON1 = new TH1F("hON1", "hON1", 400,0, 200);
  hON1->SetLineColor(kGreen + 2); hON1->SetStats(1); 
  hON1->GetXaxis()->SetTitle("max_deviation [RAU]");
  TH1F *hOFF1 = new TH1F("hOFF1", "hOFF1", 400, 0,200);
  hOFF1->SetLineColor(kRed); hOFF1->SetStats(1);
  hOFF1->GetXaxis()->SetTitle("max_deviation [RAU]");

  snapshots->SetBranchAddress("randomTime", &randomTime);
  //snapshots->SetBranchAddress("snap", &snapshot);
  snapshots->SetBranchAddress("mpsCoda", &mpsCoda);
  snapshots->SetBranchAddress("numSamples", &numSamples);
  snapshots->SetBranchAddress("snapClock", &snapClock);
  snapshots->SetBranchAddress("beamState", &beamState);
  snapshots->SetBranchAddress("laserState", &laserState);  
  // snapshots->SetBranchAddress("snap",&(snap[0]));

  float_t snap[2000];
  double sum = 0;
  snapshots->SetBranchAddress("snap",&(snap[0]));
   
  for(Int_t i = 0; i < snapshots->GetEntries(); i++){
  snapshots->GetEntry(i);   
  std::vector<double> avgs;
  double all_avg = 0;
  double avg;
  
   for( int s = 0; s < 300; s++) {
     if( beamState == 1){
     all_avg += snap[s];
       }
     avg = 0;
     bool did_sum=false;

   for(int is = 0; is < 10 && s+10 < 300; is++) {
     if(beamState==1){
     avg+= snap[s+is];
      }
      did_sum = true;
      
    }
    if(did_sum) {
      avgs.push_back(avg/10.);
    }
   
   }
  
  all_avg/=300;
  // cout << all_avg << endl;

  double max = *max_element(avgs.begin(), avgs.end());
  //  cout << max << endl;
   

  double max_dev =0;
   max_dev = max - all_avg;
   //  cout << max_dev <<endl;
  
      if((laserState == 0 || laserState == 1) && randomTime == 1 && beamState ==1){
 hON1->Fill(max_dev);
 }
 else if((laserState == 2 || laserState == 3) && randomTime == 1 && beamState ==1){
 hOFF1->Fill(max_dev);
 }
 
 if(max_dev<30) { // Include the sum of this snapshot in the global sum
    sum = all_avg;
  

 if((laserState == 0 || laserState == 1) && randomTime == 1 && beamState ==1){
 hON->Fill(sum);
 }
 else if((laserState == 2 || laserState == 3) && randomTime == 1 && beamState ==1){
 hOFF->Fill(sum);
 }
  }
  }

 gStyle->SetOptStat();
  gStyle->SetStatFormat("6.6g");


  TCanvas *c = new TCanvas("c", "c", 700, 700);
  //   c->Divide(1,2);
  c->SetLogy();
  c->cd();
  hON->Draw();
  // c->cd(2);
  // hOFF->Draw();

   TCanvas *c1 = new TCanvas("c1", "c1", 700, 700);
// c->Divide(1,2);
    c1->SetLogy();
  c1->cd();
  // hON->Draw();
  // c->cd(2);
  hOFF->Draw();

 TCanvas *c2 = new TCanvas("c2", "c2", 700, 700);
  //   c->Divide(1,2);
 //  c2->SetLogy();
  c2->cd();
  hON1->Draw();
  // c->cd(2);
  // hOFF->Draw();

   TCanvas *c3 = new TCanvas("c3", "c3", 700, 700);
// c->Divide(1,2);
   // c3->SetLogy();
  c3->cd();
  // hON->Draw();
  // c->cd(2);
  hOFF1->Draw();
 
 

  // TPaveText *ptON = new TPaveText(0.78, 0.82, 0.9, 0.9, "blNDC");
  // ptON->SetFillColor(0); ptON->SetBorderSize(1);
  // TPaveText *ptOFF = new TPaveText(0.78, 0.74, 0.9, 0.82, "blNDC");
  // ptOFF->SetFillColor(0); ptOFF->SetBorderSize(1);
  // ptON->AddText("----Laser ON----")->SetTextColor(kGreen + 2);
  // ptON->AddText(Form("Mean: %.2f  %.2f +/- %.2f", hON->GetEntryNumber(),hON->GetMean(), hON->GetMeanError()))->SetTextColor(kGreen + 2);
  // ptOFF->AddText("----Laser OFF----")->SetTextColor(kRed);
  // ptOFF->AddText(Form("Mean: %.2f +/- %.2f", hOFF->GetMean(), hOFF->GetMeanError()))->SetTextColor(kRed);
  // ptON->Draw("same");
  // ptOFF->Draw("same");

 
}

void snapshotPedestal_test1(){
  //  Int_t runNum1 = 4364; //Specifically chosen for stability
  //  Int_t runNum2 = 4614;
  // Int_t runNum3 = 4554;
  //Int_t runNum4 = 4430;
  //Int_t runNum5 = 4345;

  
   makeRunPedestalPlots();
  // makeRunPedestalPlots(runNum2);
  // makeRunPedestalPlots(runNum3);
  //makeRunPedestalPlots(runNum4);
  //makeRunPedestalPlots(runNum5);
}
