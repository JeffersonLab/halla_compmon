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

void makeRunPedestalPlots(Int_t runNum){
  printf("Plotting pedestals for run %i...\n", runNum);
  TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *snapshots = (TTree *)f->Get("snapshots");

   float snapshot[2000];
  int randomTime, mpsCoda, numSamples, snapClock, beamState, laserState;
  int randoms = 0;
 TH1F *hON = new TH1F(Form("hON_run%i", runNum), Form("Run %i Pedestal Mean", runNum), 50, 3770, 3800);
  hON->SetLineColor(kGreen + 2); hON->SetStats(1); 
  hON->GetXaxis()->SetTitle("Pedestal Mean [RAU]");
  TH1F *hOFF = new TH1F(Form("hOFF_run%i", runNum), Form("Run %i Pedestal Mean", runNum), 50, 3770, 3800);
  hOFF->SetLineColor(kRed); hOFF->SetStats(1);
  hOFF->GetXaxis()->SetTitle("Pedestal Mean [RAU]"); 

  TH2F *h = new TH2F("h","h", 100,-0.001,0.001,100,-0.001,0.001); 
  snapshots->SetBranchAddress("randomTime", &randomTime);
  //snapshots->SetBranchAddress("snap", &snapshot);
  snapshots->SetBranchAddress("mpsCoda", &mpsCoda);
  snapshots->SetBranchAddress("numSamples", &numSamples);
  snapshots->SetBranchAddress("snapClock", &snapClock);
  snapshots->SetBranchAddress("beamState", &beamState);
  snapshots->SetBranchAddress("laserState", &laserState);  
  // snapshots->SetBranchAddress("snap",&(snap[0]));

  float_t snap[2000];
  snapshots->SetBranchAddress("snap",&(snap[0]));
  // for(Int_t i = 0; i < snapshots->GetEntries(); i++){
  snapshots->GetEntry(10);   

  std::vector<double> avgs;
  double all_avg = 0;
  double avg;
  for( int s = 0; s < 300; s++) {
   if( beamState == 1 && randomTime == 1 && snapClock > 1.6e6 && snapClock < 100 );
    all_avg += snap[s];
    avg = 0;
    bool did_sum=false;
    for(int is = 0; is < 10 && s+10 < 300; is++) {
        if( beamState == 1 && randomTime == 1 && snapClock>1.6e6 && snapClock< 100);
      avg+= snap[s+is];
      did_sum = true;
    }
    if(did_sum) {
      avgs.push_back(avg/10.);
    }
   
  }
  
  all_avg/=300;
  cout << all_avg << endl;

 double  rolling_avg = 0;
  for (int i = 0; i < avgs.size(); i++) {
    rolling_avg += avgs[i];
    
  }
   rolling_avg/=289;
   //   h->Fill(all_avg, all_avg-rolling_avg);

   double sum = 0;
for(Int_t i = 0; i < snapshots->GetEntries(); i++){
  snapshots->GetEntry(i);
  if(all_avg-rolling_avg>2) continue;
  sum += snapshot[i];
  if((laserState == 0 || laserState == 1) && randomTime == 1) hON->Fill(sum);     else if((laserState == 2 || laserState == 3) && randomTime == 1) hOFF->Fill(sum);
 }

   // cout << rolling_avg << endl;

  // for (int i = 0; i < avgs.size(); i++) {
  // std::cout  << avgs[i] << " ";
  // }
  // std::cout << std::endl;
  //}
  // all_avg will be average of full snapshot
  // avgs will be vector of rolling averages (size of 289)

   //  std::vector<double> good_avgs;
//      for( auto tmp_avg: avgs ) {  // Loop over all entries in avgs, and store each entry in tmp_avg
//        if ((all_avg - tmp_avg) > 0.001)   // Setup criteria for "good average"
//      good_avgs.push_back(tmp_avg); // Store in vector with only "good averages"
//    }

  
//      for(auto value : good_avgs)
//      if((laserState == 0 || laserState == 1) && randomTime == 1) hON->Fill(value);
//      else if((laserState == 2 || laserState == 3) && randomTime == 1) hOFF->Fill(value);
   }




//  gStyle->SetOptStat();
//   gStyle->SetStatFormat("6.6g");


   TCanvas *c = new TCanvas(Form("cRun%i", runNum), Form("Run %i Pedestals", runNum), 700, 700);
  //   c->Divide(1,2);
  c->SetLogy();
   c->cd();
  h->Draw();
  // c->cd(2);
  // hOFF->Draw();

   TCanvas *c1 = new TCanvas(Form("c1Run%i", runNum), Form("Run %i Pedestals", runNum), 700, 700);
// c->Divide(1,2);
  c1->SetLogy();
  c1->cd();
  // hON->Draw();
  // c->cd(2);
  hOFF->Draw();
 
 

//   // TPaveText *ptON = new TPaveText(0.78, 0.82, 0.9, 0.9, "blNDC");
//   // ptON->SetFillColor(0); ptON->SetBorderSize(1);
//   // TPaveText *ptOFF = new TPaveText(0.78, 0.74, 0.9, 0.82, "blNDC");
//   // ptOFF->SetFillColor(0); ptOFF->SetBorderSize(1);
//   // ptON->AddText("----Laser ON----")->SetTextColor(kGreen + 2);
//   // ptON->AddText(Form("Mean: %.2f  %.2f +/- %.2f", hON->GetEntryNumber(),hON->GetMean(), hON->GetMeanError()))->SetTextColor(kGreen + 2);
//   // ptOFF->AddText("----Laser OFF----")->SetTextColor(kRed);
//   // ptOFF->AddText(Form("Mean: %.2f +/- %.2f", hOFF->GetMean(), hOFF->GetMeanError()))->SetTextColor(kRed);
//   // ptON->Draw("same");
//   // ptOFF->Draw("same");

 
}

void snapshotPedestal_amali(){
  Int_t runNum1 = 4364; //Specifically chosen for stability
  //  Int_t runNum2 = 4614;
  // Int_t runNum3 = 4554;
  //Int_t runNum4 = 4430;
  //Int_t runNum5 = 4345;

  
  makeRunPedestalPlots(runNum1);
  // makeRunPedestalPlots(runNum2);
  // makeRunPedestalPlots(runNum3);
  //makeRunPedestalPlots(runNum4);
  //makeRunPedestalPlots(runNum5);
}
