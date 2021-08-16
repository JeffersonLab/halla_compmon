#include "../online/utils.h"
#include <vector>

using namespace std;

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
  int pedSamp = 30; double pre_ped = 0; double post_ped = 0;
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
    sum += (pre_ped)*1.0/(pedSamp) - snapshot[i];
    if(snapshot[i] < minY){minY = snapshot[i]; minX = i;}
    if(snapshot[i] > maxY){maxY = snapshot[i];}
  }
  vector<int> statsVec; statsVec.push_back(minY); statsVec.push_back(maxY);
  statsVec.push_back(minX); statsVec.push_back(sum); statsVec.push_back(min_ped);
  statsVec.push_back(max_ped); statsVec.push_back(pre_ped*1000.0/pedSamp);
  statsVec.push_back(post_ped*1000.0/pedSamp);
  return statsVec;
}

void snapshotSumPlots(TChain* snapshots, Int_t runNum){
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

  TH2F *hNeg = new TH2F("hNeg", Form("Negative Sums in Run %i", runNum), 200, 0, 1.6e6, 200, -800e3, 0);

  int accepted = 0;
  for(int i = 0; i < snapshots->GetEntries(); i++){
    snapshots->GetEntry(i);
    vector<int> stats; stats = get_stats(snapshot, numSamples);
    int minY = stats[0]; int maxY = stats[1]; int minX = stats[2]; int sum = stats[3];
    int min_ped = stats[4]; int max_ped = stats[5]; double pre_ped = stats[6]*1.0/1000.0; double post_ped = stats[7]*1.0/1000.0;
    vector<bool> cuts = accept(snapshot, numSamples, minY, maxY, minX, min_ped, max_ped, pre_ped, post_ped, snapClock);
    bool allCuts = true;
    if(randomTime == 1){randoms++; continue;}
    if(snapClock > 1.6e6 || sum > 0e3){continue;}
    for(int j = 0; j < 5; j++){
      //if(j >= 0){continue;}
      allCuts = allCuts && cuts[j];
      if(not cuts[j]){rej_cuts[j]++;}
    }
    hNeg->Fill(snapClock, sum);
  }

  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", "Negative Pulse Sums Canvas", 1200, 800);
  c->cd();
  hNeg->Draw("colz");
}

void snapshotSumVerify(Int_t runNum){
  vector<TChain *> chains = loadChain(runNum);
  TChain *snapshots = chains[6];
  snapshotSumPlots(snapshots, runNum);
}
