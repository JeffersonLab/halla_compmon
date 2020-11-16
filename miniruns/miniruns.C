#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <fstream>
#include <string>

using namespace std;

int translate_to_index(TString conf){
  TString on("L_ON"); TString off("L_OFF"); TString beamOff("B_OFF");
  if(conf.EqualTo(on)) return 0;
  else if(conf.EqualTo(off)) return 1;
  else if (conf.EqualTo(beamOff)) return 2;
  else return -1;
}

TString translate_to_string(int index){
  TString on("LASER ON"); TString off("LASER OFF"); TString beamOff("BEAM OFF");
  TString states[3] = {on, off, beamOff};
  return states[index];
}

vector<vector<vector<int>>> identify_run_states(TTree* quartetwise){
  vector<vector<vector<int>>> laserTimes;
  vector<vector<int>> laserONtimes; vector<vector<int>> laserOFFtimes; vector<vector<int>> beamOFFtimes;
  laserTimes.push_back(laserONtimes); laserTimes.push_back(laserOFFtimes); laserTimes.push_back(beamOFFtimes);  

  Int_t laserState, beamState, firstMPSnumber;
  Double_t posAcc, negAcc, posSamp, negSamp;
  //Float_t bcm;
  TString curState(""); TString prevState(""); 
  TString on("L_ON"); TString off("L_OFF"); TString beamOff("B_OFF");

  //quartetwise->SetBranchAddress("bcm", &bcm);
  quartetwise->SetBranchAddress("laserState", &laserState); quartetwise->SetBranchAddress("beamState", &beamState);
  quartetwise->SetBranchAddress("PosHelAcc0", &posAcc); quartetwise->SetBranchAddress("PosHelNSamples0", &posSamp);
  quartetwise->SetBranchAddress("NegHelAcc0", &negAcc); quartetwise->SetBranchAddress("NegHelNSamples0", &negSamp);
  quartetwise->SetBranchAddress("firstMPSnumber", &firstMPSnumber);
  Int_t max_mps = 0;

  for(int i = 0; i < quartetwise->GetEntries(); i++){
    quartetwise->GetEntry(i);
    TString thisState("");
    if(beamState==1){
      if(laserState==0 || laserState==1){thisState = on;}
      else{thisState = off;}
    }
    else{thisState = beamOff;}
    if(not curState.EqualTo(thisState)){
      if(translate_to_index(curState) >= 0){
        laserTimes[translate_to_index(curState)].back().push_back(firstMPSnumber - 1);
      }
      vector<int> period;
      laserTimes[translate_to_index(thisState)].push_back(period);
      laserTimes[translate_to_index(thisState)].back().push_back(firstMPSnumber);
      prevState = curState; curState = thisState;
    }
    if(firstMPSnumber > max_mps){max_mps = firstMPSnumber;}
  }
  
  laserTimes[translate_to_index(curState)].back().push_back(max_mps + 7);
  return laserTimes;
}

vector<vector<int>> refine_run_states(vector<vector<vector<int>>> run_states){
  vector<vector<int>> new_run_states;
  int limits[3] = {0, 0, 0};
  while(limits[0] < run_states[0].size() || limits[1] < run_states[1].size() || limits[2] < run_states[2].size()){
    int min_index = -1;
    if(limits[0] < run_states[0].size() && limits[1] < run_states[1].size() && limits[2] < run_states[2].size()){
      if(     run_states[0][limits[0]][0] < run_states[1][limits[1]][0] && run_states[0][limits[0]][0] < run_states[2][limits[2]][0]){min_index = 0;}
      else if(run_states[1][limits[1]][0] < run_states[0][limits[0]][0] && run_states[1][limits[1]][0] < run_states[2][limits[2]][0]){min_index = 1;}
      else if(run_states[2][limits[2]][0] < run_states[0][limits[0]][0] && run_states[2][limits[2]][0] < run_states[1][limits[1]][0]){min_index = 2;}
    }
    else if(limits[0] < run_states[0].size() && limits[1] < run_states[1].size()){
      if(     run_states[0][limits[0]][0] < run_states[1][limits[1]][0]){min_index = 0;}
      else if(run_states[1][limits[1]][0] < run_states[0][limits[0]][0]){min_index = 1;}
    }
    else if(limits[0] < run_states[0].size() && limits[2] < run_states[2].size()){
      if(     run_states[0][limits[0]][0] < run_states[2][limits[2]][0]){min_index = 0;}
      else if(run_states[2][limits[2]][0] < run_states[0][limits[0]][0]){min_index = 2;}
    }
    else if(limits[1] < run_states[1].size() && limits[2] < run_states[2].size()){
      if(     run_states[1][limits[1]][0] < run_states[2][limits[2]][0]){min_index = 1;}
      else if(run_states[2][limits[2]][0] < run_states[1][limits[1]][0]){min_index = 2;}
    }
    else if(limits[0] < run_states[0].size()){min_index = 0;}
    else if(limits[1] < run_states[1].size()){min_index = 1;}
    else if(limits[2] < run_states[2].size()){min_index = 2;}
    vector<int> period; period.push_back(min_index); period.push_back(run_states[min_index][limits[min_index]][0]); 
    period.push_back(run_states[min_index][limits[min_index]][1]);
    new_run_states.push_back(period);
    limits[min_index]++;
  }
  return new_run_states;
}

vector<vector<int>> identify_miniruns(vector<vector<int>> run_periods){
  vector<vector<int>> miniruns;
  int laser_on_evts = 0; int laser_off_evts = 0;
  int minirun_start = -1;
  for(int period = 0; period < run_periods.size(); period++){
    if(run_periods[period][0] == 0){laser_on_evts += run_periods[period][2] - run_periods[period][1];}
    else if(run_periods[period][0] == 1){laser_off_evts += run_periods[period][2] - run_periods[period][1];}
    if(minirun_start == -1){minirun_start = run_periods[period][1];}
    if(laser_on_evts > 6*240*60 && laser_off_evts > 3*240*60){
      vector<int> minirun; minirun.push_back(minirun_start); minirun.push_back(run_periods[period][2]);
      miniruns.push_back(minirun);
      minirun_start = -1; laser_on_evts = 0; laser_off_evts = 0;
    }
  }
  return miniruns;
}

void miniruns(TString group_fname){
  string run_num_str;
  ifstream infile(Form("%s/%s.list", getenv("COMPMON_SNAILS"), group_fname.Data()));
  while(getline(infile, run_num_str)){
    int run_num = atoi(run_num_str.c_str());
    TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), run_num));
    TTree *mpswise = (TTree *)f->Get("mpswise");
    TTree *quartetwise = (TTree *)f->Get("quartetwise");
    vector<vector<int>> miniruns = identify_miniruns(refine_run_states(identify_run_states(quartetwise)));
    f->Close();

    ofstream output;
    output.open(Form("%s/minirun_%i.dat", getenv("COMPMON_MINIRUNS"), run_num));
    for(int minirun = 0; minirun < miniruns.size(); minirun++){
      output<<miniruns[minirun][0]<<","<<miniruns[minirun][1]<<endl;
    }
    output.close();
    cout<<"Created file "<<Form("%s/minirun_%i.dat", getenv("COMPMON_MINIRUNS"), run_num)<<endl;
  }
  infile.close();
}
