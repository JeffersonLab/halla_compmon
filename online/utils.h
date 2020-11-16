#ifndef UTILS_H
#define UTILS_H

#include <TCanvas.h>
#include <TPad.h>
#include <TObject.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TPad.h>
#include <TString.h>
#include <THStack.h>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;

void resetChain(TChain *ch){
  if(ch){
    delete ch;
  }
}

vector<TChain *> loadChain(Int_t runnum){
  TChain *mpswise = 0; TChain *quartetwise = 0;
  TChain *pulserwise = 0; TChain *triggerwise = 0;
  TChain *epicswise = 0; TChain *runwise = 0;
  TChain *snapshots = 0;

  resetChain(mpswise); resetChain(quartetwise);
  resetChain(pulserwise); resetChain(triggerwise);
  resetChain(epicswise); resetChain(runwise);
  resetChain(snapshots);

  mpswise = new TChain("mpswise");
  quartetwise = new TChain("quartetwise");
  pulserwise = new TChain("pulserwise");
  triggerwise = new TChain("triggerwise");
  epicswise = new TChain("epicswise");
  runwise = new TChain("runwise");
  snapshots = new TChain("snapshots");
  std::vector<TChain*> chains;
  chains.push_back(mpswise);
  chains.push_back(quartetwise);
  chains.push_back(pulserwise);
  chains.push_back(triggerwise);
  chains.push_back(epicswise);
  chains.push_back(runwise);
  chains.push_back(snapshots);

  TString filesPre = Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),runnum);
  int nfiles = 0;
  for(size_t ch = 0; ch < chains.size(); ch++) {
    nfiles = chains[ch]->Add(filesPre+".root");
    nfiles += chains[ch]->Add(filesPre+"_*.root");

    if(nfiles<=0) {
	    std::cerr << "Looked for files under: " << filesPre+".root" << std::endl;
	    std::cerr << "Found no files to plot!" << std::endl;
	    return chains;
    }
  }
  return chains;
}

Double_t get_analyzing_power(int run_num){
  if(run_num < 3800){return 0.1105;}
  else if(run_num >= 4232 && run_num <= 4929){return 0.01655915;}
  else if(run_num >= 4930){return 0.036017934;}
  else{
    printf("Run doesn't have a defined analyzing power.\n");
    return 0.0;
  }
}

Double_t get_analyzing_power_error(int run_num){
  return 0.0;
}

TString translate(int beam, int las){
  TString beam_str("ON");
  TString laser_str("ON");
  if(beam == 0){beam_str = "OFF";}
  if(las == 0){laser_str = "OFF";}
  TString title = Form("beam %s, laser %s", beam_str.Data(), laser_str.Data());
  return title;
}

std::vector<TString> hist_stats(TH1F* h){
  std::vector<TString> stats;
  
  stats.push_back(Form("Entries: %i", (Int_t)h->GetEntries()));
  stats.push_back(Form("Mean: %.4f +/- %.4f", h->GetMean(), h->GetMeanError()));
  stats.push_back(Form("Std Dev: %.4f +/- %.4f", h->GetRMS(), h->GetRMSError()));

  return stats;
}

TString date_time_string(TFile *infile, const char* suffix=""){
  Int_t day      = (Int_t)((TH1F *)infile->Get(Form("h_day%s", suffix)))->GetMean();
  Int_t month    = (Int_t)((TH1F *)infile->Get(Form("h_month%s", suffix)))->GetMean();
  Int_t date_num = (Int_t)((TH1F *)infile->Get(Form("h_date_num%s", suffix)))->GetMean();
  Int_t year     = (Int_t)((TH1F *)infile->Get(Form("h_year%s", suffix)))->GetMean();
  Int_t hour     = (Int_t)((TH1F *)infile->Get(Form("h_hour%s", suffix)))->GetMean();
  Int_t minute   = (Int_t)((TH1F *)infile->Get(Form("h_minute%s", suffix)))->GetMean();
  Int_t second   = (Int_t)((TH1F *)infile->Get(Form("h_second%s", suffix)))->GetMean();

  TString month_str(Form("%i", month)); TString date_num_str(Form("%i", date_num)); TString hour_str(Form("%i", hour));
  TString minute_str(Form("%i", minute)); TString second_str(Form("%i", second));
  if(month < 10){month_str = Form("0%i", month);}
  if(date_num < 10){date_num_str = Form("0%i", date_num);}
  if(hour < 10){hour_str = Form("0%i", hour);}
  if(minute < 10){minute_str = Form("0%i", minute);}
  if(second < 10){second_str = Form("0%i", second);}
  return Form("%i-%s-%s %s:%s:%s", year, month_str.Data(), date_num_str.Data(), hour_str.Data(), minute_str.Data(), second_str.Data());
}

TString get_current_target(TFile *infile){
  Float_t avgTargPos = (Float_t)((TH1F *)infile->Get("epics_targetPos"))->GetMean();
  Int_t targetPos[17] = {0, 13163050, 16224300, 20324764, 23394068, 26463372, 29532678, 32601984, 35671288,
                            38751872, 41832456, 44913040, 47993624, 51074208, 54147540, 57226912, 60306288};
  TString targetNames[17] = {"HOME", "Ca48", "Ca40", "Carbon Hole", "D9-208Pb10-D10", "D7-208Pb9-D8", "D5-208Pb8-D6", "D3-208Pb7-D4", "D1-208Pb6-D2",
                             "DG-208Pb5-D20", "DE-208Pb4-DF", "DC-208Pb3-DD", "DA-208Pb2-DB", "Carbon 1%", "C-208Pb-C", "DI-Pb-DJ", "C-Pb-C"};
  for(int i = 0; i < 17; i++){
    if(TMath::Abs(avgTargPos - targetPos[i]) < 10000){
      return targetNames[i];
    }
  }

  return "UNK";
}

bool isCloseTo(Float_t x, Float_t y){return TMath::Abs(x - y)/x < 0.001;}

bool hasLowError(Float_t mean, Float_t err){return TMath::Abs(err/mean) < 1e-4;}

TString get_wien_state(TFile *infile, Int_t run_num){
  if(run_num < 4700 || 
    !(TH1F *)infile->Get("epics_VWienAngle") ||
    !(TH1F *)infile->Get("epics_HWienAngle") ||
    !(TH1F *)infile->Get("epics_PhiFG")){return "UNK";}
  Float_t vWienAngle = ((TH1F *)infile->Get("epics_VWienAngle"))->GetMean();
  Float_t hWienAngle = ((TH1F *)infile->Get("epics_HWienAngle"))->GetMean();
  Float_t phiFG = ((TH1F *)infile->Get("epics_PhiFG"))->GetMean();
  if(isCloseTo(vWienAngle, 88.0) && isCloseTo(hWienAngle, -29.64) && isCloseTo(phiFG, 89.956)){return "FLIP-RIGHT";}
  else if(isCloseTo(vWienAngle, -90.6) && isCloseTo(hWienAngle, -29.64) && isCloseTo(phiFG, 88.028)){return "FLIP-LEFT";}
  else{return "UNK";}
}

Int_t numValidLaserCycles(Int_t run_num){
  ifstream infile(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), run_num));
  string read_str;
  Int_t nValidCycles = 0;
  while(getline(infile, read_str)){
    nValidCycles++;
  }
  return nValidCycles;
}

std::vector<TString> calc_polarization(TFile *infile, int run_num, int accum, int beamCut=1){
  TString tree_name("quartetwise");
  if(!infile->IsOpen()){printf("Infile is bad here...\n");}
  TH1F *hPoshOn  = (TH1F *)infile->Get(Form("%s_beam%i_las1_posH%i", tree_name.Data(), beamCut, accum));
  TH1F *hNeghOn  = (TH1F *)infile->Get(Form("%s_beam%i_las1_negH%i", tree_name.Data(), beamCut, accum));
  TH1F *hDiffOn  = (TH1F *)infile->Get(Form("%s_beam%i_las1_diff%i", tree_name.Data(), beamCut, accum));
  TH1F *hSummOn  = (TH1F *)infile->Get(Form("%s_beam%i_las1_summ%i", tree_name.Data(), beamCut, accum));
  TH1F *hAsymOn  = (TH1F *)infile->Get(Form("%s_beam%i_las1_asym%isubOn", tree_name.Data(), beamCut, accum));
  TH1F *hPoshOff = (TH1F *)infile->Get(Form("%s_beam%i_las0_posH%i", tree_name.Data(), beamCut, accum));
  TH1F *hNeghOff = (TH1F *)infile->Get(Form("%s_beam%i_las0_negH%i", tree_name.Data(), beamCut, accum));
  TH1F *hDiffOff = (TH1F *)infile->Get(Form("%s_beam%i_las0_diff%i", tree_name.Data(), beamCut, accum));
  TH1F *hSummOff = (TH1F *)infile->Get(Form("%s_beam%i_las0_summ%i", tree_name.Data(), beamCut, accum));
  TH1F *hAsymOff = (TH1F *)infile->Get(Form("%s_beam%i_las0_asym%isubOff", tree_name.Data(), beamCut, accum));

  Double_t vPoshOn  = hPoshOn->GetMean();  Double_t vPoshOnE  = hPoshOn->GetMeanError();
  Double_t vNeghOn  = hNeghOn->GetMean();  Double_t vNeghOnE  = hNeghOn->GetMeanError();
  Double_t vDiffOn  = hDiffOn->GetMean();  Double_t vDiffOnE  = hDiffOn->GetMeanError();
  Double_t vSummOn  = hSummOn->GetMean();  Double_t vSummOnE  = hSummOn->GetMeanError();
  Double_t vAsymOnHist  = hAsymOn->GetMean();  Double_t vAsymOnHistE  = hAsymOn->GetMeanError();
  Double_t vPoshOff = hPoshOff->GetMean(); Double_t vPoshOffE = hPoshOff->GetMeanError();
  Double_t vNeghOff = hNeghOff->GetMean(); Double_t vNeghOffE = hNeghOff->GetMeanError();
  Double_t vDiffOff = hDiffOff->GetMean(); Double_t vDiffOffE = hDiffOff->GetMeanError();
  Double_t vSummOff = hSummOff->GetMean(); Double_t vSummOffE = hSummOff->GetMeanError();
  Double_t vAsymOffHist = hAsymOff->GetMean(); Double_t vAsymOffHistE = hAsymOff->GetMeanError();

  Double_t vAsymOn    = 1000*vDiffOn/(vSummOn - vSummOff); Double_t vAsymOff = 1000*vDiffOff/(vSummOn - vSummOff);
  Double_t vSummMeanE = TMath::Sqrt(TMath::Power(vSummOnE, 2) + TMath::Power(vSummOffE, 2));
  Double_t vAsymOnE   = vAsymOn*TMath::Sqrt(TMath::Power(vDiffOnE/vDiffOn, 2)   + TMath::Power(vSummMeanE/(vSummOn - vSummOff), 2));
  Double_t vAsymOffE  = vAsymOff*TMath::Sqrt(TMath::Power(vDiffOffE/vDiffOff, 2) + TMath::Power(vSummMeanE/(vSummOn - vSummOff), 2));
  Double_t vAsymSub   = vAsymOn - vAsymOff;
  Double_t vAsymSubE  = TMath::Sqrt(TMath::Power(vAsymOnE, 2) + TMath::Power(vAsymOffE, 2));
  Double_t pol = vAsymOn/get_analyzing_power(run_num); Double_t polE = vAsymOnE/get_analyzing_power(run_num);

  vector<TString> results; results.push_back("--------Histogram Vars--------\n");
  results.push_back(Form("Positive Hel (ON): %f +/- %f\n",  vPoshOn,  vPoshOnE));  results.push_back(Form("Negative Hel (ON): %f +/- %f\n",  vNeghOn,  vNeghOnE));
  results.push_back(Form("Positive Hel (OFF): %f +/- %f\n", vPoshOff, vPoshOffE)); results.push_back(Form("Negative Hel (OFF): %f +/- %f\n", vNeghOff, vNeghOffE));
  results.push_back("--------Sums and Diffs--------\n");
  results.push_back(Form("On Diff: %f +/- %f\n", vDiffOn, vDiffOnE));   results.push_back(Form("Off Diff: %f +/- %f\n", vDiffOff, vDiffOffE));
  results.push_back(Form("On Sum: %f +/- %f\n",  vSummOn,  vSummOnE));  results.push_back(Form("Off Sum: %f +/- %f\n",  vSummOff,  vSummOffE));
  results.push_back("--------Asymmetries--------\n");
  results.push_back(Form("On Asym: %f +/- %f (ppt)\n", vAsymOn, vAsymOnE)); results.push_back(Form("Off Asym: %f +/- %f (ppt)\n", vAsymOff, vAsymOffE));
  results.push_back("--------Bkgd Subtraction--------\n");
  results.push_back(Form("BkSub Sum: %f +/- %f\n", vSummOn - vSummOff, vSummMeanE));
  results.push_back(Form("BkSub Asym: %f +/- %f (ppt)\n", vAsymSub, vAsymSubE));
  //results.push_back("--------Polarization--------\n");
  //results.push_back(Form("Polarization: %f +/- %f (%f%%)", 0.1*pol, 0.1*polE, 100*polE/TMath::Abs(pol)));
  
  return results;
}

void essential_stats_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TH1F *acc0_all = (TH1F *)infile->Get("mpswise_acc0");
  TH1F *acc0_bON_lON = (TH1F *)infile->Get("mpswise_beam1_las1_acc0");
  TH1F *acc0_bON_lOFF = (TH1F *)infile->Get("mpswise_beam1_las0_acc0");
  TH1F *acc0_bOFF = (TH1F *)infile->Get("mpswise_beam0_acc0");
  TH1F *h_ihwp = (TH1F *)infile->Get("mpswise_ihwp");
  TH1F *h_qw1 = (TH1F *)infile->Get("epics_qw1");
  TH1F *h_hw1 = (TH1F *)infile->Get("epics_hw1");
  TString ihwp_state("UNK");
  if(h_ihwp->GetMean() < 0.05) ihwp_state = "OUT";
  else if(h_ihwp->GetMean() > 0.95) ihwp_state = "IN";
  float flip_freq = 240.0;
  if(run_num > 4700){flip_freq = 120.0;}

  Double_t seconds = acc0_all->GetEntries()*1.0/(flip_freq); Int_t iSeconds = (Int_t)seconds;
  Double_t minutes = seconds/60.0; Int_t iMinutes = (Int_t)minutes;
  Double_t hours = minutes/60.0; Int_t iHours = (Int_t)hours;
  Int_t realMinutes = (Int_t)((hours - iHours)*60.0);
  Int_t realSeconds = (Int_t)((minutes - iMinutes)*60.0);

  Int_t lON_evts = acc0_bON_lON->GetEntries();
  Int_t lOFF_evts = acc0_bON_lOFF->GetEntries();
  Int_t bOFF_evts = acc0_bOFF->GetEntries();
  Int_t tot_evts = lON_evts + lOFF_evts + bOFF_evts;

  Float_t qw1_deg = fmod((-0.004*h_qw1->GetMean() - 1390.768), 360);
  if(run_num > 4700){qw1_deg = fmod((-0.004*h_qw1->GetMean()), 360);}
  Float_t hw1_deg = fmod((-0.004*h_hw1->GetMean()), 360);
  
  TPaveText *essStats = new TPaveStats(0.0, 0.0, 1.0, 1.0, "blNDC");
  essStats->AddText(Form("--------Essential Statistics: Run %i--------", run_num));
  essStats->AddText(Form("First Epics Time: %s", date_time_string(infile).Data()));
  essStats->AddText(Form("Last Epics Time: %s", date_time_string(infile, "_last").Data()));
  essStats->AddText(Form("Run Time: %i hour(s), %i minute(s), %i second(s)", iHours, realMinutes, realSeconds));
  essStats->AddText(Form("Current IHWP State: %s", ihwp_state.Data()));
  essStats->AddText(Form("QW1 Pos: %.1f deg, HW1 Pos: %.1f deg", qw1_deg, hw1_deg));
  essStats->AddText(Form("# of Events: %i", (int)acc0_all->GetEntries()));
  essStats->AddText(Form("# of Beam ON, Laser ON Events: %i, (%.2f%%)", lON_evts, lON_evts*100.0/tot_evts));
  essStats->AddText(Form("# of Beam ON, Laser OFF Events: %i, (%.2f%%)", lOFF_evts, lOFF_evts*100.0/tot_evts));
  essStats->AddText(Form("# of Beam OFF Events: %i, (%.2f%%)", bOFF_evts, bOFF_evts*100.0/tot_evts));
  essStats->SetBorderSize(0);
  essStats->SetFillColor(0);
  myPad->cd(); essStats->Draw();

  if(printfile)
    myPad->Print(Form("%s/ess_stats_1.png", output_path.Data()), "png");
}

void essential_stats_2_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString target = get_current_target(infile);
  TH1F *h_ihwp = (TH1F *)infile->Get("mpswise_ihwp");
  TString ihwp_state("UNK");
  if(h_ihwp->GetMean() < 0.05) ihwp_state = "OUT";
  else if(h_ihwp->GetMean() > 0.95) ihwp_state = "IN";
  TString wien_state = get_wien_state(infile, run_num);

  Float_t Ax_mean = ((TH1F *)infile->Get("epics_bpmAx"))->GetMean(); Float_t Ax_err = ((TH1F *)infile->Get("epics_bpmAx"))->GetMeanError();
  Float_t Ay_mean = ((TH1F *)infile->Get("epics_bpmAy"))->GetMean(); Float_t Ay_err = ((TH1F *)infile->Get("epics_bpmAy"))->GetMeanError();
  Float_t Bx_mean = ((TH1F *)infile->Get("epics_bpmBx"))->GetMean(); Float_t Bx_err = ((TH1F *)infile->Get("epics_bpmBx"))->GetMeanError();
  Float_t By_mean = ((TH1F *)infile->Get("epics_bpmBy"))->GetMean(); Float_t By_err = ((TH1F *)infile->Get("epics_bpmBy"))->GetMeanError();
  Float_t tableX = ((TH1F *)infile->Get("epics_tablePosX"))->GetMean(); 
  Float_t tableX_err = ((TH1F *)infile->Get("epics_tablePosX"))->GetMeanError();
  Float_t tableY = ((TH1F *)infile->Get("epics_tablePosY"))->GetMean();
  Float_t tableY_err = ((TH1F *)infile->Get("epics_tablePosY"))->GetMeanError();
  Int_t qw1 = (Int_t)((TH1F *)infile->Get("epics_qw1"))->GetMean();
  Int_t hw1 = (Int_t)((TH1F *)infile->Get("epics_hw1"))->GetMean();
  Float_t qw1_deg = fmod((-0.004*qw1 - 1390.768), 360);
  if(run_num > 4700){qw1_deg = fmod((-0.004*qw1), 360);}
  Float_t hw1_deg = fmod((-0.004*hw1), 360);
  Int_t loThresh = (Int_t)((TH1F *)infile->Get("runwise_FADC_ithrnear"))->GetMean();
  Int_t hiThresh = (Int_t)((TH1F *)infile->Get("runwise_FADC_ithrfar"))->GetMean();

  TPaveText *essStats = new TPaveStats(0.0, 0.0, 1.0, 1.0, "blNDC");
  essStats->AddText(Form("--------Essential Statistics (EPICS): Run %i--------", run_num));
  essStats->AddText(Form("IHWP State: %s    Wien State: %s", ihwp_state.Data(), wien_state.Data()));
  essStats->AddText(Form("Target for this run: %s", target.Data()));
  essStats->AddText(Form("# Valid Laser Cycles Found: %i", numValidLaserCycles(run_num)));
  essStats->AddText(Form("BPM2Ax: %.3f +/- %.3f    BPM2Ay: %.3f +/- %.3f", Ax_mean, Ax_err, Ay_mean, Ay_err));
  essStats->AddText(Form("BPM2Bx: %.3f +/- %.3f    BPM2By: %.3f +/- %.3f", Bx_mean, Bx_err, By_mean, By_err));
  essStats->AddText(Form("Table X: %.2f mm    Table Y: %.2f mm", tableX, tableY));
  essStats->AddText(Form("QW1 Pos RB: %i    HW1 Pos RB: %i", qw1, hw1));
  essStats->AddText(Form("QW1 Angle: %.1f    HW1 Angle: %.1f", qw1_deg, hw1_deg));
  essStats->AddText(Form("FADC Low Thresh: %i    FADC Hi Thresh: %i", loThresh, hiThresh));
  essStats->SetBorderSize(0);
  essStats->SetFillColor(0);
  myPad->cd(); essStats->Draw();

  if(printfile)
    myPad->Print(Form("%s/ess_stats_2.png", output_path.Data()), "png");
}

void snapshots_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString tree_name("snapshots");
  TString data_name("sumPeakHeight"); TString ped_name("pedestals");
  TString acc_name("accepted");
  TString rej_names[5] = {"rej_peak", "rej_sat", "rej_corrPed", "rej_narrPed", "rej_time"};
  TString rej_titles[5] = {"Peak Height Cut", "Saturation Cut", "Correct Pedestal Cut", "Narrow Pedestal Cut", "Time Cut"};
  gStyle->SetOptStat(0);
  //myPad->Divide(1, 2);
  //myPad->cd(1); 
  myPad->cd();
  TGraph *g = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), data_name.Data()));
  TGraph *g_acc = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), acc_name.Data()));
  g->SetMarkerColor(kRed);
  g->Draw("ap");

  Double_t pct = 100.0*g_acc->GetMean(2);
  TString info = Form("%.2f%% snapshots accepted", pct);
  TPaveText *pt = new TPaveText(0.18, 0.65, 0.45, 0.88, "blNDC");
  pt->AddText(info.Data());
  for(int i = 0; i < 5; i++){
    //if(i >= 0){continue;}
    TGraph *g_rej = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), rej_names[i].Data()));
    pt->AddText(Form("%.2f%% %s", g_rej->GetMean(2)*100.0, rej_titles[i].Data()));
  }
  pt->SetBorderSize(0);
  pt->SetFillColor(0);
  pt->Draw();

  //myPad->cd(2);
  //TGraph *g_ped = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), ped_name.Data()));
  //g_ped->Draw("ap");

  if(printfile)
    myPad->Print(Form("%s/snapshots_1.png", output_path.Data()), "png");
}

void snapshots_pad_2(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString tree_name("snapshots");
  TString data_name("sumPeakHeight"); TString ped_name("pedestals");
  myPad->Divide(1, 2);
  myPad->cd(1);
  TGraph *g_ped = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), ped_name.Data()));
  g_ped->Draw("ap");
  myPad->cd(2);
  TH2F *h_ped = (TH2F *)infile->Get(Form("%s_%s_histo", tree_name.Data(), ped_name.Data()));
  h_ped->Draw("colz");

  if(printfile)
    myPad->Print(Form("%s/snapshots_2.png", output_path.Data()), "png");
}

void breakdown_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, TString tree_name, TString data_name, bool printfile=false){
  myPad->Divide(2, 2);

  TH1F *h_all = (TH1F *)infile->Get(Form("%s_%s", tree_name.Data(), data_name.Data())); myPad->cd(1); h_all->Draw();
  std::vector<TString> s1 = hist_stats(h_all);
  TPaveText *pt1 = new TPaveText(0.75, 0.7, 0.95, 0.9, "blNDC"); pt1->SetBorderSize(1); pt1->SetFillColor(0);
  for(int i = 0; i < s1.size(); i++){pt1->AddText(s1[i].Data())->SetTextColor(kBlue);} pt1->Draw("same"); 
  TH1F *h_beamOff = (TH1F *)infile->Get(Form("%s_beam0_%s", tree_name.Data(), data_name.Data())); myPad->cd(2); h_beamOff->Draw();
  std::vector<TString> s2 = hist_stats(h_beamOff);
  TPaveText *pt2 = new TPaveText(0.75, 0.7, 0.95, 0.9, "blNDC"); pt2->SetBorderSize(1); pt2->SetFillColor(0);
  for(int i = 0; i < s2.size(); i++){pt2->AddText(s2[i].Data())->SetTextColor(kBlue);} pt2->Draw("same"); 

  TH1F *h_l1_h0 = (TH1F *)infile->Get(Form("%s_beam1_las1_hel0_%s", tree_name.Data(), data_name.Data()));
  TH1F *h_l1_h1 = (TH1F *)infile->Get(Form("%s_beam1_las1_hel1_%s", tree_name.Data(), data_name.Data()));
  h_l1_h0->SetLineColor(kGreen - 3); h_l1_h1->SetLineColor(kRed); myPad->cd(3);
  THStack *hs3 = new THStack(Form("%s_beam1_las1_%s_stack", tree_name.Data(), data_name.Data()), 
      Form("Run %i %s, %s: Beam ON, Laser ON", run_num, tree_name.Data(), data_name.Data()));
  hs3->Add(h_l1_h0); hs3->Add(h_l1_h1); hs3->Draw("nostack");
  TLegend *leg3 = new TLegend(0.75, 0.75, 0.95, 0.9, "", "NDC");
  leg3->AddEntry(h_l1_h0, "Helicity==0"); leg3->AddEntry(h_l1_h1, "Helicity==1"); leg3->Draw("same");
  std::vector<TString> s_l1_h0 = hist_stats(h_l1_h0); std::vector<TString> s_l1_h1 = hist_stats(h_l1_h1);
  TPaveText *pt3_0 = new TPaveText(0.75, 0.55, 0.95, 0.75, "blNDC"); pt3_0->SetBorderSize(1); pt3_0->SetFillColor(0);
  TPaveText *pt3_1 = new TPaveText(0.75, 0.35, 0.95, 0.55, "blNDC"); pt3_1->SetBorderSize(1); pt3_1->SetFillColor(0);
  for(int i = 0; i < s_l1_h0.size(); i++){pt3_0->AddText(s_l1_h0[i].Data())->SetTextColor(kGreen - 3);}
  for(int i = 0; i < s_l1_h1.size(); i++){pt3_1->AddText(s_l1_h1[i].Data())->SetTextColor(kRed);}
  pt3_0->Draw("same"); pt3_1->Draw("same");

  TH1F *h_l0_h0 = (TH1F *)infile->Get(Form("%s_beam1_las0_hel0_%s", tree_name.Data(), data_name.Data()));
  TH1F *h_l0_h1 = (TH1F *)infile->Get(Form("%s_beam1_las0_hel1_%s", tree_name.Data(), data_name.Data()));
  h_l0_h0->SetLineColor(kGreen - 3); h_l0_h1->SetLineColor(kRed); myPad->cd(4);
  THStack *hs4 = new THStack(Form("%s_beam1_las0_%s_stack", tree_name.Data(), data_name.Data()), 
      Form("Run %i %s, %s: Beam ON, Laser OFF", run_num, tree_name.Data(), data_name.Data()));
  hs4->Add(h_l0_h0); hs4->Add(h_l0_h1); hs4->Draw("nostack");
  TLegend *leg4 = new TLegend(0.75, 0.75, 0.95, 0.9, "", "NDC");
  leg4->AddEntry(h_l0_h0, "Helicity==0"); leg4->AddEntry(h_l0_h1, "Helicity==1"); leg4->Draw("same");
  std::vector<TString> s_l0_h0 = hist_stats(h_l0_h0); std::vector<TString> s_l0_h1 = hist_stats(h_l0_h1);
  TPaveText *pt4_0 = new TPaveText(0.75, 0.55, 0.95, 0.75, "blNDC"); pt4_0->SetBorderSize(1); pt4_0->SetFillColor(0);
  TPaveText *pt4_1 = new TPaveText(0.75, 0.35, 0.95, 0.55, "blNDC"); pt4_1->SetBorderSize(1); pt4_1->SetFillColor(0);
  for(int i = 0; i < s_l0_h0.size(); i++){pt4_0->AddText(s_l0_h0[i].Data())->SetTextColor(3);}
  for(int i = 0; i < s_l0_h1.size(); i++){pt4_1->AddText(s_l0_h1[i].Data())->SetTextColor(2);}
  pt4_0->Draw("same"); pt4_1->Draw("same");

  if(printfile)
    myPad->Print(Form("%s/%s.png", output_path.Data(), data_name.Data()), "png");
}

void acc0_time_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString tree_name("mpswise");
  TString data_name("acc0_time");
  myPad->Divide(2, 2);

  gStyle->SetOptStat(0);
  for(int las = 0; las < 2; las++){
    TGraph *g_pos0 = (TGraph *)infile->Get(Form("%s_beam1_las%i_hel0_%s", tree_name.Data(), las, data_name.Data()));
    TGraph *g_pos1 = (TGraph *)infile->Get(Form("%s_beam1_las%i_hel1_%s", tree_name.Data(), las, data_name.Data()));
    TGraph *g_pos0_OFF = (TGraph *)infile->Get(Form("%s_beam0_las%i_hel0_%s", tree_name.Data(), las, data_name.Data()));
    TGraph *g_pos1_OFF = (TGraph *)infile->Get(Form("%s_beam0_las%i_hel1_%s", tree_name.Data(), las, data_name.Data()));
    myPad->cd(1);
    if(las == 0){
      g_pos0->SetTitle(Form("Run %i %s, Acc0/NAcc0:mpsCoda, all events", run_num, tree_name.Data()));
      g_pos0->Draw("ap");
    }
    else{g_pos0->Draw("p");}
    g_pos1->Draw("p");
    if(las == 0){
      myPad->cd(2);
      g_pos0_OFF->SetTitle(Form("Run %i %s, Acc0/NAcc0:mpsCoda, beam OFF", run_num, tree_name.Data())); g_pos0_OFF->Draw("ap");
      g_pos1_OFF->SetTitle(Form("Run %i %s, Acc0/NAcc0:mpsCoda, %s", run_num, tree_name.Data(), translate(1, las).Data())); g_pos1_OFF->Draw("p");
    }
    else{
      myPad->cd(2);
      g_pos0_OFF->SetTitle(Form("Run %i %s, Acc0/NAcc0:mpsCoda, beam OFF", run_num, tree_name.Data())); g_pos0_OFF->Draw("p");
      g_pos1_OFF->SetTitle(Form("Run %i %s, Acc0/NAcc0:mpsCoda, %s", run_num, tree_name.Data(), translate(1, las).Data())); g_pos1_OFF->Draw("p");
    }
    myPad->cd(3 + (int)!las);
    g_pos0->SetTitle(Form("Run %i %s, Acc0/NAcc0:mpsCoda, %s", run_num, tree_name.Data(), translate(1, las).Data()));
    g_pos0->Draw("ap"); g_pos1->Draw("p");
  }
  if(printfile)
    myPad->Print(Form("%s/%s.png", output_path.Data(), data_name.Data()), "png");
}

void quartet_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, TString data_name, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 2);

  gStyle->SetOptStat();
  TH1F *h_all = (TH1F *)infile->Get(Form("%s_%s", tree_name.Data(), data_name.Data())); myPad->cd(1); h_all->Draw();
  TH1F *h_beamOff = (TH1F *)infile->Get(Form("%s_beam0_%s", tree_name.Data(), data_name.Data())); myPad->cd(2); h_beamOff->Draw();
  TH1F *h_l1 = (TH1F *)infile->Get(Form("%s_beam1_las1_%s", tree_name.Data(), data_name.Data())); myPad->cd(3); h_l1->Draw();
  TH1F *h_l0 = (TH1F *)infile->Get(Form("%s_beam1_las0_%s", tree_name.Data(), data_name.Data())); myPad->cd(4); h_l0->Draw();

  if(printfile)
    myPad->Print(Form("%s/quart_%s.png", output_path.Data(), data_name.Data()), "png");
}

void asym_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, int accum=0, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 2); gStyle->SetOptStat(0);

  TString msmts[3] = {Form("posH%i", accum), Form("negH%i", accum), Form("asym%isub", accum)};
  TString axes[3] = {Form("posH%i (sRAU)", accum), Form("negH%i (sRAU)", accum), Form("asym%isub (ppt)", accum)};
  for(int i = 0;  i < 3; i++){
    myPad->cd(i + 1);
    TString onName = ""; TString offName = "";
    if(i != 2){
      onName = Form("%s_beam1_las1_%s", tree_name.Data(), msmts[i].Data());
      offName = Form("%s_beam1_las0_%s", tree_name.Data(), msmts[i].Data());
    }
    else{
      onName = Form("%s_beam1_las1_%sOn", tree_name.Data(), msmts[i].Data());
      offName = Form("%s_beam1_las0_%sOff", tree_name.Data(), msmts[i].Data());
    }
    TH1F *hON  = (TH1F *)infile->Get(onName.Data());
    TH1F *hOFF = (TH1F *)infile->Get(offName.Data());
    hON->SetLineColor(kGreen - 3); hOFF->SetLineColor(kRed);
    hON->GetXaxis()->SetTitle(axes[i].Data()); hOFF->GetXaxis()->SetTitle(axes[i].Data());
    hON->GetYaxis()->SetRangeUser(0, 1.2*hON->GetMaximum()); hOFF->GetYaxis()->SetRangeUser(0, 1.2*hOFF->GetMaximum());
    THStack *hs = new THStack(Form("quartet_beam1_%s_stack", msmts[i].Data()), Form("Run %i %s, %s: Beam ON", run_num, tree_name.Data(), msmts[i].Data()));
    hs->Add(hON); hs->Add(hOFF); hs->Draw("nostack");
    TLegend *leg = new TLegend(0.20, 0.80, 0.42, 0.92, "", "NDC");
    leg->AddEntry(hON, "Laser ON"); leg->AddEntry(hOFF, "Laser OFF"); leg->Draw();
    std::vector<TString> sON = hist_stats(hON); std::vector<TString> sOFF = hist_stats(hOFF);
    TPaveText *ptON  = new TPaveText(0.42, 0.80, 0.64, 0.92, "blNDC"); ptON->SetBorderSize(1);  ptON->SetFillColor(0);
    TPaveText *ptOFF = new TPaveText(0.64, 0.80, 0.86, 0.92, "blNDC"); ptOFF->SetBorderSize(1); ptOFF->SetFillColor(0);
    for(int i = 0; i < sON.size();  i++){ptON->AddText(sON[i].Data())->SetTextColor(kGreen + 1);}
    for(int i = 0; i < sOFF.size(); i++){ptOFF->AddText(sOFF[i].Data())->SetTextColor(kRed);}
    ptON->Draw(); ptOFF->Draw();
    //hs->GetXaxis()->SetTitle(axes[i].Data()); myPad->Modified();
  }
  myPad->cd(4);
  vector<TString> results = calc_polarization(infile, run_num, accum);
  TPaveText *pt1 = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
  for(int i = 0; i < results.size(); i++){pt1->AddText(results[i].Data());}
  pt1->SetBorderSize(0); pt1->SetFillColor(0);
  pt1->Draw();

  if(printfile)
    myPad->Print(Form("%s/asym_acc%i_hists.png", output_path.Data(), accum), "png");
}

void asym_graph_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, int accum=0, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 2);

  gStyle->SetOptStat(0);
  TString msmts[3] = {Form("posH%i", accum), Form("negH%i", accum), Form("asym%isub", accum)};
  for(int i = 0; i < 3; i++){
    myPad->cd(i + 1);
    TString onName = ""; TString offName = "";
    if(i != 2){
      onName = Form("%s_beam1_las1_%s_time", tree_name.Data(), msmts[i].Data());
      offName = Form("%s_beam1_las0_%s_time", tree_name.Data(), msmts[i].Data());
    }
    else{
      onName = Form("%s_beam1_las1_%sOn_time", tree_name.Data(), msmts[i].Data());
      offName = Form("%s_beam1_las0_%sOff_time", tree_name.Data(), msmts[i].Data());
    }
    TGraph *gON  = (TGraph *)infile->Get(onName.Data());
    TGraph *gOFF = (TGraph *)infile->Get(Form("%s_beam1_las0_%s_time", tree_name.Data(), msmts[i].Data()));
    gON->SetMarkerColor(kGreen - 3); gOFF->SetMarkerColor(kRed);
    gON->SetTitle(Form("Run %i %s, %s vs time", run_num, tree_name.Data(), msmts[i].Data()));
    gOFF->SetTitle(Form("Run %i %s, %s vs time", run_num, tree_name.Data(), msmts[i].Data()));
    gON->Draw("ap"); gOFF->Draw("p && same");
  }

  myPad->cd(4);
  TGraph *g_b0_l0 = (TGraph *)infile->Get(Form("%s_beam0_las0_bcm_time", tree_name.Data())); TGraph *g_b0_l1 = (TGraph *)infile->Get(Form("%s_beam0_las1_bcm_time", tree_name.Data()));
  TGraph *g_b1_l0 = (TGraph *)infile->Get(Form("%s_beam1_las0_bcm_time", tree_name.Data())); TGraph *g_b1_l1 = (TGraph *)infile->Get(Form("%s_beam1_las1_bcm_time", tree_name.Data()));
  g_b0_l0->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data())); g_b0_l1->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data()));
  g_b1_l0->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data())); g_b1_l1->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data()));
  g_b0_l0->SetMarkerColor(kRed); g_b0_l1->SetMarkerColor(kRed);
  g_b1_l0->Draw("ap"); g_b1_l1->Draw("p && same");
  g_b0_l0->Draw("* && same"); g_b0_l1->Draw("* && same"); 

  if(printfile)
    myPad->Print(Form("%s/q_acc%i_graphs.png", output_path.Data(), accum), "png");
}

void beam_off_quartet_pad(TFile *infile, int run_num, TString output_path, TPad *myPad){
  TString tree_name("quartetwise");
  myPad->Divide(2, 2);

  TString msmts[4] = {"posH0", "negH0", "asym0"};
  for(int i = 0;  i < 3; i++){
    myPad->cd(i + 1);
    TH1F *h_bOFF = (TH1F *)infile->Get(Form("%s_beam0_%s", tree_name.Data(), msmts[i].Data()));
    h_bOFF->Draw();
  }
  std::vector<TString> results = calc_polarization(infile, run_num, 0, 0);
  myPad->cd(4);
  TPaveText *pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
  for(int i = 0; i < results.size(); i++){pt->AddText(results[i].Data());}
  pt->SetBorderSize(0); pt->SetFillColor(0);
  pt->Draw();
}

void detector_asyms(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 4);

  gStyle->SetOptStat();
  TString msmts[7] = {"asym0", "USbg1", "USbg2", "DSbg1", "DSbg2", "horizFing", "vertFing"};
  TString titles[7] = {"Acc0 Asym", "USbg1 Asym", "USbg2 Asym", "DSbg1 Asym", "DSbg2 Asym", "Top Finger Asym", "Side Finger Asym"};
  Float_t means[7]; Float_t errs[7];
  for(int i = 0; i < 7; i++){
    myPad->cd(i + 1);
    TH1F *h = (TH1F *)infile->Get(Form("%s_beam1_las0_%s", tree_name.Data(), msmts[i].Data()));
    h->SetTitle(Form("Run %i %s, %s: Beam ON, Laser OFF", run_num, tree_name.Data(), titles[i].Data()));
    h->GetXaxis()->SetTitle(titles[i].Data()); h->SetLineColor(kBlue);
    h->Draw();
    means[i] = h->GetMean(); errs[i] = h->GetMeanError();
  }

  myPad->cd(8);
  TPaveText *pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
  for(int i = 0; i < 7; i++){
    pt->AddText(Form("%s (ppt): %.2f +/- %.2f", titles[i].Data(), means[i]*1000.0, errs[i]*1000.0));
    if(TMath::Abs(errs[i]/means[i]) < 0.5) ((TText*)pt->GetListOfLines()->Last())->SetTextColor(kRed);
  }
  pt->SetBorderSize(0); pt->SetFillColor(0);
  pt->Draw();

  if(printfile)
    myPad->Print(Form("%s/back_asyms.png", output_path.Data()), "png");
}

void detector_rates(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 3);
  
  gStyle->SetOptStat(0);
  TString msmts[6] = {"USbg1", "USbg2", "DSbg1", "DSbg2", "vertFing", "horizFing"};
  for(int i = 0; i < 6; i++){
    myPad->cd(i + 1);
    
    TGraph *g_l1_b1 = (TGraph *)infile->Get(Form("%s_beam1_las1_%s_time", tree_name.Data(), msmts[i].Data()));
    TGraph *g_l0_b1 = (TGraph *)infile->Get(Form("%s_beam1_las0_%s_time", tree_name.Data(), msmts[i].Data()));
    TGraph *g_l1_b0 = (TGraph *)infile->Get(Form("%s_beam0_las1_%s_time", tree_name.Data(), msmts[i].Data()));
    TGraph *g_l0_b0 = (TGraph *)infile->Get(Form("%s_beam0_las0_%s_time", tree_name.Data(), msmts[i].Data()));

    g_l1_b1->SetMarkerColor(kBlack); g_l1_b1->SetMarkerStyle(6);
    g_l0_b1->SetMarkerColor(kBlack); g_l0_b1->SetMarkerStyle(6);
    g_l1_b0->SetMarkerColor(kRed); g_l1_b0->SetMarkerStyle(3);
    g_l0_b0->SetMarkerColor(kRed); g_l0_b0->SetMarkerStyle(3);

    TString title = Form("Run %i, %s: %s vs time", run_num, tree_name.Data(), msmts[i].Data());
    g_l1_b1->SetTitle(title.Data()); g_l1_b0->SetTitle(title.Data());
    g_l0_b1->SetTitle(title.Data()); g_l0_b0->SetTitle(title.Data());

    g_l1_b1->Draw("ap"); g_l1_b0->Draw("p && same");
    g_l0_b1->Draw("p && same"); g_l0_b0->Draw("p && same");
  }

  if(printfile)
    myPad->Print(Form("%s/back_rates.png", output_path.Data()), "png");
}

#endif
