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

Double_t get_analyzing_power(int run_num){
  if(run_num < 3800){return 0.1105;}
  else if(run_num >= 4232){return 0.01655915;}
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

TString date_time_string(TFile *infile){
  Int_t day      = (Int_t)((TH1F *)infile->Get("h_day"))->GetMean();
  Int_t month    = (Int_t)((TH1F *)infile->Get("h_month"))->GetMean();
  Int_t date_num = (Int_t)((TH1F *)infile->Get("h_date_num"))->GetMean();
  Int_t year     = (Int_t)((TH1F *)infile->Get("h_year"))->GetMean();
  Int_t hour     = (Int_t)((TH1F *)infile->Get("h_hour"))->GetMean();
  Int_t minute   = (Int_t)((TH1F *)infile->Get("h_minute"))->GetMean();
  Int_t second   = (Int_t)((TH1F *)infile->Get("h_second"))->GetMean();

  TString month_str(Form("%i", month)); TString date_num_str(Form("%i", date_num)); TString hour_str(Form("%i", hour));
  TString minute_str(Form("%i", minute)); TString second_str(Form("%i", second));
  if(month < 10){month_str = Form("0%i", month);}
  if(date_num < 10){date_num_str = Form("0%i", date_num);}
  if(hour < 10){hour_str = Form("0%i", hour);}
  if(minute < 10){minute_str = Form("0%i", minute);}
  if(second < 10){second_str = Form("0%i", second);}
  return Form("%i-%s-%s %s:%s:%s", year, month_str.Data(), date_num_str.Data(), hour_str.Data(), minute_str.Data(), second_str.Data());
}

std::vector<TString> calc_polarization(std::vector<Double_t> factors, int run_num){
  Double_t vOnPos = factors[0]; Double_t vOnPosE = factors[1];
  Double_t vOffPos = factors[2]; Double_t vOffPosE = factors[3];
  Double_t vOnNeg = factors[4]; Double_t vOnNegE = factors[5];
  Double_t vOffNeg = factors[6]; Double_t vOffNegE = factors[7];
  Double_t vOnDiff = factors[8]; Double_t vOnDiffE = factors[9];
  Double_t vOffDiff = factors[10]; Double_t vOffDiffE = factors[11];
  Double_t vOnSum = factors[12]; Double_t vOnSumE = factors[13];
  Double_t vOffSum = factors[14]; Double_t vOffSumE = factors[15];

  //Double_t vOnSum = vOnPos + vOnNeg; Double_t vOnSumE = TMath::Sqrt(TMath::Power(vOnPosE, 2) + TMath::Power(vOnNegE, 2)); 
  //Double_t vOnDiff = vOnPos - vOnNeg; Double_t vOnDiffE = TMath::Sqrt(TMath::Power(vOnPosE, 2) + TMath::Power(vOnNegE, 2));
  //Double_t vOffSum = vOffPos + vOffNeg; Double_t vOffSumE = TMath::Sqrt(TMath::Power(vOffPosE, 2) + TMath::Power(vOffNegE, 2)); 
  //Double_t vOffDiff = vOffPos - vOffNeg; Double_t vOffDiffE = TMath::Sqrt(TMath::Power(vOffPosE, 2) + TMath::Power(vOffNegE, 2));

  Double_t vOnAsym = vOnDiff/vOnSum; Double_t vOnAsymE = vOnAsym*TMath::Sqrt(TMath::Power(vOnSumE/vOnSum, 2) + TMath::Power(vOnDiffE/vOnDiff, 2));
  Double_t vOffAsym = vOffDiff/vOffSum; Double_t vOffAsymE = vOffAsym*TMath::Sqrt(TMath::Power(vOffSumE/vOffSum, 2) + TMath::Power(vOffDiffE/vOffDiff, 2));

  Double_t vSubSum = vOnSum - vOffSum; Double_t vSubSumE = TMath::Sqrt(TMath::Power(vOnSumE, 2) + TMath::Power(vOffSumE, 2));
  Double_t vSubDiff = vOnDiff - vOffDiff; Double_t vSubDiffE = TMath::Sqrt(TMath::Power(vOnDiffE, 2) + TMath::Power(vOffDiffE, 2));
  //Double_t vSubDiff = vOnDiff; Double_t vSubDiffE = vOnDiffE;
  Double_t vSubAsym = vSubDiff/vSubSum; Double_t vSubAsymE = vSubAsym*TMath::Sqrt(TMath::Power(vSubSumE/vSubSum, 2) + TMath::Power(vSubDiffE/vSubDiff, 2));

  Double_t pol = vSubAsym/get_analyzing_power(run_num); Double_t polE = vSubAsymE/get_analyzing_power(run_num);
  /** TString results("--------Histogram Vars--------\n"); 
  results += Form("    Positive Hel (ON): %f +/- %f\n",  vOnPos,  vOnPosE);  results += Form("    Negative Hel (ON): %f +/- %f\n",  vOnNeg,  vOnNegE);
  results += Form("    Positive Hel (OFF): %f +/- %f\n", vOffPos, vOffPosE); results += Form("    Negative Hel (OFF): %f +/- %f\n", vOffNeg, vOffNegE);
  results += "--------Sums and Diffs--------\n";
  results += Form("    On Diff: %f +/- %f\n", vOnDiff, vOnDiffE); results += Form("    Off Diff: %f +/- %f\n", vOffDiff, vOffDiffE);
  results += Form("    On Sum: %f +/- %f\n",  vOnSum,  vOnSumE);  results += Form("    Off Sum: %f +/- %f\n",  vOffSum,  vOffSumE);
  results += "--------Asymmetries--------\n";
  results += Form("    On Asym: %f +/- %f\n", vOnAsym, vOnAsymE); results += Form("    Off Asym: %f +/- %f\n", vOffAsym, vOffAsymE);
  results += "--------Background Subtracted--------\n";
  results += Form("    Back Sub Diff: %f +/- %f\n", vSubDiff, vSubDiffE); results += Form("    Back Sub Sum: %f +/- %f\n", vSubSum, vSubSumE);
  results += Form("    Back Sub Asym: %f +/- %f\n", vSubAsym, vSubAsymE);
  results += "--------Polarization--------\n";
  results += Form("    Polarization: %f +/- %f", 100*pol, 100*polE); **/
  std::vector<TString> results; results.push_back("--------Histogram Vars--------\n"); 
  results.push_back(Form("    Positive Hel (ON): %f +/- %f\n",  vOnPos,  vOnPosE));  results.push_back(Form("    Negative Hel (ON): %f +/- %f\n",  vOnNeg,  vOnNegE));
  results.push_back(Form("    Positive Hel (OFF): %f +/- %f\n", vOffPos, vOffPosE)); results.push_back(Form("    Negative Hel (OFF): %f +/- %f\n", vOffNeg, vOffNegE));
  results.push_back("--------Sums and Diffs--------\n");
  results.push_back(Form("    On Diff: %f +/- %f\n", vOnDiff, vOnDiffE)); results.push_back(Form("    Off Diff: %f +/- %f\n", vOffDiff, vOffDiffE));
  results.push_back(Form("    On Sum: %f +/- %f\n",  vOnSum,  vOnSumE));  results.push_back(Form("    Off Sum: %f +/- %f\n",  vOffSum,  vOffSumE));
  results.push_back("--------Asymmetries--------\n");
  results.push_back(Form("    On Asym: %f +/- %f\n", vOnAsym, vOnAsymE)); results.push_back(Form("    Off Asym: %f +/- %f\n", vOffAsym, vOffAsymE));
  results.push_back("--------Background Subtracted--------\n");
  results.push_back(Form("    Back Sub Diff: %f +/- %f\n", vSubDiff, vSubDiffE)); results.push_back(Form("    Back Sub Sum: %f +/- %f\n", vSubSum, vSubSumE));
  results.push_back(Form("    Back Sub Asym: %f +/- %f\n", vSubAsym, vSubAsymE));
  results.push_back("--------Polarization--------\n");
  results.push_back(Form("    Polarization: %f +/- %f (%f%%)", 100*pol, 100*polE, 100*polE/pol));
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

  Double_t seconds = acc0_all->GetEntries()*1.0/(240.0); Int_t iSeconds = (Int_t)seconds;
  Double_t minutes = seconds/60.0; Int_t iMinutes = (Int_t)minutes;
  Double_t hours = minutes/60.0; Int_t iHours = (Int_t)hours;
  Int_t realMinutes = (Int_t)((hours - iHours)*60.0);
  Int_t realSeconds = (Int_t)((minutes - iMinutes)*60.0);

  Int_t lON_evts = acc0_bON_lON->GetEntries();
  Int_t lOFF_evts = acc0_bON_lOFF->GetEntries();
  Int_t bOFF_evts = acc0_bOFF->GetEntries();
  Int_t tot_evts = lON_evts + lOFF_evts + bOFF_evts;

  Float_t qw1_deg = fmod((-0.004*h_qw1->GetMean() - 1390.768), 360);
  Float_t hw1_deg = fmod((-0.004*h_hw1->GetMean()), 360);
  
  TPaveText *essStats = new TPaveStats(0.0, 0.0, 1.0, 1.0, "blNDC");
  essStats->AddText(Form("--------Essential Statistics: Run %i--------", run_num));
  essStats->AddText(Form("First Epics Time: %s", date_time_string(infile).Data()));
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
    myPad->Print(Form("%s/ess_stats.pdf", output_path.Data()), "pdf");
}

void snapshots_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, bool printfile=false){
  TString tree_name("snapshots");
  TString data_name("sumPeakHeight");
  TString acc_name("accepted");
  TString rej_names[5] = {"rej_peak", "rej_sat", "rej_corrPed", "rej_narrPed", "rej_time"};
  TString rej_titles[5] = {"Peak Height Cut", "Saturation Cut", "Correct Pedestal Cut", "Narrow Pedestal Cut", "Time Cut"};
  gStyle->SetOptStat(0);
  myPad->cd();
  TGraph *g = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), data_name.Data()));
  TGraph *g_acc = (TGraph *)infile->Get(Form("%s_%s", tree_name.Data(), acc_name.Data()));
  g->SetMarkerColor(2);
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

  if(printfile)
    myPad->Print(Form("%s/snapshots.pdf", output_path.Data()), "pdf");
}

void breakdown_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, TString tree_name, TString data_name, bool printfile=false){
  myPad->Divide(2, 2);

  TH1F *h_all = (TH1F *)infile->Get(Form("%s_%s", tree_name.Data(), data_name.Data())); myPad->cd(1); h_all->Draw();
  std::vector<TString> s1 = hist_stats(h_all);
  TPaveText *pt1 = new TPaveText(0.75, 0.7, 0.95, 0.9, "blNDC"); pt1->SetBorderSize(1); pt1->SetFillColor(0);
  for(int i = 0; i < s1.size(); i++){pt1->AddText(s1[i].Data())->SetTextColor(4);} pt1->Draw("same"); 
  TH1F *h_beamOff = (TH1F *)infile->Get(Form("%s_beam0_%s", tree_name.Data(), data_name.Data())); myPad->cd(2); h_beamOff->Draw();
  std::vector<TString> s2 = hist_stats(h_beamOff);
  TPaveText *pt2 = new TPaveText(0.75, 0.7, 0.95, 0.9, "blNDC"); pt2->SetBorderSize(1); pt2->SetFillColor(0);
  for(int i = 0; i < s2.size(); i++){pt2->AddText(s2[i].Data())->SetTextColor(4);} pt2->Draw("same"); 

  TH1F *h_l1_h0 = (TH1F *)infile->Get(Form("%s_beam1_las1_hel0_%s", tree_name.Data(), data_name.Data()));
  TH1F *h_l1_h1 = (TH1F *)infile->Get(Form("%s_beam1_las1_hel1_%s", tree_name.Data(), data_name.Data()));
  h_l1_h0->SetLineColor(3); h_l1_h1->SetLineColor(2); myPad->cd(3);
  THStack *hs3 = new THStack(Form("%s_beam1_las1_%s_stack", tree_name.Data(), data_name.Data()), 
      Form("Run %i %s, %s: Beam ON, Laser ON", run_num, tree_name.Data(), data_name.Data()));
  hs3->Add(h_l1_h0); hs3->Add(h_l1_h1); hs3->Draw("nostack");
  TLegend *leg3 = new TLegend(0.75, 0.75, 0.95, 0.9, "", "NDC");
  leg3->AddEntry(h_l1_h0, "Helicity==0"); leg3->AddEntry(h_l1_h1, "Helicity==1"); leg3->Draw("same");
  std::vector<TString> s_l1_h0 = hist_stats(h_l1_h0); std::vector<TString> s_l1_h1 = hist_stats(h_l1_h1);
  TPaveText *pt3_0 = new TPaveText(0.75, 0.55, 0.95, 0.75, "blNDC"); pt3_0->SetBorderSize(1); pt3_0->SetFillColor(0);
  TPaveText *pt3_1 = new TPaveText(0.75, 0.35, 0.95, 0.55, "blNDC"); pt3_1->SetBorderSize(1); pt3_1->SetFillColor(0);
  for(int i = 0; i < s_l1_h0.size(); i++){pt3_0->AddText(s_l1_h0[i].Data())->SetTextColor(3);}
  for(int i = 0; i < s_l1_h1.size(); i++){pt3_1->AddText(s_l1_h1[i].Data())->SetTextColor(2);}
  pt3_0->Draw("same"); pt3_1->Draw("same");

  TH1F *h_l0_h0 = (TH1F *)infile->Get(Form("%s_beam1_las0_hel0_%s", tree_name.Data(), data_name.Data()));
  TH1F *h_l0_h1 = (TH1F *)infile->Get(Form("%s_beam1_las0_hel1_%s", tree_name.Data(), data_name.Data()));
  h_l0_h0->SetLineColor(3); h_l0_h1->SetLineColor(2); myPad->cd(4);
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
    myPad->Print(Form("%s/%s.pdf", output_path.Data(), data_name.Data()), "pdf");
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
    myPad->Print(Form("%s/q%s.pdf", output_path.Data(), data_name.Data()), "pdf");
}

void asym_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, int accum=0, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 3);

  std::vector<double> factors;
  TString msmts[5] = {Form("posH%i", accum), Form("negH%i", accum), Form("diff%i", accum), Form("summ%i", accum), Form("asym%i", accum)};
  for(int i = 0;  i < 5; i++){
    myPad->cd(i + 1);
    TH1F *hON  = (TH1F *)infile->Get(Form("%s_beam1_las1_%s", tree_name.Data(), msmts[i].Data()));
    TH1F *hOFF = (TH1F *)infile->Get(Form("%s_beam1_las0_%s", tree_name.Data(), msmts[i].Data()));
    if(i < 4){
      factors.push_back(hON->GetMean()); factors.push_back(hON->GetMeanError());
      factors.push_back(hOFF->GetMean()); factors.push_back(hOFF->GetMeanError());
    }
    hON->SetLineColor(3); hOFF->SetLineColor(2);
    THStack *hs = new THStack(Form("quartet_beam1_%s_stack", msmts[i].Data()), Form("Run %i %s, %s: Beam ON", run_num, tree_name.Data(), msmts[i].Data()));
    hs->Add(hON); hs->Add(hOFF); hs->Draw("nostack");
    TLegend *leg = new TLegend(0.75, 0.75, 0.99, 0.9, "", "NDC");
    leg->AddEntry(hON, "Laser ON"); leg->AddEntry(hOFF, "Laser OFF"); leg->Draw("same");
    std::vector<TString> sON = hist_stats(hON); std::vector<TString> sOFF = hist_stats(hOFF);
    TPaveText *ptON  = new TPaveText(0.75, 0.55, 0.99, 0.75, "blNDC"); ptON->SetBorderSize(1);  ptON->SetFillColor(0);
    TPaveText *ptOFF = new TPaveText(0.75, 0.35, 0.99, 0.55, "blNDC"); ptOFF->SetBorderSize(1); ptOFF->SetFillColor(0);
    for(int i = 0; i < sON.size();  i++){ptON->AddText(sON[i].Data())->SetTextColor(3);}
    for(int i = 0; i < sOFF.size(); i++){ptOFF->AddText(sOFF[i].Data())->SetTextColor(2);}
    ptON->Draw("same"); ptOFF->Draw("same");
  }
  myPad->cd(6);
  std::vector<TString> results = calc_polarization(factors, run_num);
  TPaveText *pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
  for(int i = 0; i < results.size(); i++){pt->AddText(results[i].Data());}
  pt->SetBorderSize(0); pt->SetFillColor(0);
  pt->Draw();

  if(printfile)
    myPad->Print(Form("%s/asym_acc%i_hists.pdf", output_path.Data(), accum), "pdf");
}

void asym_graph_pad(TFile *infile, int run_num, TString output_path, TPad *myPad, int accum=0, bool printfile=false){
  TString tree_name("quartetwise");
  myPad->Divide(2, 3);

  gStyle->SetOptStat(0);
  TString msmts[5] = {Form("posH%i", accum), Form("negH%i", accum), Form("diff%i", accum), Form("summ%i", accum), Form("asym%i", accum)};
  for(int i = 0; i < 5; i++){
    myPad->cd(i + 1);
    TGraph *gON  = (TGraph *)infile->Get(Form("%s_beam1_las1_%s_time", tree_name.Data(), msmts[i].Data()));
    TGraph *gOFF = (TGraph *)infile->Get(Form("%s_beam1_las0_%s_time", tree_name.Data(), msmts[i].Data()));
    gON->SetMarkerColor(3); gOFF->SetMarkerColor(2);
    gON->SetTitle(Form("Run %i %s, %s vs time", run_num, tree_name.Data(), msmts[i].Data()));
    gOFF->SetTitle(Form("Run %i %s, %s vs time", run_num, tree_name.Data(), msmts[i].Data()));
    gON->Draw("ap"); gOFF->Draw("p && same");
  }

  myPad->cd(6);
  TGraph *g_b0_l0 = (TGraph *)infile->Get(Form("%s_beam0_las0_bcm_time", tree_name.Data())); TGraph *g_b0_l1 = (TGraph *)infile->Get(Form("%s_beam0_las1_bcm_time", tree_name.Data()));
  TGraph *g_b1_l0 = (TGraph *)infile->Get(Form("%s_beam1_las0_bcm_time", tree_name.Data())); TGraph *g_b1_l1 = (TGraph *)infile->Get(Form("%s_beam1_las1_bcm_time", tree_name.Data()));
  g_b0_l0->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data())); g_b0_l1->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data()));
  g_b1_l0->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data())); g_b1_l1->SetTitle(Form("Run %i %s, bcm vs time", run_num, tree_name.Data()));
  g_b0_l0->SetMarkerColor(2); g_b0_l1->SetMarkerColor(2);
  g_b1_l0->Draw("ap"); g_b1_l1->Draw("p && same");
  g_b0_l0->Draw("* && same"); g_b0_l1->Draw("* && same"); 

  if(printfile)
    myPad->Print(Form("%s/q_acc%i_graphs.pdf", output_path.Data(), accum), "pdf");
}

void beam_off_quartet_pad(TFile *infile, int run_num, TString output_path, TPad *myPad){
  TString tree_name("quartetwise");
  myPad->Divide(2, 3);

  gStyle->SetOptStat(); std::vector<double> factors;
  TString msmts[5] = {"posH0", "negH0", "diff0", "summ0", "asym0"};
  for(int i = 0;  i < 5; i++){
    myPad->cd(i + 1);
    TH1F *h_bOFF = (TH1F *)infile->Get(Form("%s_beam0_%s", tree_name.Data(), msmts[i].Data()));
    if(i < 3){
      factors.push_back(h_bOFF->GetMean()); factors.push_back(h_bOFF->GetMeanError());
      factors.push_back(0.0); factors.push_back(0.0);
    }
    else if(i < 4){
      factors.push_back(19.17); factors.push_back(0.001599);
      factors.push_back(0.0); factors.push_back(0.0);
    }
    h_bOFF->Draw();
  }
  myPad->cd(6);
  std::vector<TString> results = calc_polarization(factors, run_num);
  TPaveText *pt = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
  for(int i = 0; i < results.size(); i++){pt->AddText(results[i].Data());}
  pt->SetBorderSize(0); pt->SetFillColor(0);
  pt->Draw();
}

#endif
