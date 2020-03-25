#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

Double_t get_analyzing_power(int run_num){
  if(run_num < 3800){return 0.1105;}
  else if(run_num >= 4232 && run_num <= 4929){return 0.01655915;}
  else if(run_num >= 4930){return 0.036017934;}
  else{
    printf("Run doesn't have a defined analyzing power.\n");
    return 0.0;
  }
}

vector<vector<int>> findDivisions(ifstream &infile, int run_num, int cycleMode){
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
    if((limits[1] - limits[0]) < 3*60*240 && !cycleMode) continue;
    miniruns_in_this_run++;
    vector<int> minirun; minirun.push_back(run_num); minirun.push_back(miniruns_in_this_run);
    minirun.push_back(limits[0]); minirun.push_back(limits[1]);
    if(cycleMode){
      minirun.push_back(limits[2]); minirun.push_back(limits[3]);
      minirun.push_back(limits[4]); minirun.push_back(limits[5]);
    }
    run.push_back(minirun);
  }
  return run;
}

void backgroundPlots(TCanvas *cPosAcc, TCanvas *cNegAcc, vector<vector<vector<int>>> all_runs, TString snail_fname, int acc, int n_miniruns){
  TH1D *hPosBef = new TH1D(Form("hPosBef_acc%i", acc), Form("Laser PosHel Background by cycle: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hPosBef->GetXaxis()->SetTitle("cycle"); hPosBef->GetYaxis()->SetTitle(Form("PosHelAcc%i", acc)); hPosBef->SetMarkerStyle(22); hPosBef->SetMarkerColor(kRed);
  TH1D *hPosAft = new TH1D(Form("hPosAft_acc%i", acc), Form("Laser PosHel Background by cycle: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hPosAft->GetXaxis()->SetTitle("cycle"); hPosAft->GetYaxis()->SetTitle(Form("PosHelAcc%i", acc)); hPosAft->SetMarkerStyle(23); hPosAft->SetMarkerColor(kRed + 2);
  TH1D *hNegBef = new TH1D(Form("hNegBef_acc%i", acc), Form("Laser NegHel Background by cycle: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hNegBef->GetXaxis()->SetTitle("cycle"); hNegBef->GetYaxis()->SetTitle(Form("NegHelAcc%i", acc)); hNegBef->SetMarkerStyle(22); hNegBef->SetMarkerColor(kRed);
  TH1D *hNegAft = new TH1D(Form("hNegAft_acc%i", acc), Form("Laser NegHel Background by cycle: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hNegAft->GetXaxis()->SetTitle("cycle"); hNegAft->GetYaxis()->SetTitle(Form("NegHelAcc%i", acc)); hNegAft->SetMarkerStyle(23); hNegAft->SetMarkerColor(kRed + 2);

  int tot_miniruns = 1;
  for(int run = 0; run < all_runs.size(); run++){
    for(int minirun = 0; minirun < all_runs[run].size(); minirun++){
      int factor = 1;
      TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), all_runs[run][minirun][0]));
      TTree *quartetwise = (TTree *)f->Get("quartetwise");
      TString posH("PosHelAcc0/PosHelNSamples0"); 
      TString negH("NegHelAcc0/NegHelNSamples0");
      TString hOFF1_NamePosH(Form("hOFF1_acc%i_PosH", acc)); TString hOFF2_NamePosH(Form("hOFF2_acc%i_PosH", acc)); 
      TString hOFF1_NameNegH(Form("hOFF1_acc%i_NegH", acc)); TString hOFF2_NameNegH(Form("hOFF2_acc%i_NegH", acc));
      if(acc > 0){
        posH = Form("PosHelAcc%i", acc);
        negH = Form("NegHelAcc%i", acc);
      }
      TString off_cuts("(laserState==2 || laserState==3) && beamState==1");
      TString meas[2] = {posH, negH};
      TString names[4] = {hOFF1_NamePosH, hOFF2_NamePosH, hOFF1_NameNegH, hOFF2_NameNegH};
      int inds[4] = {2, 3, 6, 7};
      for(int m = 0; m < 2; m++){
        for(int las = 0; las < 2; las++){
          quartetwise->Project(names[2*m + las].Data(), meas[m].Data(), 
              Form("%s && (firstMPSnumber>=%i && firstMPSnumber<=%i)", off_cuts.Data(), all_runs[run][minirun][inds[2*las]], all_runs[run][minirun][inds[2*las + 1]]));
        }
      }
      
      TH1F *hOFF1_PosH = (TH1F *)gDirectory->Get(hOFF1_NamePosH.Data()); TH1F *hOFF2_PosH = (TH1F *)gDirectory->Get(hOFF2_NamePosH.Data());
      TH1F *hOFF1_NegH = (TH1F *)gDirectory->Get(hOFF1_NameNegH.Data()); TH1F *hOFF2_NegH = (TH1F *)gDirectory->Get(hOFF2_NameNegH.Data());

      Double_t vOFF1Pos = hOFF1_PosH->GetMean(); Double_t vOFF1PosE = hOFF1_PosH->GetMeanError();
      Double_t vOFF2Pos = hOFF2_PosH->GetMean(); Double_t vOFF2PosE = hOFF2_PosH->GetMeanError();
      Double_t vOFF1Neg = hOFF1_NegH->GetMean(); Double_t vOFF1NegE = hOFF1_NegH->GetMeanError();
      Double_t vOFF2Neg = hOFF2_NegH->GetMean(); Double_t vOFF2NegE = hOFF2_NegH->GetMeanError();

      cout<<"Analyzing background of run "<<all_runs[run][minirun][0]<<", cycle: "<<all_runs[run][minirun][1]<<"; Pos Mean 1: "<<vOFF1Pos<<"; Pos Mean 2: "<<vOFF2Pos<<endl;
      
      hPosBef->SetBinContent(tot_miniruns, vOFF1Pos); hPosBef->SetBinError(tot_miniruns, vOFF1PosE);
      hPosAft->SetBinContent(tot_miniruns, vOFF2Pos); hPosAft->SetBinError(tot_miniruns, vOFF2PosE);
      hNegBef->SetBinContent(tot_miniruns, vOFF1Neg); hNegBef->SetBinError(tot_miniruns, vOFF1NegE);
      hNegAft->SetBinContent(tot_miniruns, vOFF2Neg); hNegAft->SetBinError(tot_miniruns, vOFF2NegE);
      hPosBef->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));
      hPosAft->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));
      hNegBef->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));
      hNegAft->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));

      tot_miniruns++;
    }
  }
  TLegend *lPos = new TLegend(0.75, 0.75, 0.95, 0.9, "", "NDC");
  lPos->AddEntry(hPosBef, "Pre-Off");
  lPos->AddEntry(hPosAft, "Post-Off");
  TLegend *lNeg = new TLegend(0.75, 0.75, 0.95, 0.9, "", "NDC");
  lNeg->AddEntry(hNegBef, "Pre-Off");
  lNeg->AddEntry(hNegAft, "Post-Off");

  gStyle->SetOptStat(kFALSE);
  cPosAcc->cd();
  hPosBef->GetXaxis()->LabelsOption("v");
  hPosAft->GetXaxis()->LabelsOption("v");
  hPosBef->Draw("P");
  hPosAft->Draw("P && same");
  lPos->Draw("same");
  
  cNegAcc->cd();
  hNegBef->GetXaxis()->LabelsOption("v");
  hNegAft->GetXaxis()->LabelsOption("v");
  hNegBef->Draw("P");
  hNegAft->Draw("P && same");
  lNeg->Draw("same");
}

void minirunAnalysis(TString snail_fname, int acc){
  TCanvas *cPol = new TCanvas(Form("cPol_acc%i", acc), "Polarization Average", 1200, 800);
  TPad *pPol1 = new TPad(Form("pPol1_acc%i", acc), "Pol Avg", 0.0, 0.0, 0.7, 1.0);
  TPad *pPol2 = new TPad(Form("pPol2_acc%i", acc), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
  pPol1->SetGridx(1); pPol1->SetGridy(1);
  pPol1->Draw(); pPol2->Draw();
  TCanvas *cAsym = new TCanvas(Form("cAsym_acc%i", acc), "Asymmetry Average", 1200, 800);
  cAsym->cd()->SetGridx(1); cAsym->cd()->SetGridy(1);

  string run_num_str;
  ifstream infile(Form("%s/%s.list", getenv("COMPMON_SNAILS"), snail_fname.Data()));
  vector<vector<vector<int>>> all_runs; int n_miniruns = 0;
  while(getline(infile, run_num_str)){
    vector<vector<int>> run;
    int run_num = atoi(run_num_str.c_str());
    ifstream minirun_infile(Form("%s/minirun_%i.dat", getenv("COMPMON_MINIRUNS"), run_num));
    run = findDivisions(minirun_infile, run_num, 0);
    all_runs.push_back(run);
    n_miniruns += run.size();
  }

  TH1D *hPol = new TH1D(Form("hPol_acc%i", acc), Form("Polarization by minirun: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hPol->GetXaxis()->SetTitle("minirun"); hPol->GetYaxis()->SetTitle("polarization");
  TH1D *hAsym = new TH1D(Form("hAsym_acc%i", acc), Form("Asymmetry by minirun: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hAsym->GetXaxis()->SetTitle("minirun"); hAsym->GetYaxis()->SetTitle("asymmetry / ppt");
  TH1D *hPol_Pull = new TH1D(Form("hPol_Pull_acc%i", acc), "Polarization Pull Plot", 40, -8, 8);

  TF1 *fconstPol = new TF1(Form("fconstPol_acc%i", acc), "pol0");
  TF1 *fconstAsym = new TF1(Form("fConstAsym_acc%i", acc), "pol0");
  TPaveText *ptPol = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TPaveText *ptAsym = new TPaveText(0.78, 0.62, 0.98, 0.75, "blNDC");
  
  TF1 *fconstOnDiff = new TF1(Form("fConstOnDiff_acc%i", acc), "pol0");
  TPaveText *ptOnDiff = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TF1 *fconstOffDiff = new TF1(Form("fConstOffDiff_acc%i", acc), "pol0");
  TPaveText *ptOffDiff = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TF1 *fconstOnSum = new TF1(Form("fConstOnSum_acc%i", acc), "pol0");
  TPaveText *ptOnSum = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TF1 *fconstOffSum = new TF1(Form("fConstOffSum_acc%i", acc), "pol0");
  TPaveText *ptOffSum = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");

  int tot_miniruns = 1;
  for(int run = 0; run < all_runs.size(); run++){
    for(int minirun = 0; minirun < all_runs[run].size(); minirun++){
      int factor = 1; TString ihwp_state("UNK");
      TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), all_runs[run][minirun][0]));
      TTree *quartetwise = (TTree *)f->Get("quartetwise");
      TString diff("PosHelAcc0/PosHelNSamples0 - NegHelAcc0/NegHelNSamples0"); 
      TString summ("PosHelAcc0/PosHelNSamples0 + NegHelAcc0/NegHelNSamples0");
      TString hON_NameDiff(Form("hON_acc%i_Diff", acc)); TString hOFF_NameDiff(Form("hOFF_acc%i_Diff", acc)); 
      TString hON_NameSumm(Form("hON_acc%i_Summ", acc)); TString hOFF_NameSumm(Form("hOFF_acc%i_Summ", acc));
      if(acc > 0){
        diff = Form("PosHelAcc%i - NegHelAcc%i", acc, acc);
        summ = Form("PosHelAcc%i + NegHelAcc%i", acc, acc);
      }
      TString on_cuts("(laserState==0 || laserState==1) && beamState==1");
      TString off_cuts("(laserState==2 || laserState==3) && beamState==1");
      TString meas[2] = {diff, summ};
      TString lasCuts[2] = {on_cuts, off_cuts};
      TString names[4] = {hON_NameDiff, hOFF_NameDiff, hON_NameSumm, hOFF_NameSumm};
      for(int m = 0; m < 2; m++){
        for(int las = 0; las < 2; las++){
          quartetwise->Project(names[2*m + las].Data(), meas[m].Data(),
              Form("%s && firstMPSnumber>=%i && firstMPSnumber<=%i", lasCuts[las].Data(), all_runs[run][minirun][2], all_runs[run][minirun][3]));
        }
      }
      quartetwise->Draw("epics_ihwp_in>>h_ihwp", "", "goff");
      TH1F *h_ihwp = (TH1F *)gDirectory->Get("h_ihwp");
      Double_t ihwp_in = h_ihwp->GetMean();
      if(ihwp_in >= 0.95){factor = 1; ihwp_state = "IN";}
      else if(ihwp_in <= 0.05){factor = -1; ihwp_state = "OUT";}
      else{cout<<"IHWP State changed mid run! Happened during run "<<all_runs[run][minirun][0]<<"!"<<endl; exit(0);}

      TH1F *hON_Diff = (TH1F *)gDirectory->Get(hON_NameDiff.Data()); TH1F *hOFF_Diff = (TH1F *)gDirectory->Get(hOFF_NameDiff.Data());
      TH1F *hON_Summ = (TH1F *)gDirectory->Get(hON_NameSumm.Data()); TH1F *hOFF_Summ = (TH1F *)gDirectory->Get(hOFF_NameSumm.Data());

      Double_t vOnDiff = hON_Diff->GetMean(); Double_t vOnDiffE = hON_Diff->GetMeanError();
      Double_t vOnSum = hON_Summ->GetMean(); Double_t vOnSumE = hON_Summ->GetMeanError();
      Double_t vOffDiff = hOFF_Diff->GetMean(); Double_t vOffDiffE = hOFF_Diff->GetMeanError();
      Double_t vOffSum = hOFF_Summ->GetMean(); Double_t vOffSumE = hOFF_Summ->GetMeanError();
      Double_t vOnAsym = vOnDiff/vOnSum;    Double_t vOnAsymE = vOnAsym*TMath::Sqrt(TMath::Power(vOnSumE/vOnSum, 2) + TMath::Power(vOnDiffE/vOnDiff, 2));
      Double_t vOffAsym = vOffDiff/vOffSum; Double_t vOffAsymE = vOffAsym*TMath::Sqrt(TMath::Power(vOffSumE/vOffSum, 2) + TMath::Power(vOffDiffE/vOffDiff, 2));
      Double_t vSubSum = vOnSum - vOffSum; Double_t vSubSumE = TMath::Sqrt(TMath::Power(vOnSumE, 2) + TMath::Power(vOffSumE, 2));
      Double_t vSubDiff = vOnDiff - vOffDiff; Double_t vSubDiffE = TMath::Sqrt(TMath::Power(vOnDiffE, 2) + TMath::Power(vOffDiffE, 2));
      Double_t vSubAsym = vSubDiff/vSubSum; Double_t vSubAsymE = vSubAsym*TMath::Sqrt(TMath::Power(vSubSumE/vSubSum, 2) + TMath::Power(vSubDiffE/vSubDiff, 2));
      Double_t pol = vSubAsym/get_analyzing_power(all_runs[run][minirun][0]); Double_t polE = vSubAsymE/get_analyzing_power(all_runs[run][minirun][0]);

      cout<<"Run: "<<all_runs[run][minirun][0]<<", Minirun: "<<all_runs[run][minirun][1]<<"; Start Evt: "<<all_runs[run][minirun][2]<<"; End Evt: "<<all_runs[run][minirun][3]<<"; IHWP: "<<ihwp_state.Data()<<endl;
      hAsym->SetBinContent(tot_miniruns, 1000*vSubAsym); hAsym->SetBinError(tot_miniruns, 1000*vSubAsymE);
      hAsym->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));
      hPol->SetBinContent(tot_miniruns, 100*pol); hPol->SetBinError(tot_miniruns, 100*polE);
      hPol->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));

      tot_miniruns++;
    }
  }
  
  pPol1->cd();
  hPol->SetStats(kFALSE);
  hPol->GetXaxis()->LabelsOption("v");
  hPol->Fit(fconstPol, "Q", "", 0, tot_miniruns);
  ptPol->AddText(Form("--------Polarization--------"));
  ptPol->AddText(Form("%.3f%% +/- %.3f",fconstPol->GetParameter(0), fconstPol->GetParError(0)));
  ptPol->AddText(Form("Rel. Err: %.3f%%", fconstPol->GetParError(0)*100.0/fconstPol->GetParameter(0)));
  //ptPol->AddText(Form("Assumed An. Pow: %.5f", get_analyzing_power(4300)));
  ptPol->AddText(Form("Chi^2 / ndf: %f / %d", fconstPol->GetChisquare(), fconstPol->GetNDF()));
  ptPol->SetBorderSize(1); ptPol->SetFillColor(0);
  hPol->Draw("P");
  ptPol->Draw();

  for(int i = 1; i <= n_miniruns; i++){
    hPol_Pull->Fill((hPol->GetBinContent(i) - fconstPol->Eval(hPol->GetBinCenter(i))) / hPol->GetBinError(i));
  }
  pPol2->cd(); hPol_Pull->Draw();

  cAsym->cd();
  hAsym->GetXaxis()->LabelsOption("v");
  hAsym->Fit(fconstAsym, "Q", "", 0, tot_miniruns);
  ptAsym->AddText(Form("--------Asymmetry--------"));
  ptAsym->AddText(Form("%.3f +/- %.3f (ppt)",fconstAsym->GetParameter(0), fconstAsym->GetParError(0)));
  ptAsym->AddText(Form("Rel. Err: %.3f%%", fconstAsym->GetParError(0)*100.0/fconstAsym->GetParameter(0)));
  ptAsym->AddText(Form("Chi^2 / ndf: %f / %d", fconstAsym->GetChisquare(), fconstAsym->GetNDF()));
  ptAsym->SetBorderSize(1); ptAsym->SetFillColor(0);
  hAsym->Draw("P");
  ptAsym->Draw();

  TString output_dir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snail_fname.Data()));
  cPol->Print(Form("%s/polarization_acc%i.pdf(", output_dir.Data(), acc), "pdf");
  cAsym->Print(Form("%s/polarization_acc%i.pdf)", output_dir.Data(), acc), "pdf");
}

void cycleAnalysis(TString snail_fname, int acc){
  TCanvas *cPol = new TCanvas(Form("cPol_acc%i", acc), "Polarization Average", 1200, 800);
  TPad *pPol1 = new TPad(Form("pPol1_acc%i", acc), "Pol Avg", 0.0, 0.0, 0.7, 1.0);
  TPad *pPol2 = new TPad(Form("pPol2_acc%i", acc), "Pol Pull Plot", 0.7, 0.0, 1.0, 1.0);
  pPol1->SetGridx(1); pPol1->SetGridy(1);
  pPol1->Draw(); pPol2->Draw();
  TCanvas *cAsym = new TCanvas(Form("cAsym_acc%i", acc), "Asymmetry Average", 1200, 800);
  cAsym->cd()->SetGridx(1); cAsym->cd()->SetGridy(1);
  TCanvas *cPos = new TCanvas(Form("cPos_acc%i", acc), "Pos Hel Back Average", 1200, 800);
  cPos->cd()->SetGridx(1); cPos->cd()->SetGridy(1);
  TCanvas *cNeg = new TCanvas(Form("cNeg_acc%i", acc), "Neg Hel Back Average", 1200, 800);
  cNeg->cd()->SetGridx(1); cNeg->cd()->SetGridy(1);

  string run_num_str;
  ifstream infile(Form("%s/%s.list", getenv("COMPMON_SNAILS"), snail_fname.Data()));
  vector<vector<vector<int>>> all_runs; int n_miniruns = 0;
  while(getline(infile, run_num_str)){
    vector<vector<int>> run;
    int run_num = atoi(run_num_str.c_str());
    ifstream cycle_infile(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), run_num));
    run = findDivisions(cycle_infile, run_num, 1);
    all_runs.push_back(run);
    n_miniruns += run.size();
  }

  TH1D *hPol = new TH1D(Form("hPol_acc%i", acc), Form("Polarization by cycle: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hPol->GetXaxis()->SetTitle("cycle"); hPol->GetYaxis()->SetTitle("polarization");
  TH1D *hAsym = new TH1D(Form("hAsym_acc%i", acc), Form("Asymmetry by cycle: %s (Acc%i)", snail_fname.Data(), acc), n_miniruns, 0, n_miniruns);
  hAsym->GetXaxis()->SetTitle("cycle"); hAsym->GetYaxis()->SetTitle("asymmetry / ppt");
  TH1D *hPol_Pull = new TH1D(Form("hPol_Pull_acc%i", acc), "Polarization Pull Plot", 40, -8, 8);

  TF1 *fconstPol = new TF1(Form("fconstPol_acc%i", acc), "pol0");
  TF1 *fconstAsym = new TF1(Form("fConstAsym_acc%i", acc), "pol0");
  TPaveText *ptPol = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TPaveText *ptAsym = new TPaveText(0.78, 0.62, 0.98, 0.75, "blNDC");
  
  TF1 *fconstOnDiff = new TF1(Form("fConstOnDiff_acc%i", acc), "pol0");
  TPaveText *ptOnDiff = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TF1 *fconstOffDiff = new TF1(Form("fConstOffDiff_acc%i", acc), "pol0");
  TPaveText *ptOffDiff = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TF1 *fconstOnSum = new TF1(Form("fConstOnSum_acc%i", acc), "pol0");
  TPaveText *ptOnSum = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");
  TF1 *fconstOffSum = new TF1(Form("fConstOffSum_acc%i", acc), "pol0");
  TPaveText *ptOffSum = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");

  int tot_miniruns = 1; bool same_ihwp_state = true; int ihwp_const = -1;
  for(int run = 0; run < all_runs.size(); run++){
    for(int minirun = 0; minirun < all_runs[run].size(); minirun++){
      int factor = 1; TString ihwp_state("UNK");
      TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), all_runs[run][minirun][0]));
      TTree *quartetwise = (TTree *)f->Get("quartetwise");
      TString diff("PosHelAcc0/PosHelNSamples0 - NegHelAcc0/NegHelNSamples0"); 
      TString summ("PosHelAcc0/PosHelNSamples0 + NegHelAcc0/NegHelNSamples0");
      TString hON_NameDiff(Form("hON_acc%i_Diff", acc)); TString hOFF_NameDiff(Form("hOFF_acc%i_Diff", acc)); 
      TString hON_NameSumm(Form("hON_acc%i_Summ", acc)); TString hOFF_NameSumm(Form("hOFF_acc%i_Summ", acc));
      if(acc > 0){
        diff = Form("PosHelAcc%i - NegHelAcc%i", acc, acc);
        summ = Form("PosHelAcc%i + NegHelAcc%i", acc, acc);
      }
      TString on_cuts("(laserState==0 || laserState==1) && beamState==1");
      TString off_cuts("(laserState==2 || laserState==3) && beamState==1");
      TString meas[2] = {diff, summ};
      TString lasCuts[2] = {on_cuts, off_cuts};
      TString names[4] = {hON_NameDiff, hOFF_NameDiff, hON_NameSumm, hOFF_NameSumm};
      for(int m = 0; m < 2; m++){
        for(int las = 0; las < 2; las++){
          quartetwise->Project(names[2*m + las].Data(), meas[m].Data(), 
              Form("%s && ((firstMPSnumber>=%i && firstMPSnumber<=%i) || (firstMPSnumber>=%i && firstMPSnumber<=%i) || (firstMPSnumber>=%i && firstMPSnumber<=%i))",
              lasCuts[las].Data(), all_runs[run][minirun][2], all_runs[run][minirun][3], all_runs[run][minirun][4], all_runs[run][minirun][5], all_runs[run][minirun][6], all_runs[run][minirun][7]));
        }
      }
      quartetwise->Draw("epics_ihwp_in>>h_ihwp", "", "goff");
      TH1F *h_ihwp = (TH1F *)gDirectory->Get("h_ihwp");
      Double_t ihwp_in = h_ihwp->GetMean();
      if(ihwp_in >= 0.95){factor = 1; ihwp_state = "IN";}
      else if(ihwp_in <= 0.05){factor = -1; ihwp_state = "OUT";}
      else{cout<<"IHWP State changed mid run! Happened during run "<<all_runs[run][minirun][0]<<"!"<<endl; exit(0);}
      if(run == 0 && minirun == 0){
        ihwp_const = (int)(0.5*(factor + 1));
      }
      else{
        if(ihwp_const == 0 && factor == 1){same_ihwp_state = false;}
        else if(ihwp_const == 1 && factor == -1){same_ihwp_state = false;}
      }

      TH1F *hON_Diff = (TH1F *)gDirectory->Get(hON_NameDiff.Data()); TH1F *hOFF_Diff = (TH1F *)gDirectory->Get(hOFF_NameDiff.Data());
      TH1F *hON_Summ = (TH1F *)gDirectory->Get(hON_NameSumm.Data()); TH1F *hOFF_Summ = (TH1F *)gDirectory->Get(hOFF_NameSumm.Data());

      Double_t vOnDiff = hON_Diff->GetMean(); Double_t vOnDiffE = hON_Diff->GetMeanError();
      Double_t vOnSum = hON_Summ->GetMean(); Double_t vOnSumE = hON_Summ->GetMeanError();
      Double_t vOffDiff = hOFF_Diff->GetMean(); Double_t vOffDiffE = hOFF_Diff->GetMeanError();
      Double_t vOffSum = hOFF_Summ->GetMean(); Double_t vOffSumE = hOFF_Summ->GetMeanError();
      Double_t vOnAsym = vOnDiff/vOnSum;    Double_t vOnAsymE = vOnAsym*TMath::Sqrt(TMath::Power(vOnSumE/vOnSum, 2) + TMath::Power(vOnDiffE/vOnDiff, 2));
      Double_t vOffAsym = vOffDiff/vOffSum; Double_t vOffAsymE = vOffAsym*TMath::Sqrt(TMath::Power(vOffSumE/vOffSum, 2) + TMath::Power(vOffDiffE/vOffDiff, 2));
      Double_t vSubSum = vOnSum - vOffSum; Double_t vSubSumE = TMath::Sqrt(TMath::Power(vOnSumE, 2) + TMath::Power(vOffSumE, 2));
      Double_t vSubDiff = vOnDiff - vOffDiff; Double_t vSubDiffE = TMath::Sqrt(TMath::Power(vOnDiffE, 2) + TMath::Power(vOffDiffE, 2));
      Double_t vSubAsym = vSubDiff/vSubSum; Double_t vSubAsymE = vSubAsym*TMath::Sqrt(TMath::Power(vSubSumE/vSubSum, 2) + TMath::Power(vSubDiffE/vSubDiff, 2));
      Double_t pol = vSubAsym/get_analyzing_power(all_runs[run][minirun][0]); Double_t polE = vSubAsymE/get_analyzing_power(all_runs[run][minirun][0]);
      cout<<"Run: "<<all_runs[run][minirun][0]<<", Laser Cycle: "<<all_runs[run][minirun][1]<<"; Start Evt: "<<all_runs[run][minirun][2]<<"; End Evt: "<<all_runs[run][minirun][7]<<"; IHWP Factor: "<<factor<<"; Asym: "<<vSubAsym<<endl;
      hAsym->SetBinContent(tot_miniruns, 1000*vSubAsym); hAsym->SetBinError(tot_miniruns, 1000*vSubAsymE);
      hAsym->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));
      hPol->SetBinContent(tot_miniruns, 100*pol); hPol->SetBinError(tot_miniruns, 100*polE);
      hPol->GetXaxis()->SetBinLabel(tot_miniruns, Form("%i.%i", all_runs[run][minirun][0], all_runs[run][minirun][1]));

      tot_miniruns++;
    }
  }
  
  if(same_ihwp_state){
    TString state = (ihwp_const == 1) ? "IHWP IN" : "IHWP OUT";
    hPol->SetTitle(Form("Polarization by cycle: %s (Acc%i, %s)", snail_fname.Data(), acc, state.Data()));
    hAsym->SetTitle(Form("Asymmetry by cycle: %s (Acc%i, %s)", snail_fname.Data(), acc, state.Data()));
  }
  else{
    printf("Could not see a common ihwp state!\n");
  }
  pPol1->cd();
  hPol->SetStats(kFALSE);
  hPol->GetXaxis()->LabelsOption("v");
  hPol->Fit(fconstPol, "Q", "", 0, tot_miniruns);
  ptPol->AddText(Form("--------Polarization--------"));
  ptPol->AddText(Form("%.3f%% +/- %.3f",fconstPol->GetParameter(0), fconstPol->GetParError(0)));
  ptPol->AddText(Form("Rel. Err: %.3f%%", fconstPol->GetParError(0)*100.0/fconstPol->GetParameter(0)));
  //ptPol->AddText(Form("Assumed An. Pow: %.5f", get_analyzing_power(4300)));
  ptPol->AddText(Form("Chi^2 / ndf: %f / %d", fconstPol->GetChisquare(), fconstPol->GetNDF()));
  ptPol->SetBorderSize(1); ptPol->SetFillColor(0);
  hPol->Draw("P");
  ptPol->Draw();

  for(int i = 1; i <= n_miniruns; i++){
    hPol_Pull->Fill((hPol->GetBinContent(i) - fconstPol->Eval(hPol->GetBinCenter(i))) / hPol->GetBinError(i));
  }
  pPol2->cd(); hPol_Pull->Draw();

  cAsym->cd();
  hAsym->GetXaxis()->LabelsOption("v");
  hAsym->Fit(fconstAsym, "Q", "", 0, tot_miniruns);
  ptAsym->AddText(Form("--------Asymmetry--------"));
  ptAsym->AddText(Form("%.3f +/- %.3f (ppt)",fconstAsym->GetParameter(0), fconstAsym->GetParError(0)));
  ptAsym->AddText(Form("Rel. Err: %.3f%%", fconstAsym->GetParError(0)*100.0/fconstAsym->GetParameter(0)));
  ptAsym->AddText(Form("Chi^2 / ndf: %f / %d", fconstAsym->GetChisquare(), fconstAsym->GetNDF()));
  ptAsym->SetBorderSize(1); ptAsym->SetFillColor(0);
  hAsym->Draw("P");
  ptAsym->Draw();

  backgroundPlots(cPos, cNeg, all_runs, snail_fname, acc, n_miniruns);

  TString output_dir(Form("%s/snails/%s", getenv("COMPMON_WEB"), snail_fname.Data()));
  cPol->Print(Form("%s/polarization_acc%i_cycles.pdf(", output_dir.Data(), acc), "pdf");
  cAsym->Print(Form("%s/polarization_acc%i_cycles.pdf", output_dir.Data(), acc), "pdf");
  cPos->Print(Form("%s/polarization_acc%i_cycles.pdf",  output_dir.Data(), acc), "pdf");
  cNeg->Print(Form("%s/polarization_acc%i_cycles.pdf)", output_dir.Data(), acc), "pdf");
}

/**
===================================================
||                                               ||
||  |\    /|       /\       ========== |\    ||  ||
||  ||\  /||      //\\          ||     ||\   ||  ||
||  ||\\//||     //  \\         ||     ||\\  ||  ||
||  || \/ ||    //====\\        ||     || \\ ||  ||
||  ||    ||   //      \\       ||     ||  \\||  ||
||  ||    ||  //        \\      ||     ||   \||  ||
||  ||    || //          \\ ========== ||    \|  ||
||                                               ||
===================================================
**/
void aggregate(TString snail_fname, int acc=0, int cycleMode=0){
  if(cycleMode)
    cycleAnalysis(snail_fname, acc);
  else
    minirunAnalysis(snail_fname, acc);
}
