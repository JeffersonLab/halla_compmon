/**
Code Commissioned by AJZ late June 2019

This code assumes you have a rootfile produced by dataQualityCheck_simple.C in that exact format
Code WILL segfault if this is not true
**/

#include <TCanvas.h>
#include <TObject.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TPad.h>
#include <TString.h>
#include <vector>

#include "utils.h"

using namespace std;

void convert(int run_num, TString output_path){
  gSystem->Exec(Form("convert %s/ess_stats*.png %s/ess_stats.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("convert %s/snapshots.png %s/snapshots.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("convert %s/sums.png %s/sums.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("convert %s/acc0*.png %s/acc0.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("convert %s/quart_*.png %s/quartet.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("convert %s/*_acc*.png %s/asymmetries.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("convert %s/back_*.png %s/backgrounds.pdf", output_path.Data(), output_path.Data()));
  gSystem->Exec(Form("rm -rf %s/*.png", output_path.Data()));
}

void writeToPDF(int run_num){
  TString output_path(Form("%s/runs/Run%i/", getenv("COMPMON_WEB"), run_num));
  TString input_plotfile(Form("%s/compton_online_run_%i.root", getenv("COMPMON_PLOTFILES"), run_num));

  TFile *infile = TFile::Open(input_plotfile.Data());

  TCanvas *c_all = new TCanvas("c_all", "Wrapper Canvas", 1200, 800);
  TPad *c_ess_stats = new TPad("c_ess_stats", "Essential Stats", 0, 0, 1, 1);
  TPad *c_ess_stats_2 = new TPad("c_ess_stats_2", "Essential Stats 2", 0, 0, 1, 1);
  TPad *c_snap = new TPad("c_snap", "Snapshots", 0, 0, 1, 1);
  TPad *c_sums = new TPad("c_sums", "Triggered Sums", 0, 0, 1, 1);
  TPad *c_acc0 = new TPad("c_acc0", "Acc0 Histos", 0, 0, 1, 1);
  TPad *c_acc0_time = new TPad("c_acc0_time", "Acc0 Graphs", 0, 0, 1, 1);
  TPad *c_posH0 = new TPad("c_posH0", "Pos Hs", 0, 0, 1, 1);
  TPad *c_negH0 = new TPad("c_negH0", "Neg Hs", 0, 0, 1, 1);
  TPad *c_diff0 = new TPad("c_diff0", "Diffs", 0, 0, 1, 1);
  TPad *c_summ0 = new TPad("c_summ0", "Summs", 0, 0, 1, 1);
  TPad *c_asym0 = new TPad("c_asym0", "Asyms", 0, 0, 1, 1);
  TPad *c_pol0 = new TPad("c_pol0", "Polarization 0", 0, 0, 1, 1);
  TPad *c_pol0_graphs = new TPad("c_pol0_graphs", "Polarization Graphs 0", 0, 0, 1, 1);
  TPad *c_pol4 = new TPad("c_pol4", "Polarization 4", 0, 0, 1, 1);
  TPad *c_pol4_graphs = new TPad("c_pol4_graphs", "Polarization Graphs 4", 0, 0, 1, 1);
  TPad *c_back_asyms = new TPad("c_back_asyms", "Background Detector Asyms", 0, 0, 1, 1);
  TPad *c_back_rates = new TPad("c_back_rates", "Background Detector Rates", 0, 0, 1, 1);

  essential_stats_pad(infile, run_num, output_path, c_ess_stats, true); //ess_stats.png
  essential_stats_2_pad(infile, run_num, output_path, c_ess_stats_2, true); //ess_stats_2.png
  snapshots_pad(infile, run_num, output_path, c_snap, true); //snapshots.png
  breakdown_pad(infile, run_num, output_path, c_sums, "triggerwise", "sums", true); //sums.png
  breakdown_pad(infile, run_num, output_path, c_acc0, "mpswise", "acc0", true); //acc0.png
  acc0_time_pad(infile, run_num, output_path, c_acc0_time, true); // acc0_time.png
  quartet_pad(infile, run_num, output_path, c_posH0, "posH0", true); //posH0.png
  quartet_pad(infile, run_num, output_path, c_negH0, "negH0", true); //negH0.png
  quartet_pad(infile, run_num, output_path, c_diff0, "diff0", true); //diff0.png
  quartet_pad(infile, run_num, output_path, c_summ0, "summ0", true); //summ0.png
  quartet_pad(infile, run_num, output_path, c_asym0, "asym0", true); //asym0.png
  asym_pad(infile, run_num, output_path, c_pol0, 0, true); //asym_acc0_hists.png
  asym_pad(infile, run_num, output_path, c_pol4, 4, true); //asym_acc4_hists.png
  asym_graph_pad(infile, run_num, output_path, c_pol0_graphs, 0, true); //q_acc0_graphs.png
  asym_graph_pad(infile, run_num, output_path, c_pol4_graphs, 4, true); //q_acc4_graphs.png
  detector_asyms(infile, run_num, output_path, c_back_asyms, true); //back_asyms.png
  detector_rates(infile, run_num, output_path, c_back_rates, true); //back_rates.png

  convert(run_num, output_path);
}
