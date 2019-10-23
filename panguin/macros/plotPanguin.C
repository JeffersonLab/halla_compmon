#include <TCanvas.h>
#include <TPad.h>
#include <TObject.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TPad.h>
#include <TString.h>
#include <THStack.h>
#include <vector>

#include "../../utils.h"

using namespace std;

void plotPanguin(int run_num, int id){
  //TTree *runwise = (TTree *)gDirectory->Get("runwise");
  //Int_t run_num_base; runwise->SetBranchAddress("runNumber", &run_num_base);
  //runwise->GetEntry(0); int run_num = run_num_base;
  TString output_path(Form("%s/runs/Run%i/", getenv("COMPMON_WEB"), run_num));
  TString compmon_out(getenv("COMPMON_PLOTFILES"));
  TFile *f = TFile::Open(Form("%s/compton_online_run_%i.root", compmon_out.Data(), run_num));
  TPad *myPad = new TPad("myPad", "Acc0 Pad", 0, 0, 1, 1); myPad->Draw(); myPad->cd();
  gStyle->SetStatH(0.8); gStyle->SetStatW(0.2);

  if(id == 1){essential_stats_pad(f, run_num, output_path, myPad);}
  else if(id == 2){ snapshots_pad(f, run_num, output_path, myPad);}
  else if(id == 3) {breakdown_pad(f, run_num, output_path, myPad, "triggerwise", "sums");}
  else if(id == 4) {breakdown_pad(f, run_num, output_path, myPad, "mpswise", "acc0");}
  else if(id == 5) {acc0_time_pad(f, run_num, output_path, myPad);}
  else if(id == 6) {quartet_pad(f, run_num, output_path, myPad, "posH0");}
  else if(id == 7) {quartet_pad(f, run_num, output_path, myPad, "negH0");}
  else if(id == 8) {quartet_pad(f, run_num, output_path, myPad, "diff0");}
  else if(id == 9) {quartet_pad(f, run_num, output_path, myPad, "summ0");}
  else if(id == 10){quartet_pad(f, run_num, output_path, myPad, "asym0");}
  else if(id == 11){asym_pad(f, run_num, output_path, myPad);}
  else if(id == 12){asym_graph_pad(f, run_num, output_path, myPad);}
  else if(id == 13){asym_pad(f, run_num, output_path, myPad, 4);}
  else if(id == 14){asym_graph_pad(f, run_num, output_path, myPad, 4);}
  else if(id == 15){beam_off_quartet_pad(f, run_num, output_path, myPad);}
  else{printf("Nothing to plot.");}
}
