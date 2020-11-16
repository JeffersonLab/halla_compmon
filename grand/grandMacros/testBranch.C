#include "makePlots.h"
#include "../vars.h"

void testBranch(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/%sGrandCompton.root", getenv("COMPMON_GRAND"), expt.Data()));
  TTree *run = (TTree *)f->Get("run");
  
  Int_t runNum;
  FitPolVar test;
  run->SetBranchAddress("runNum", &runNum);
  run->SetBranchAddress("Asym0LasOff", &test);

  for(Int_t i = 0; i < run->GetEntries(); i++){
    run->GetEntry(i);
    printf("Run %i - Asym: %0.2f, Err: %0.2f, Chi2: %f, NDF: %i\n", runNum, 
            1000*test.mean, 1000*test.meanErr, test.Chi2, test.NDF);
  }
}
