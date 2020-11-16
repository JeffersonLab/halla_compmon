#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "../vars.h"

using namespace std;

void printPedInfo(Int_t selRunNum, Int_t selCycNum){
  TFile *f = TFile::Open(Form("%s/prexGrandCompton.root", getenv("COMPMON_GRAND")));
  TTree *cyc = (TTree *)f->Get("cyc");
  
  Int_t runNum, cycNum;
  DataVar fOff, lOff;
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("PedestalMeanFirstOff", &fOff);
  cyc->SetBranchAddress("PedestalMeanLastOff", &lOff);

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    if(runNum == selRunNum && cycNum == selCycNum){
      printf("Run %i, Cycle%i:\n", runNum, cycNum);
      printf("  First Off Mean: %.4f\n", fOff.mean);
      printf("  First Off Err: %.4f\n", fOff.meanErr);
      printf("  First Off RMS: %.4f\n", fOff.rms);
      printf("  Last Off Mean: %.4f\n", lOff.mean);
      printf("  Last Off Err: %.4f\n", lOff.meanErr);
      printf("  Last Off RMS: %.4f\n", lOff.rms);
    }
  }
}
