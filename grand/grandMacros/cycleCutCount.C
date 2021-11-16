#include "../grandOnline/makePlots.h"

using namespace std;


const Int_t nCuts = 8;
Int_t cutCounts[nCuts] = {0, 0, 0, 0, 0, 0, 0, 0};


void cycleCutCount(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TFile *f = TFile::Open(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()));
  TTree *cyc = (TTree *)f->Get("cyc");

  Int_t runNum, cycNum, cycleCut;
  
  cyc->SetBranchAddress("runNum", &runNum);
  cyc->SetBranchAddress("cycleNum", &cycNum);
  cyc->SetBranchAddress("CycleCut", &cycleCut);

  for(Int_t i = 0; i < cyc->GetEntries(); i++){
    cyc->GetEntry(i);
    for(Int_t j = 0; j < nCuts; j++){
      Int_t cutSel = (Int_t)TMath::Power(2, j + 1) - 1;
      Int_t check = (Int_t)TMath::Power(2, j);
      if((cycleCut&cutSel)==check){
        cutCounts[j]++;
        printf("Found a cycle that hasn't been cut yet: %i.%i with CycleCut=%i and cutSel=%i\n", runNum, cycNum, cycleCut, cutSel);
        break;
      }
    }
  }

  printf("\n\n");

  for(Int_t i = 0; i < nCuts; i++){
    printf("Cut #%i has %i cycles\n", i+1, cutCounts[i]);
  }
}
