#ifndef COMPTON_H
#define COMPTON_H
#include <TChain.h>

void loadTree(int runnum, TChain *chain)
{
  TString filesPre = Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),runnum);
  chain->Add(filesPre+".root");
  chain->Add(filesPre+"_*.root");
}

#endif // COMPTON_H
