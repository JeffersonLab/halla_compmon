#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

#include "buildGrandRootfile.h"
#include "plot.h"
#include "burst.h"

using namespace std;


void buildRunRootfile(Int_t runNum){
  TFile *runOut = new TFile(Form("%s/Run%i_Plots.root", getenv("COMPMON_RUNPLOTS"), runNum), "RECREATE");
  //TFile *infile = new TFile(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum), "READ");
  vector<TChain *> runChains = loadChain(runNum);
  TChain *mpswise = runChains[0]; TChain *quartetwise = runChains[1];
  TChain *epicswise = runChains[4]; TChain *triggerwise = runChains[3];

  vector<vector<int>> cycles = findCycles(runNum);
  for(Int_t c = 0; c < cycles.size(); c++){
    printf("Plotting cycle %i/%i...\n", c + 1, (Int_t)cycles.size());
    printf("  Plotting mpswise vars...\n");
    cycMPSPlots(mpswise, runNum, c+1, cycles[c], runOut);
    printf("  Plotting quartetwise vars...\n");
    cycQrtPlots(quartetwise, runNum, c+1, cycles[c], runOut);
    printf("  Plotting triggerwise vars...\n");
    cycTrgPlots(triggerwise, runNum, c+1, cycles[c], runOut);
    printf("  Plotting burst vars...\n");
    cycBurstPlots(quartetwise, runNum, c+1, runOut);
  }
  printf("  Plotting run mps vars...\n");
  runMPSPlots(mpswise, runNum, runOut);
  runEpicsPlots(epicswise, mpswise, runNum, runOut);

  //runOut->Write(); 
  runOut->Close();
  printf("...Done!\n");
}
