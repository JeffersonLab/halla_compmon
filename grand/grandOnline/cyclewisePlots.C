#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "THStack.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "../grandConstruction/buildGrandRootfile.h"
#include "makePlots.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

using namespace std;

void cyclewisePlots(Int_t prexOrCrex){
  TString expt = experimentCode(prexOrCrex);
  TString fname = Form("%sGrandCompton.root", expt.Data());
  const Int_t nPolVars = 6;
  const Int_t nAsymVars = 11;
  const Int_t nDataVars = 24;
  const Int_t nOthVars = 8;
  TString polVars[nPolVars] = {"Pol0", "Asym0", "Asym0NGC", "Asym0LasOn", "Asym0LasOff", "SigSubSum0"};
  TString asymVars[nAsymVars] = {"AsymCentralRateLasOn", "AsymCentralRateLasOff", "AsymHFingerRateLasOn", "AsymHFingerRateLasOff",
                                 "AsymVFingerRateLasOn", "AsymVFingerRateLasOff", "AsymUSbg1", "AsymUSbg2",
                                 "AsymDSbg1", "AsymDSbg2", "AsymBCM"};
  TString dataVars[nDataVars] = {"Acc0LasOn", "Acc0LasOff1", "Acc0LasOff2", "Acc0BeamOff", 
                                 "Acc4LasOn", "Acc4LasOff1", "Acc4LasOff2", "Acc4BeamOff",
                                 "LaserPower", "BeamCurrent", "bpmAx", "bpmAy", "bpmBx", "bpmBy", 
                                 "CentralRateLasOn", "CentralRateLasOff", "USbg1", "USbg2", "DSbg1", "DSbg2",
                                 "HFingerLasOn", "HFingerLasOff", "VFingerLasOn", "VFingerLasOff"};
  //TString othVars[nOthVars] = {"runTime", "numCycles", "tablePosX", "tablePosY",
  //                             "qw1", "hw1", "qw2", "targetPos"};
  //Float_t polYmins[nPolVars] = {0.99, 0.99, 0.99, 0.9, 0.9, 0.9};
  //Float_t polYmaxs[nPolVars] = {1.01, 1.01, 1.01, 1.1, 1.1, 1.1};
  //Bool_t polSign[nPolVars] = {true, true, true, false, false, false};
  //Bool_t othFloat[nOthVars] = {true, false, true, true, true, true, true, true};

  Int_t msmtNum = 0;
  for(Int_t i = 0; i < nPolVars; i++){
    plotPolCyc(fname.Data(), polVars[i].Data(), msmtNum++, 0.99, 1.01, true);
    plotPolCyc(fname.Data(), polVars[i].Data(), msmtNum++, 0.90, 1.10, false);
  }
  for(Int_t i = 0; i < nAsymVars; i++){
    plotStandardCyc(fname.Data(), asymVars[i].Data(), msmtNum++, 0.9, 1.1);
  }
  for(Int_t i = 0; i < nDataVars; i++){
    plotStandardCyc(fname.Data(), dataVars[i].Data(), msmtNum++, 0.9, 1.1);
  }
  //for(Int_t i = 0; i < nOthVars; i++){
  //  plotStandardRun(fname.Data(), othVars[i].Data(), msmtNum++, 0.9, 1.1, othFloat[i]);
  //}

  gSystem->Exec(Form("pdfunite %s/grandOnline/plots/msmt*.pdf %s/aggregates/%sGrandCyclewise.pdf", 
                getenv("COMPMON_GRAND"), getenv("COMPMON_WEB"), expt.Data()));
  gSystem->Exec(Form("rm -f %s/grandOnline/plots/msmt*.pdf", getenv("COMPMON_GRAND")));
}
