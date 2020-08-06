#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"

#include <vector>

using namespace std;

//Before mighty Facebook, how would I have found Tom Hanks?
void beamStatePedestalComp(){
  Int_t runNum = 4537;
  Int_t mpsStart = 1164e3;
  Int_t mpsEnd = 1178e3;
  Int_t mpsMiddle = (mpsStart + mpsEnd)/2;

  TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *triggerwise = (TTree *)f->Get("triggerwise");
  TString correctPedCut1("2.0*abs(sumPre - sumPost)/(sumPre + sumPost) < 0.03");
  TString correctPedCut2("(sumPre + sumPost)/2.0 < 3900");
  TString narrowPedCut("abs(sumPre - sumPost) < 50");
  TString mpsCut = Form("mpsCoda>=%i && mpsCoda<=%i", mpsStart, mpsEnd);
  TString firstMPSCut = Form("mpsCoda>=%i && mpsCoda<=%i", mpsStart, mpsMiddle);
  TString lastMPSCut = Form("mpsCoda>=%i && mpsCoda<=%i", mpsMiddle, mpsEnd);
  TString beamOnCut("beamState==1");
  TString beamOffCut("beamState==0");
  TString laserCut("(laserState==0 || laserState==1)");

  TString pedCuts = Form("%s && %s && %s", correctPedCut1.Data(), correctPedCut2.Data(), narrowPedCut.Data());
  TString firstSectionCut = Form("%s && %s && %s", firstMPSCut.Data(), beamOnCut.Data(), laserCut.Data());
  TString secondSectionCut = Form("%s && %s && %s", mpsCut.Data(), beamOffCut.Data(), laserCut.Data());
  TString thirdSectionCut = Form("%s && %s && %s", lastMPSCut.Data(), beamOnCut.Data(), laserCut.Data());

  TH1F *hFirst = new TH1F("hFirst", "First Beam On Section", 250, 3770, 3800);
  TH1F *hSecond = new TH1F("hSecond", "Beam Off Section", 250, 3770, 3800);
  TH1F *hThird = new TH1F("hThird", "Second Beam On Section", 250, 3770, 3800);
  TF1 *fFirst1 = new TF1("fFirst1", "gaus"); TF1 *fFirst2 = new TF1("fFirst2", "gaus");
  TF1 *fSecond1 = new TF1("fSecond1", "gaus"); TF1 *fSecond2 = new TF1("fSecond2", "gaus");
  TF1 *fThird1 = new TF1("fThird1", "gaus"); TF1 *fThird2 = new TF1("fThird2", "gaus");

  TCanvas *cFirstSec = new TCanvas("cFirstSec", "First Section Fit", 1200, 800);
  TCanvas *cSecondSec = new TCanvas("cSecondSec", "Second Section Fit", 1200, 800);
  TCanvas *cThirdSec = new TCanvas("cThirdSec", "Third Section Fit", 1200, 800);

  cFirstSec->cd();
  triggerwise->Project("hFirst", "(sumPre + sumPost)/2.0", Form("%s && %s", pedCuts.Data(), firstSectionCut.Data()));
  Float_t firstMean1 = hFirst->GetMean(); Float_t firstRMS1 = hFirst->GetRMS();
  hFirst->Draw();
  hFirst->Fit("fFirst1", "", "goff", firstMean1 - firstRMS1, firstMean1 + firstRMS1);
  Float_t firstMean2 = fFirst1->GetParameter(1); Float_t firstRMS2 = fFirst1->GetParameter(2);
  printf("Mean and RMS: %.2f, %.2f\n", firstMean2, firstRMS2);
  hFirst->Fit("fFirst2", "", "", firstMean2 - firstRMS2, firstMean2 + firstRMS2);

  cSecondSec->cd();
  triggerwise->Project("hSecond", "(sumPre + sumPost)/2.0", Form("%s && %s", pedCuts.Data(), secondSectionCut.Data()));
  Float_t secondMean1 = hSecond->GetMean(); Float_t secondRMS1 = hSecond->GetRMS();
  hSecond->Draw();
  hSecond->Fit("fSecond1", "", "goff", secondMean1 - secondRMS1, secondMean1 + secondRMS1);
  Float_t secondMean2 = fSecond1->GetParameter(1); Float_t secondRMS2 = fSecond1->GetParameter(2);
  printf("Mean and RMS: %.2f, %.2f\n", secondMean2, secondRMS2);
  hSecond->Fit("fSecond2", "", "", secondMean2 - secondRMS2, secondMean2 + secondRMS2);


  cThirdSec->cd();
  triggerwise->Project("hThird", "(sumPre + sumPost)/2.0", Form("%s && %s", pedCuts.Data(), thirdSectionCut.Data()));
  Float_t thirdMean1 = hThird->GetMean(); Float_t thirdRMS1 = hThird->GetRMS();
  hThird->Draw();
  hThird->Fit("fThird1", "", "goff", thirdMean1 - thirdRMS1, thirdMean1 + thirdRMS1);
  Float_t thirdMean2 = fThird1->GetParameter(1); Float_t thirdRMS2 = fThird1->GetParameter(2);
  printf("Mean and RMS: %.2f, %.2f\n", thirdMean2, thirdRMS2);
  hThird->Fit("fThird2", "", "", thirdMean2 - thirdRMS2, thirdMean2 + thirdRMS2);
  
  printf("First fit mean and error: %.2f +/- %.2f\n", fFirst2->GetParameter(1), fFirst2->GetParError(1));
  printf("Second fit mean and error: %.2f +/- %.2f\n", fSecond2->GetParameter(1), fSecond2->GetParError(1));
  printf("Third fit mean and error: %.2f +/- %.2f\n", fThird2->GetParameter(1), fThird2->GetParError(1));
}
