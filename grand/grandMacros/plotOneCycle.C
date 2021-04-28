#include "../grandOnline/makePlots.h"
#include "../grandConstruction/buildGrandRootfile.h"

#include <vector>

using namespace std;

void fillGraph(TGraph *g, TTree *quartetwise, Int_t start, Int_t end){
  Double_t posAcc, posSamp, negAcc, negSamp;
  Int_t mpsNum, laserState, beamState;

  quartetwise->SetBranchAddress("firstMPSnumber", &mpsNum);
  quartetwise->SetBranchAddress("laserState", &laserState);
  quartetwise->SetBranchAddress("beamState", &beamState);
  quartetwise->SetBranchAddress("PosHelAcc0", &posAcc); quartetwise->SetBranchAddress("PosHelNSamples0", &posSamp);
  quartetwise->SetBranchAddress("NegHelAcc0", &negAcc); quartetwise->SetBranchAddress("NegHelNSamples0", &negSamp);

  Int_t nPts = 0;
  for(Int_t i = 0; i < quartetwise->GetEntries(); i++){
    quartetwise->GetEntry(i);
    if(mpsNum>=start && mpsNum<=end && laserState<4 && beamState==1){
      g->SetPoint(nPts++, mpsNum, posAcc/posSamp + negAcc/negSamp);
    }
  }
}

void plotOneCycle(Int_t runNum, Int_t cycNum){
  TFile *f = TFile::Open(Form("%s/Run%i_Plots.root", getenv("COMPMON_RUNPLOTS"), runNum));
  TFile *infile = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  TTree *quartetwise = (TTree *)infile->Get("quartetwise");
  
  TH1F *diffOn = (TH1F *)f->Get(Form("h%i.%i_DiffAcc0LasOn", runNum, cycNum));
  TH1F *diffOff1 = (TH1F *)f->Get(Form("h%i.%i_DiffAcc0LasOff1", runNum, cycNum));
  TH1F *diffOff2 = (TH1F *)f->Get(Form("h%i.%i_DiffAcc0LasOff2", runNum, cycNum));
  TH1F *summOn = (TH1F *)f->Get(Form("h%i.%i_SumAcc0LasOn", runNum, cycNum));
  TH1F *summOff1 = (TH1F *)f->Get(Form("h%i.%i_SumAcc0LasOff1", runNum, cycNum));
  TH1F *summOff2 = (TH1F *)f->Get(Form("h%i.%i_SumAcc0LasOff2", runNum, cycNum));
  TH1F *asymOn = (TH1F *)f->Get(Form("h%i.%i_AsymAcc0LasOn", runNum, cycNum));
  TH1F *asymOff1 = (TH1F *)f->Get(Form("h%i.%i_AsymAcc0LasOff1", runNum, cycNum));
  TH1F *asymOff2 = (TH1F *)f->Get(Form("h%i.%i_AsymAcc0LasOff2", runNum, cycNum));

  Int_t width = 2;
  diffOn->SetLineColor(kGreen + 2); diffOff1->SetLineColor(kRed); diffOff2->SetLineColor(kRed + 2);
  summOn->SetLineColor(kGreen + 2); summOff1->SetLineColor(kRed); summOff2->SetLineColor(kRed + 2);
  asymOn->SetLineColor(kGreen + 2); asymOff1->SetLineColor(kRed); asymOff2->SetLineColor(kRed + 2);
  diffOn->SetLineWidth(width); diffOff1->SetLineWidth(width); diffOff2->SetLineWidth(width);
  summOn->SetLineWidth(width); summOff1->SetLineWidth(width); summOff2->SetLineWidth(width);
  asymOn->SetLineWidth(width); asymOff1->SetLineWidth(width); asymOff2->SetLineWidth(width);

  THStack *hsDiff = new THStack("hsDiff", Form("Run %i, Cycle %i: Helicity Correlated Differences", runNum, cycNum));
  hsDiff->Add(diffOn); hsDiff->Add(diffOff1); hsDiff->Add(diffOff2);
  THStack *hsSumm = new THStack("hsSumm", Form("Run %i, Cycle %i: Helicity Pattern Sums", runNum, cycNum));
  hsSumm->Add(summOn); hsSumm->Add(summOff1); hsSumm->Add(summOff2);
  THStack *hsAsym = new THStack("hsAsym", Form("Run %i, Cycle %i: Helicity-Correlated Asymmetries", runNum, cycNum));
  hsAsym->Add(asymOn); hsAsym->Add(asymOff1); hsAsym->Add(asymOff2);
  
  Float_t xmin = 0.15; Float_t xmax = 0.4; Float_t ymin = 0.75; Float_t ymax = 0.9;
  TString lasOn("Laser ON"); TString lasOff1("Laser OFF 1"); TString lasOff2("Laser OFF 2");
  TLegend *lDiff = new TLegend(xmin, ymin, xmax, ymax);
  TLegend *lSumm = new TLegend(xmin, ymin, xmax, ymax);
  TLegend *lAsym = new TLegend(xmin, ymin, xmax, ymax);
  lDiff->AddEntry(diffOn, lasOn.Data()); lDiff->AddEntry(diffOff1, lasOff1.Data()); lDiff->AddEntry(diffOff2, lasOff2.Data());
  lSumm->AddEntry(summOn, lasOn.Data()); lSumm->AddEntry(summOff1, lasOff1.Data()); lSumm->AddEntry(summOff2, lasOff2.Data());
  lAsym->AddEntry(asymOn, lasOn.Data()); lAsym->AddEntry(asymOff1, lasOff1.Data()); lAsym->AddEntry(asymOff2, lasOff2.Data());

  vector<int> cycLimits = findCycles(runNum)[cycNum-1];

  TGraph *g = new TGraph(); g->SetName("gSummMPS");
  g->SetTitle(Form("Run %i, Cycle %i: Pattern Sums vs Time", runNum, cycNum));
  g->SetMarkerStyle(7);
  fillGraph(g, quartetwise, cycLimits[0], cycLimits[5]);

  TCanvas *c = new TCanvas("c", "Cycle Canavs", 1200, 1000);
  c->Divide(2, 2);
  c->cd(1); hsDiff->Draw("nostack"); lDiff->Draw("same");
  c->cd(2); hsSumm->Draw("nostack"); lSumm->Draw("same");
  c->cd(3); hsAsym->Draw("nostack"); lAsym->Draw("same");
  c->cd(4); g->Draw("ap");
}
