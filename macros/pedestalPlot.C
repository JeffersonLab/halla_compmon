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

vector<vector<Int_t>> findCycles(Int_t runNum){
  ifstream infile(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), runNum));
  string readStr;
  Int_t cyclesInThisRun = 0;
  vector<vector<Int_t>> run;
  while(getline(infile, readStr)){
    vector<Int_t> limits; stringstream ss(readStr);
    for(Int_t i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    cyclesInThisRun++;
    vector<Int_t> cycle;
    for(Int_t i = 0; i < limits.size(); i++)
      cycle.push_back(limits[i]);
    run.push_back(cycle);
  }
  printf("Looked in %s/cycles_%i.dat\n", getenv("COMPMON_MINIRUNS"), runNum);
  printf("Found %i cycles\n", (Int_t)run.size());
  return run;
}

void pedestalPlot(Int_t runNum){
  TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum));
  vector<vector<Int_t>> cycles = findCycles(runNum);
  TTree *triggerwise = (TTree *)f->Get("triggerwise");

  vector<Int_t> accMPSOff; vector<Float_t> pedsOff; vector<Float_t> pedErrsOff;
  vector<Int_t> accMPSOn; vector<Float_t> pedsOn; vector<Float_t> pedErrsOn;
  TString correctPedCut1("2.0*abs(sumPre - sumPost)/(sumPre + sumPost) < 0.03");
  TString correctPedCut2("(sumPre + sumPost)/2.0 < 3900");
  TString narrowPedCut("abs(sumPre - sumPost) < 50");

  TCanvas *cFitOn = new TCanvas("cFitOn", "Pedestal Fit On", 1200, 800);
  TCanvas *cFitOff = new TCanvas("cFitOff", "Pedestal Fit Off", 1200, 800);
  
  for(Int_t i = 0; i < cycles.size(); i++){
    TString hNameOn = Form("hPedOn_%i", i);    
    TString hNameOff = Form("hPedOff_%i", i);
    TString fNameOn1 = Form("fPedOn1_%i", i); TString fNameOn2 = Form("fPedOn2_%i", i);
    TString fNameOff1 = Form("fPedOff1_%i", i); TString fNameOff2 = Form("fPedOff2_%i", i);
    TH1F *hPedOn = new TH1F(hNameOn.Data(), "Pedestal Measurement", 250, 3770, 3800);
    hPedOn->GetXaxis()->SetTitle("Pedestal (RAU)");
    TH1F *hPedOff = new TH1F(hNameOff.Data(), "Pedestal Measurement", 250, 3770, 3800);
    hPedOff->GetXaxis()->SetTitle("Pedestal (RAU)");
    TF1 *fPedOn1 = new TF1(fNameOn1.Data(), "gaus");
    TF1 *fPedOn2 = new TF1(fNameOn2.Data(), "gaus");
    TF1 *fPedOff1 = new TF1(fNameOff1.Data(), "gaus");
    TF1 *fPedOff2 = new TF1(fNameOff2.Data(), "gaus");
    
    //TString onCut = Form("(mpsCoda>=%i && mpsCoda<=%i) && (laserState==0 || laserState==1)", 
    //                      cycles[i][2], cycles[i][3]);
    //TString offCut = Form("((mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i)) && (laserState==2 || laserState==3)", 
    //                      cycles[i][0], cycles[i][1], cycles[i][4], cycles[i][5]);
    TString onCut = Form("((mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i)) && (laserState==2 || laserState==3) && sumIsRandom==0",
                          cycles[i][0], cycles[i][1], cycles[i][4], cycles[i][5]);
    TString offCut = Form("((mpsCoda>=%i && mpsCoda<=%i) || (mpsCoda>=%i && mpsCoda<=%i)) && (laserState==2 || laserState==3) && sumIsRandom==1",
                          cycles[i][0], cycles[i][1], cycles[i][4], cycles[i][5]);
    TString allCutsOn = Form("%s && %s && %s && %s && beamState==1", onCut.Data(), correctPedCut1.Data(),
                              correctPedCut2.Data(), narrowPedCut.Data());
    TString allCutsOff = Form("%s && %s && %s && %s && beamState==1", offCut.Data(), correctPedCut1.Data(),
                              correctPedCut2.Data(), narrowPedCut.Data());    
    //TString allCuts = Form("%s && %s && beamState==1", mpsCut.Data(), correctPedCut1.Data());
    triggerwise->Project(hNameOff.Data(), "(sumPre + sumPost)/2.0", allCutsOff.Data());
    triggerwise->Project(hNameOn.Data(), "(sumPre + sumPost)/2.0", allCutsOn.Data());

    cFitOn->cd();
    hPedOn->Draw();
    Float_t startMeanOn1 = hPedOn->GetMean(); Float_t startRMSOn1 = hPedOn->GetRMS();
    hPedOn->Fit(fNameOn1.Data(), "", "goff", startMeanOn1 - startRMSOn1, startMeanOn1 + startRMSOn1);
    Float_t startMeanOn2 = fPedOn1->GetParameter(1); Float_t startRMSOn2 = fPedOn1->GetParameter(2);
    hPedOn->Fit(fNameOn2.Data(), "", "", startMeanOn2 - startRMSOn2, startMeanOn2 + startRMSOn2);
    //accMPSOn.push_back((cycles[i][2] + cycles[i][3])/2);
    accMPSOn.push_back((cycles[i][0] + cycles[i][1] + cycles[i][4] + cycles[i][5])/4);
    pedsOn.push_back(fPedOn2->GetParameter(1)); pedErrsOn.push_back(fPedOn2->GetParError(1));

    cFitOff->cd();
    hPedOff->Draw();
    Float_t startMeanOff1 = hPedOff->GetMean(); Float_t startRMSOff1 = hPedOff->GetRMS();
    hPedOff->Fit(fNameOff1.Data(), "", "goff", startMeanOff1 - startRMSOff1, startMeanOff1 + startRMSOff1);
    Float_t startMeanOff2 = fPedOff1->GetParameter(1); Float_t startRMSOff2 = fPedOff1->GetParameter(2);
    hPedOff->Fit(fNameOff2.Data(), "", "", startMeanOff2 - startRMSOff2, startMeanOff2 + startRMSOff2);
    accMPSOff.push_back((cycles[i][0] + cycles[i][1] + cycles[i][4] + cycles[i][5])/4);
    pedsOff.push_back(fPedOff2->GetParameter(1)); pedErrsOff.push_back(fPedOff2->GetParError(1));

    //cFitOn->SaveAs(Form("FitOn_%i.pdf", i)); cFitOff->SaveAs(Form("FitOff_%i.pdf", i));
  }

  TH2F *h2Ped = new TH2F("h2Ped", "Calculated Pedestal Comparison", 
                          400, 0, cycles[(Int_t)cycles.size() - 1][5], 200, 3770, 3800);
  h2Ped->GetXaxis()->SetTitle("mpsCoda"); h2Ped->GetYaxis()->SetTitle("Pedestal (RAU)");
  triggerwise->Project("h2Ped", "(sumPre + sumPost)/2.0:mpsCoda", Form("%s && %s && %s && beamState==1",
                        correctPedCut1.Data(), correctPedCut2.Data(), narrowPedCut.Data()), "goff");

  TCanvas *cGraph = new TCanvas("cGraph", "Pedestal Graph", 1200, 800);
  TCanvas *cComp = new TCanvas("cComp", "Pedestal Calc and Graph Comparison", 1200, 800);
  gStyle->SetOptStat(0);
  
  TGraphErrors *gPedOn = new TGraphErrors(); gPedOn->SetName("gPedOn");
  TGraphErrors *gPedOff = new TGraphErrors(); gPedOff->SetName("gPedOff");
  gPedOn->SetTitle(Form("Run %i: Calculated Pedestal by MPS", runNum)); 
  gPedOn->SetMarkerStyle(21); gPedOn->SetMarkerColor(3);
  gPedOn->GetXaxis()->SetTitle("mpsCoda"); gPedOn->GetYaxis()->SetTitle("Pedestal (RAU)");
  gPedOff->SetTitle(Form("Run %i: Calculated Pedestal by MPS", runNum)); 
  gPedOff->SetMarkerStyle(21); gPedOff->SetMarkerColor(2);
  gPedOff->GetXaxis()->SetTitle("mpsCoda"); gPedOff->GetYaxis()->SetTitle("Pedestal (RAU)");
  TGraph *gPedLineOn = new TGraph(); gPedLineOn->SetName("gPedLineOn");
  gPedLineOn->SetTitle(Form("Run %i: Calculated Pedestal Comparison", runNum));
  gPedLineOn->SetMarkerStyle(21); gPedLineOn->SetMarkerColor(3);
  gPedLineOn->SetLineColor(3); gPedLineOn->SetLineWidth(4);
  TGraph *gPedLineOff = new TGraph(); gPedLineOff->SetName("gPedLineOff");
  gPedLineOff->SetTitle(Form("Run %i: Calculated Pedestal Comparison", runNum));
  gPedLineOff->SetMarkerStyle(21); gPedLineOff->SetMarkerColor(2);
  gPedLineOff->SetLineColor(2); gPedLineOff->SetLineWidth(4);

  for(Int_t i = 0; i < pedsOn.size(); i++){
    gPedOn->SetPoint(i, accMPSOn[i], pedsOn[i]);
    gPedOn->SetPointError(i, 0, pedErrsOn[i]);
    gPedLineOn->SetPoint(i, accMPSOn[i], pedsOn[i]);
  }
  for(Int_t i = 0; i < pedsOn.size(); i++){
    gPedOff->SetPoint(i, accMPSOff[i], pedsOff[i]);
    gPedOff->SetPointError(i, 0, pedErrsOff[i]);
    gPedLineOff->SetPoint(i, accMPSOff[i], pedsOff[i]);
  }

  cGraph->cd(); gPedOn->Draw("ap");
  gPedOff->Draw("p && same");
  cComp->cd(); h2Ped->Draw("colz");
  gPedLineOn->Draw("L && same");
  gPedLineOff->Draw("L && same");
}
