#include "../online/utils.h"

using namespace std;

vector<vector<int>> findCycles(int runNum){
  ifstream infile(Form("%s/cycles_%i.dat", getenv("COMPMON_CYCLES"), runNum));
  string readStr;
  int cyclesInThisRun = 0;
  vector<vector<int>> run;
  while(getline(infile, readStr)){
    vector<int> limits; stringstream ss(readStr);
    for(int i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
    cyclesInThisRun++;
    vector<int> cycle;
    for(int i = 0; i < limits.size(); i++)
      cycle.push_back(limits[i]);
    run.push_back(cycle);
  }
  printf("Looked in %s/cycles_%i.dat\n", getenv("COMPMON_CYCLES"), runNum);
  printf("Found %i cycles\n", (Int_t)run.size());
  return run;
}

void whatsThePolarization(Int_t runNum){
  TCanvas *cPol = new TCanvas("cPol", "Polarization Canvas", 1200, 800);
  TPad *pPol = new TPad("pPol", "Pol Avg", 0.0, 0.0, 1.0, 1.0);
  cPol->cd(); pPol->Draw();

  TFile *f = new TFile(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runNum), "READ");
  if(!f->IsOpen()){
    TPaveText *ptErr = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
    ptErr->AddText("Can't display polarization.\n");
    ptErr->AddText(Form("Rootfile for run %i doesn't exist.\n", runNum));
    pPol->cd(); ptErr->Draw();
    return;
  }
  TTree *quartetwise = (TTree *)f->Get("quartetwise");

  vector<vector<int>> cycles = findCycles(runNum);
  TH1F *hPol = new TH1F("hPol", Form("Run %i: Polarization Avg (Unofficial)", runNum), (Int_t)cycles.size(), 0, (Int_t)cycles.size());
  hPol->GetXaxis()->SetTitle("cycleNum"); hPol->GetYaxis()->SetTitle("Pol0 (pct)");
  hPol->SetStats(0);

  if(cycles.size() < 5){
    TPaveText *ptErr = new TPaveText(0.0, 0.0, 1.0, 1.0, "blNDC");
    ptErr->AddText("Can't display polarization.\n");
    ptErr->AddText("Too few laser cycles. (Need >= 5.)\n");
    ptErr->AddText(Form("Found %i laser cycles.\n", (Int_t)cycles.size()));
    pPol->cd(); ptErr->Draw();
    return;
  }

  TString posH("PosHelAcc0"); TString negH("NegHelAcc0");
  TString diff = Form("(%s - %s)", posH.Data(), negH.Data());
  TString summ = Form("(%s + %s)", posH.Data(), negH.Data());

  TString B1L1("beamState==1 && (laserState==0 || laserState==1)");
  TString B1L0("beamState==1 && (laserState==2 || laserState==3)");

  for(Int_t i = 0; i < cycles.size(); i++){
    printf("Plotting cycle %i...\n", i+1);
    TString hDiffOnName = Form("hDiff_On_%i", i+1);
    TString hDiffOff1Name = Form("hDiff_Off1_%i", i+1);
    TString hDiffOff2Name = Form("hDiff_Off2_%i", i+1);
    TString hSummOnName = Form("hSumm_On_%i", i+1);
    TString hSummOff1Name = Form("hSumm_Off1_%i", i+1);
    TString hSummOff2Name = Form("hSumm_Off2_%i", i+1);

    TString onCut = Form("firstMPSnumber>=%i && firstMPSnumber<=%i", cycles[i][2], cycles[i][3]);
    TString off1Cut = Form("firstMPSnumber>=%i && firstMPSnumber<=%i", cycles[i][0], cycles[i][1]);
    TString off2Cut = Form("firstMPSnumber>=%i && firstMPSnumber<=%i", cycles[i][4], cycles[i][5]);

    quartetwise->Draw(Form("%s>>%s", diff.Data(), hDiffOnName.Data()), Form("%s && %s", B1L1.Data(), onCut.Data()), "goff");
    quartetwise->Draw(Form("%s>>%s", diff.Data(), hDiffOff1Name.Data()), Form("%s && %s", B1L0.Data(), off1Cut.Data()), "goff");
    quartetwise->Draw(Form("%s>>%s", diff.Data(), hDiffOff2Name.Data()), Form("%s && %s", B1L0.Data(), off2Cut.Data()), "goff");
    quartetwise->Draw(Form("%s>>%s", summ.Data(), hSummOnName.Data()), Form("%s && %s", B1L1.Data(), onCut.Data()), "goff");
    quartetwise->Draw(Form("%s>>%s", summ.Data(), hSummOff1Name.Data()), Form("%s && %s", B1L0.Data(), off1Cut.Data()), "goff");
    quartetwise->Draw(Form("%s>>%s", summ.Data(), hSummOff2Name.Data()), Form("%s && %s", B1L0.Data(), off2Cut.Data()), "goff");

    TH1F *hDiffOn = (TH1F *)gDirectory->Get(hDiffOnName.Data());
    TH1F *hDiffOff1 = (TH1F *)gDirectory->Get(hDiffOff1Name.Data());
    TH1F *hDiffOff2 = (TH1F *)gDirectory->Get(hDiffOff2Name.Data());
    TH1F *hSummOn = (TH1F *)gDirectory->Get(hSummOnName.Data());
    TH1F *hSummOff1 = (TH1F *)gDirectory->Get(hSummOff1Name.Data());
    TH1F *hSummOff2 = (TH1F *)gDirectory->Get(hSummOff2Name.Data());

    Float_t vDiffOn = hDiffOn->GetMean(); 
    Float_t vSummOn = hSummOn->GetMean();
    Float_t vErrDiffOn = hDiffOn->GetMeanError();
    Float_t vErrSummOn = hSummOn->GetMeanError();
    Float_t vDiffOff = (hDiffOff1->GetMean() + hDiffOff2->GetMean())/2.0; 
    Float_t vSummOff = (hSummOff1->GetMean() + hSummOff2->GetMean())/2.0;
    Float_t vErrDiffOff = TMath::Sqrt(TMath::Power(hDiffOff1->GetMeanError(), 2) + TMath::Power(hDiffOff2->GetMeanError(), 2))/2.0;
    Float_t vErrSummOff = TMath::Sqrt(TMath::Power(hSummOff1->GetMeanError(), 2) + TMath::Power(hSummOff2->GetMeanError(), 2))/2.0;
    Float_t vBkSubSumm = vSummOn - vSummOff;
    Float_t vErrBkSubSumm = TMath::Sqrt(TMath::Power(vErrSummOn, 2) + TMath::Power(vErrSummOff, 2));

    Float_t vAsymOn = vDiffOn/vBkSubSumm;
    Float_t vErrAsymOn = vAsymOn*TMath::Sqrt(TMath::Power(vErrDiffOn/vDiffOn, 2) + TMath::Power(vErrBkSubSumm/vBkSubSumm, 2));
    Float_t vAsymOff = vDiffOff/vBkSubSumm;
    Float_t vErrAsymOff = vAsymOff*TMath::Sqrt(TMath::Power(vErrDiffOff/vDiffOff, 2) + TMath::Power(vErrBkSubSumm/vBkSubSumm, 2));
    Float_t vBkSubAsym = vAsymOn - vAsymOff;
    Float_t vErrBkSubAsym = TMath::Sqrt(TMath::Power(vErrAsymOn, 2) + TMath::Power(vErrAsymOff, 2));
    Float_t pol = vBkSubAsym/get_analyzing_power(runNum);
    Float_t polErr = vErrBkSubAsym/get_analyzing_power(runNum);

    hPol->SetBinContent(i+1, 100*pol);
    hPol->SetBinError(i+1, 100*polErr);
    if(i == 0 || (i+1) % 5 == 0){
      hPol->GetXaxis()->SetBinLabel(i+1, Form("%i", i+1));
    }
  }

  pPol->SetGridx(1); pPol->SetGridy(1);
  pPol->cd(); hPol->Draw();

  TF1 *polFit = new TF1("polFit", "pol0");
  hPol->Fit("polFit", "Q", "", 0, (Int_t)cycles.size());

  TPaveText *pt0 = new TPaveText(0.7, 0.75, 0.98, 0.92, "blNDC");

  Float_t par = polFit->GetParameter(0);
  Float_t parErr = polFit->GetParError(0);
  Float_t chi2 = polFit->GetChisquare();
  Int_t ndf = polFit->GetNDF();

  pt0->AddText(Form("--------Fit Results--------"));
  pt0->AddText(Form("Mean +/- Err: %.3f +/- %.3f", par, parErr));
  pt0->AddText(Form("Rel. Error: %.3f%%", 100.0*parErr/par));
  pt0->AddText(Form("Chi^2 / ndf: %f / %d", chi2, ndf));
  pt0->SetBorderSize(1); pt0->SetFillColor(0);
  pt0->Draw("same");

  cPol->Print(Form("%s/macros/plots/polarization_run%i.pdf", getenv("COMPMON_DIR"), runNum), "pdf");
}
