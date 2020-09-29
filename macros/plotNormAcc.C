void plotNormAcc(){
  const Int_t nRuns = 24;
  Int_t runs[nRuns] =  {6086, 6087, 6088, 6089, 6090, 6091, 6092, 6093,
                        6094, 6096, 6098, 6099, 6100, 6101, 6102, 6105,
                        6108, 6110, 6111, 6112, 6113, 6114, 6115, 6116};
  TGraphErrors *gON = new TGraphErrors();
  TGraphErrors *gOFF = new TGraphErrors();

  for(Int_t i = 0; i < nRuns; i++){
    TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runs[i]));
    TTree *mpswise = (TTree *)f->Get("mpswise");
    TString hNameON = Form("hON_run%i", runs[i]);
    TString hNameOFF = Form("hOFF_run%i", runs[i]);
    mpswise->Draw(Form("(Acc0/NAcc0)/bcm>>%s", hNameON.Data()), "bcm>140 && (laserState==0 || laserState==1)", "goff");
    mpswise->Draw(Form("(Acc0/NAcc0)/bcm>>%s", hNameOFF.Data()), "bcm>140 && (laserState==2 || laserState==3)", "goff");
    //mpswise->Draw(Form("(Acc0/NAcc0)>>%s", hNameON.Data()), "bcm>140 && (laserState==0 || laserState==1)", "goff");
    //mpswise->Draw(Form("(Acc0/NAcc0)>>%s", hNameOFF.Data()), "bcm>140 && (laserState==2 || laserState==3)", "goff");
    TH1F *hON = (TH1F *)gDirectory->Get(hNameON.Data());
    TH1F *hOFF = (TH1F *)gDirectory->Get(hNameOFF.Data());
    gON->SetPoint(i, runs[i], hON->GetMean()); gON->SetPointError(i, 0.0, hON->GetMeanError());
    gOFF->SetPoint(i, runs[i], hOFF->GetMean()); gOFF->SetPointError(i, 0.0, hOFF->GetMeanError());
  }

  TCanvas *c = new TCanvas("c", "Acc0 Noramlized", 1600, 600);
  c->Divide(2, 1);
  c->cd(1);
  gON->SetTitle("Compton Laser ON: (Acc0/NAcc0)/bcm");
  //gON->SetTitle("Compton Laser ON: (Acc0/NAcc0)");
  gON->GetXaxis()->SetTitle("runNum"); gON->GetYaxis()->SetTitle("(Acc0/NAcc0)/bcm");
  //gON->GetXaxis()->SetTitle("runNum"); gON->GetYaxis()->SetTitle("(Acc0/NAcc0)");
  gON->Draw("ap");
  c->cd(2);
  gOFF->SetTitle("Compton Laser OFF: (Acc0/NAcc0)/bcm");
  //gOFF->SetTitle("Compton Laser OFF: (Acc0/NAcc0)");
  gOFF->GetXaxis()->SetTitle("runNum"); gOFF->GetYaxis()->SetTitle("(Acc0/NAcc0)/bcm");
  //gOFF->GetXaxis()->SetTitle("runNum"); gOFF->GetYaxis()->SetTitle("(Acc0/NAcc0)");
  gOFF->Draw("ap");
}
