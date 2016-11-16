// Description: Script that compares beam on/off, laser on/off states of various accumulators
// Author: Larisa Thorne
// Created: 2016-02-10
// Prerequisities: run 'compmon', create link file, then open ROOT before running this script.

void plotCompareAccs(){

  // Set x, y range here:--------
  // Run 1288:
  float x0Max = 20000000; 
  float x2Max = 10000000; 
  float y0Max = 5000; 
  float y2Max = 50000; 
  // -----------------------------

  /*
  // Get run number:
  TTree *t1 = (TTree*)rootfile->Get("snapshots");
  int runNumber;
  t1->SetBranchAddress("runNumber", &runNumber);
  TH1F *runHisto = new TH1F("runHisto", "runHisto", 100, 0, 1E6);
  Int_t entries = (Int_t)t1->GetEntries();
  printf("Entries: %d\n", entries);
  for (Int_t i=0; i<entries; i++){
    t1->GetEntry(i);
    runHisto->Fill(runNumber);
  }

  int runNum = runHisto->GetMean();
  delete runHisto;*/


  //TString c1Title;
  // c1Title.Form("Run %d: Comparison of Acc0,4 beam+laser states", runNum);
  TCanvas* c1 = new TCanvas("c1", "blh");
  TFile* rootfile = gROOT->GetFile();

  rootfile->cd("mpsHistos");

  TH1F *h0_beamOff_laserOn = hM_acc0_beamOff_LaserOn;//->Clone();
  TH1F *h0_beamOff_laserOff = hM_acc0_beamOff_LaserOff;
  TH1F *h0_beamOn_laserOn = hM_acc0_LaserOn;
  TH1F *h0_beamOn_laserOff = hM_acc0_LaserOff;
  //----changed the following from 2's to 4's for tests-----
  TH1F *h2_beamOff_laserOn = hM_acc4_beamOff_LaserOn;
  TH1F *h2_beamOff_laserOff = hM_acc4_beamOff_LaserOff;
  TH1F *h2_beamOn_laserOn = hM_acc4_LaserOn;
  TH1F *h2_beamOn_laserOff = hM_acc4_LaserOff;

  c1->Divide(1,2);


  // Acc0 histograms on top panel:

  c1->cd(1);

  printf("Drawing Acc0 plots...\n");

  h0_beamOff_laserOn->SetLineColor(kBlue);
  h0_beamOff_laserOff->SetLineColor(kRed);
  h0_beamOn_laserOn->SetLineColor(kGreen);
  h0_beamOn_laserOff->SetLineColor(kMagenta);

  h0_beamOff_laserOn->GetXaxis()->SetRangeUser(-1E7., x0Max);
  h0_beamOff_laserOn->GetYaxis()->SetRangeUser(0., y0Max);

  h0_beamOff_laserOn->GetXaxis()->SetTitle("Peak [sRAU]");
  h0_beamOff_laserOn->GetYaxis()->SetTitle("Counts");
  h0_beamOff_laserOn->SetTitle("Acc0");

  h0_beamOff_laserOn->Draw();
  h0_beamOff_laserOff->Draw("sames");
  h0_beamOn_laserOn->Draw("sames");
  h0_beamOn_laserOff->Draw("sames");
  gPad->Update(); // Need this to access stats boxes. Takes a lot longer
                  // but seems less prone to random segmentation breaks.

  TPaveStats *st1 = (TPaveStats*)h0_beamOff_laserOn->FindObject("stats");
  st1->SetX1NDC(0.8);
  st1->SetX2NDC(0.95);
  st1->SetY1NDC(0.8);
  st1->SetY2NDC(0.95);
  TPaveStats *st2 = (TPaveStats*)h0_beamOff_laserOff->FindObject("stats");
  st2->SetX1NDC(0.8);
  st2->SetX2NDC(0.95);
  st2->SetY1NDC(0.6);
  st2->SetY2NDC(0.77);
  TPaveStats *st3 = (TPaveStats*)h0_beamOn_laserOn->FindObject("stats");
  st3->SetX1NDC(0.8);
  st3->SetX2NDC(0.95);
  st3->SetY1NDC(0.4);
  st3->SetY2NDC(0.57);
  TPaveStats *st4 = (TPaveStats*)h0_beamOn_laserOff->FindObject("stats");
  st4->SetX1NDC(0.8);
  st4->SetX2NDC(0.95);
  st4->SetY1NDC(0.2);
  st4->SetY2NDC(0.37); 

  legend0 = new TLegend(0.15, 0.7, 0.3, 0.87);
  legend0->SetHeader("Legend");
  legend0->AddEntry(h0_beamOff_laserOn, "Beam off, Laser on", "l");
  legend0->AddEntry(h0_beamOff_laserOff, "Beam off, Laser off", "l");
  legend0->AddEntry(h0_beamOn_laserOn, "Beam on, Laser on", "l");
  legend0->AddEntry(h0_beamOn_laserOff, "Beam on, Laser off", "l"); 
  legend0->Draw();


  // Acc2 histograms on bottom panel:
  c1->cd(2);

  printf("Drawing Acc4 plots...\n");

  h2_beamOff_laserOn->SetLineColor(kBlue);
  h2_beamOff_laserOff->SetLineColor(kRed);
  h2_beamOn_laserOn->SetLineColor(kGreen);
  h2_beamOn_laserOff->SetLineColor(kMagenta);

  h2_beamOff_laserOn->GetXaxis()->SetRangeUser(-1E6., x2Max);
  h2_beamOff_laserOn->GetYaxis()->SetRangeUser(0., y2Max);

  h2_beamOff_laserOn->GetXaxis()->SetTitle("Peak [sRAU]");
  h2_beamOff_laserOn->GetYaxis()->SetTitle("Counts");
  h2_beamOff_laserOn->SetTitle("Acc4");

  h2_beamOff_laserOn->Draw();
  h2_beamOff_laserOff->Draw("sames");
  h2_beamOn_laserOn->Draw("sames");
  h2_beamOn_laserOff->Draw("sames");
  gPad->Update();

  TPaveStats *st1 = (TPaveStats*)h2_beamOff_laserOn->FindObject("stats");
  st1->SetX1NDC(0.8);
  st1->SetX2NDC(0.95);
  st1->SetY1NDC(0.8);
  st1->SetY2NDC(0.97);
  TPaveStats *st2 = (TPaveStats*)h2_beamOff_laserOff->FindObject("stats");
  st2->SetX1NDC(0.8);
  st2->SetX2NDC(0.95);
  st2->SetY1NDC(0.6);
  st2->SetY2NDC(0.77);
  TPaveStats *st3 = (TPaveStats*)h2_beamOn_laserOn->FindObject("stats");
  st3->SetX1NDC(0.8);
  st3->SetX2NDC(0.95);
  st3->SetY1NDC(0.4);
  st3->SetY2NDC(0.57);
  TPaveStats *st4 = (TPaveStats*)h2_beamOn_laserOff->FindObject("stats");
  st4->SetX1NDC(0.8);
  st4->SetX2NDC(0.95);
  st4->SetY1NDC(0.2);
  st4->SetY2NDC(0.37);

  legend2 = new TLegend(0.15, 0.7, 0.3, 0.87);
  legend2->SetHeader("Legend");
  legend2->AddEntry(h2_beamOff_laserOn, "Beam off, Laser on", "l");
  legend2->AddEntry(h2_beamOff_laserOff, "Beam off, Laser off", "l");
  legend2->AddEntry(h2_beamOn_laserOn, "Beam on, Laser on", "l");
  legend2->AddEntry(h2_beamOn_laserOff, "Beam on, Laser off", "l");
  legend2->Draw();

  
  printf("Done.\n");

  return;

}
