// Last edit: Larisa Thorne @ 2016-03-30
// Description: used for run analysis. Plots Compton spectrum, asymmetry and difference of sums of the two beam on states. BCM and cavity power separately.

void plotAsymSpectrum(Int_t runNum){

  //  float xMax=40000.; //set x-max for compton histo
  float xMax=32000.; // Set x max for Compton histogram
  Int_t epicsCutoff = 30; // Same as in buildQuartets: defines when Epics variables haven't changed in a long time
  Int_t onColor = kRed + 2; // For on states
  Int_t offColor = kCyan + 3; // For off states
  Int_t otherColor = kBlue + 3; // For single histogram plots
  Int_t laserOn0Color = kMagenta + 1; // For laser on anti-aligned state
  Int_t laserOn1Color = kGreen + 2; // For laser on aligned state
  Int_t laserOn2Color = kOrange + 1; // For laser off state

  char inputFile[60];
  sprintf(inputFile, "/home/compton/cornejo/rootfiles/compmon_%d.root", runNum);
  TFile *rootfile = new TFile(inputFile);
  rootfile->cd("mpsHistos");

  // Collapse spin states 0-7 into laser on/laser off, helicity/polarization 
  // LN and RP state 0
  // LP and RN state 1
  int OnState[] = {0,1,1,0,2,2,2,2};

  // Counters:
  int MPS[3] = {0, 0, 0};
  int Triggers[3] = {0, 0, 0};
  int AcceptedTriggers[3] = {0, 0, 0};

  // Sum up laser-on counts:
  TH1F* h_Trigger = hM_Trigs_Prescaled; //use if just latching scaler info available
  TH1F* h_BCMon = hM_BCM; h_BCMon->SetName("h_BCMon");
  TH1F* h_BCMoff = hM_BCM_BeamOn; h_BCMoff->SetName("h_BCMoff");
  TH1F* h_CavPowerOn = hM_Cav_Power_LaserOn; h_CavPowerOn->SetName("h_CavPowerOn");
  TH1F* h_CavPowerOff = hM_Cav_Power_LaserOff; h_CavPowerOff->SetName("h_CavPowerOff");

  int state;
  for (int i=0; i<=7; i++){
    state = OnState[i];
    if (state>=0){
      MPS[state] += hM_SpinState_BeamOn->GetBinContent(i+1);
      AcceptedTriggers[state] += hM_Trig_Accepted->GetBinContent(i+1);
      Triggers[state] += h_Trigger->GetBinContent(i+1);
    }
  }

  rootfile->cd("triggeredHistos");

  // Three histograms are spin-polarization antiparallal, parallel, and
  // laser off.  (First two might not be in that order)
  TH1F* h_laserOn[3];

  h_laserOn[0] = (TH1F*) hTrig_sums_L_N->Clone();
  h_laserOn[0]->Add(hTrig_sums_R_P);
  h_laserOn[0]->Sumw2();

  h_laserOn[1] = (TH1F*) hTrig_sums_L_P->Clone();
  h_laserOn[1]->Add(hTrig_sums_R_N);
  h_laserOn[1]->Sumw2();

  h_laserOn[2] = (TH1F*) hTrig_sums_laserOff->Clone();
  h_laserOn[2]->Sumw2();

  h_laserOn[0]->SetName("h_laserOn[0]");
  h_laserOn[1]->SetName("h_laserOn[1]");
  h_laserOn[2]->SetName("h_laserOn[2]");
  h_laserOn[0]->SetTitle("Prescaled Counts per MPS *=(2 spin states, laser off)");

  for (int i=0; i<3; i++){
    printf("State: %2d  MPS count: %5d\tTriggers: %7d\tAccepted: %6d\n",
	   i, MPS[i], Triggers[i], AcceptedTriggers[i] );
  }

  // -----------------------------------------------------------
  // ----------- Get runwise variables -------------------------
  // -----------------------------------------------------------


  // Cavity power on min:
  TTree *t1 = (TTree*)rootfile->Get("runwise");
  Float_t cavityPowerOnMin;
  t1->SetBranchAddress("cavityPowerOnMin", &cavityPowerOnMin);
  TH1F *cavPowerOn = new TH1F("cavPowerOn", "cavPowerOn", 100, 0, 1E7);
  Int_t entries = (Int_t)t1->GetEntries();
  for (int i=0; i<entries; i++){
    t1->GetEntry(i);
    cavPowerOn->Fill(cavityPowerOnMin);
  }
  Float_t cavPowerOnMin = cavPowerOn->GetMean();
  delete cavPowerOn;
  printf("\n\tcavPowerOnMin: %d\n", cavPowerOnMin);


  // Cavity power off max:
  Float_t cavityPowerOffMax;
  t1->SetBranchAddress("cavityPowerOffMax", &cavityPowerOffMax);
  TH1F *cavPowerOff = new TH1F("cavPowerOff", "cavPowerOff", 100, 0, 1E7);
  Int_t entries = (Int_t)t1->GetEntries();
  for (int i=0; i<entries; i++){
    t1->GetEntry(i);
    cavPowerOff->Fill(cavityPowerOffMax);
  }
  Float_t cavPowerOffMax = cavPowerOff->GetMean();
  delete cavPowerOff;
  printf("\tcavPowerOffMax: %d\n\n", cavPowerOffMax);


  // Plot of countEpics distribution:
  TTree *t2 = (TTree*)rootfile->Get("quartetwise");
  Int_t countEpics; // Counts how long since new Epics event
  t2->SetBranchAddress("countEpics", &countEpics);

  Int_t prevEpics = 0;
  Int_t sameEpicsCounter = 0;
  Int_t individualEpicsEntries = 0;

  TH1F *epicsCount = new TH1F("epicsCount", "epicsCount", 150, 0, 150);
  Int_t ntriggers = (Int_t)t2->GetEntries();

  for (int i=0; i<ntriggers; i++){
    t2->GetEntry(i);

    // Check if Epics variables updated within last 30 entries:
    if (countEpics == prevEpics){
      sameEpicsCounter++;
    }
    else {
      prevEpics = countEpics;
      epicsCount->Fill(sameEpicsCounter);
      sameEpicsCounter = 0;
      individualEpicsEntries++;
    }
  }


  // ---------------------------------------------------------------
  // ------------- Plot laser/beam state histograms ----------------
  // ---------------------------------------------------------------

  TString c1Title;
  c1Title.Form("Run %d: Triggered PMT Pulse Data", runNum);
  TCanvas* c1 = new TCanvas("c1", c1Title, 10, 10, 1150, 1100); 
  
  c1->Divide(2,3);
  c1->cd(1);

  h_laserOn[0]->SetLineColor(laserOn0Color); // Green
  h_laserOn[1]->SetLineColor(laserOn1Color); // Laser off?
  h_laserOn[2]->SetLineColor(laserOn2Color); // Laser on?
  h_laserOn[0]->GetXaxis()->SetTitle("Peaks [sRAU]");
  h_laserOn[0]->GetYaxis()->SetTitle("Counts");

  // Normalize per MPS and plot all three histograms
  for(int state=0; state<3; state++){
    if(MPS[state]>0 && AcceptedTriggers[state]>0){
      float scale = Triggers[state];
      scale = float(Triggers[state])/(float(MPS[state])*float(AcceptedTriggers[state]));
      printf("State: %2d\tScale factor: %8.6f\n", state, scale);
      h_laserOn[state]->Scale(scale);
      h_laserOn[state]->GetXaxis()->SetRangeUser(-1000.,xMax);
      
      if(state==0){
	h_laserOn[state]->Draw();
      }
      else{
	h_laserOn[state]->Draw("sames");
      }
    }
  }
  gPad->Update();

  TPaveStats *st1 = (TPaveStats*)h_laserOn[0]->FindObject("stats");
  st1->SetX1NDC(0.8);
  st1->SetX2NDC(0.95);
  st1->SetY1NDC(0.8);
  st1->SetY2NDC(0.95);
  TPaveStats *st2 = (TPaveStats*)h_laserOn[1]->FindObject("stats");
  st2->SetX1NDC(0.8);
  st2->SetX2NDC(0.95);
  st2->SetY1NDC(0.6);
  st2->SetY2NDC(0.77);
  TPaveStats *st3 = (TPaveStats*)h_laserOn[2]->FindObject("stats");
  st3->SetX1NDC(0.8);
  st3->SetX2NDC(0.95);
  st3->SetY1NDC(0.4);
  st3->SetY2NDC(0.57);

  hlegend = new TLegend(0.6, 0.70, 0.79, 0.89);
  hlegend->AddEntry(h_laserOn[0], "Laser on (LN, RP)", "l");
  hlegend->AddEntry(h_laserOn[1], "Laser on (LP, RN)", "l");
  hlegend->AddEntry(h_laserOn[2], "Laser off", "l");
  hlegend->Draw();

  c1->cd(2);


  // ---------------- Plot asymmetry ----------------------------

  // Build asymmetry histogram first
  TH1F* h_diff=(TH1F*)h_laserOn[1]->Clone();
  h_diff->Add(h_laserOn[0],-1.0);  // Subtract off other spin combo
  TH1F* h_sum= (TH1F*)h_laserOn[0]->Clone();
  h_sum->Add(h_laserOn[1],+1.0);  // Add on other spin combo
  h_sum->Add(h_laserOn[2],-2.0); // Subtract off laser off
  TH1F* h_asym=(TH1F*) h_diff->Clone();
  h_asym->SetTitle("Asymmetry spectrum");
  h_asym->GetXaxis()->SetTitle("Peaks [sRAU]");
  h_asym->GetYaxis()->SetTitle("Ratio");
  h_asym->SetLineColor(otherColor);
  h_asym->Divide(h_sum);
  h_asym->SetMaximum(0.2);
  h_asym->SetMinimum(-0.2);
  h_asym->Draw();

  // Draw line at zero:
  TLine *l1 = new TLine(-1E3, 0, xMax, 0);
  l1->SetLineColor(15);
  l1->Draw();

  c1->cd(3);


  // ------------- Plot difference of sums ---------------------

  TH1F* h_diffOfSums = (TH1F*)h_laserOn[0]->Clone();
  h_diffOfSums->Add(h_laserOn[1], +1.0);
  h_diffOfSums->Add(h_laserOn[2], -2.0);
  
  h_diffOfSums->GetXaxis()->SetTitle("Peaks [sRAU]");
  h_diffOfSums->GetYaxis()->SetTitle("Counts");
  Double_t diffSumsMax = h_diffOfSums->GetMaximum();
  h_diffOfSums->SetMaximum(diffSumsMax);
  h_diffOfSums->SetTitle("Sums: background");
  h_diffOfSums->SetName("h_diffOfSums");
  h_diffOfSums->SetLineColor(otherColor);
  h_diffOfSums->Draw();




  // ---------------------------------------------------------------
  // ----------- Plot BCM, cavity power histograms -----------------
  // ---------------------------------------------------------------

  //TString c2Title;
  //c2Title.Form("Run %d: BCM & Cavity power, epicsCounts", runNum);
  //TCanvas *c2 = new TCanvas("c2", c2Title, 610, 10, 500, 600);
  //c2->Divide(1,3);
  
  c1->cd(4);
  h_BCMon->SetTitle("BCM state");
  h_BCMon->SetLineColor(onColor);
  h_BCMoff->SetLineColor(offColor);
  h_BCMon->GetXaxis()->SetRangeUser(-4, 16);
  h_BCMon->GetXaxis()->SetTitle("Charge");
  h_BCMon->GetYaxis()->SetTitle("Counts");
  h_BCMon->Draw();
  h_BCMoff->Draw("sames");
  gPad->Update();

  TPaveStats *st1 = (TPaveStats*)h_BCMon->FindObject("stats");
  st1->SetX1NDC(0.8);
  st1->SetX2NDC(0.95);
  st1->SetY1NDC(0.8);
  st1->SetY2NDC(0.95);
  TPaveStats *st2 = (TPaveStats*)h_BCMoff->FindObject("stats");
  st2->SetX1NDC(0.8);
  st2->SetX2NDC(0.95);
  st2->SetY1NDC(0.6);
  st2->SetY2NDC(0.77);

  legend2 = new TLegend(0.6,0.73,0.79,0.89);
  legend2->AddEntry(h_BCMon, "On", "l");
  legend2->AddEntry(h_BCMoff, "Off", "l");
  legend2->Draw();


  c1->cd(5);

  h_CavPowerOn->SetTitle("Cavity Power state");
  h_CavPowerOn->GetXaxis()->SetRangeUser(-1E3, 4E3);
  h_CavPowerOn->GetXaxis()->SetTitle("Power");
  h_CavPowerOn->GetYaxis()->SetTitle("Counts");
  h_CavPowerOn->SetLineColor(onColor);
  h_CavPowerOff->SetLineColor(offColor);
  h_CavPowerOn->Draw();
  h_CavPowerOff->Draw("sames");
  gPad->Update();

  TPaveStats *st1 = (TPaveStats*)h_CavPowerOn->FindObject("stats");
  st1->SetX1NDC(0.8);
  st1->SetX2NDC(0.95);
  st1->SetY1NDC(0.8);
  st1->SetY2NDC(0.95);
  TPaveStats *st2 = (TPaveStats*)h_CavPowerOff->FindObject("stats");
  st2->SetX1NDC(0.8);
  st2->SetX2NDC(0.95);
  st2->SetY1NDC(0.6);
  st2->SetY2NDC(0.77);

  legend3 = new TLegend(0.6,0.73,0.79,0.89);
  legend3->AddEntry(h_CavPowerOn, "On", "l");
  legend3->AddEntry(h_CavPowerOff, "Off", "l");
  legend3->Draw();

  // Draw thresholds for on/off:
  TLine *offMax = new TLine(cavityPowerOffMax, -150, cavityPowerOffMax, 1000);
  TLine *onMin = new TLine(cavityPowerOnMin, -150, cavityPowerOnMin, 1000);

  //offMax->SetLineColor(offColor);
  //onMin->SetLineColor(onColor);
  offMax->Draw();
  onMin->Draw("same");
  TLatex text1, text2; 
  text1.SetTextFont(43); text2.SetTextFont(43);
  text1.SetTextSize(14); text2.SetTextSize(14);
  //text1.SetTextColor(offColor); text2.SetTextColor(onColor);
  text1.DrawLatex(cavityPowerOffMax - 80, -800, "P_{off,max}");
  text2.DrawLatex(cavityPowerOnMin, -800, "P_{on,min}");
  

  // -------------------- epicsCount -------------------------

  c1->cd(6);
  gStyle->SetOptStat(111111);
  
  epicsCount->GetXaxis()->SetTitle("# of events per same epics");
  epicsCount->GetYaxis()->SetTitle("Counts");
  // Make y axis log
  //epicsCount->GetXaxis()->SetRangeUser(0, 100);
  epicsCount->SetMaximum(1E3);
  c1->SetLogy(1);;
  epicsCount->SetLineColor(otherColor);
  epicsCount->Draw();


  TLine *cutoffLine = new TLine(epicsCutoff, -20, epicsCutoff, 150);
  cutoffLine->SetLineColor(15);
  cutoffLine->Draw();
  TLatex textCutoff;
  textCutoff.SetTextFont(43); textCutoff.SetTextSize(12);
  textCutoff.DrawLatex(epicsCutoff - 10, -110, "epicsCutoff");

    TString saveFileName1; saveFileName1.Form("/home/compton/lthorne/CompMon/plots/plotAsymSpectrumRuns/run%dAsyms.png", runNum);
  c1->SaveAs(saveFileName1);


  printf("\nDone.\n");

  return;

}
