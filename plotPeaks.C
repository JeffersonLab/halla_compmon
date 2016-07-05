void plotPeaks(){

  // Read snapshot info from tree and plot
  // Make histogram of pulse peaks (max deviation from pedestal)
  // Also histos channels 1 and 2 for pedestal check
  // and peak height vs integral
  // ***Use with real data, not pulser
  // TFile *rootfile=new TFile("lnkCompMon.output");

  // 2016-02-10 @ 10:30am: beamOn -> beamState. Fixed branch address issue. Works.

  TFile *rootfile = gROOT->GetFile();
  //  TTree *t1=(TTree*)rootfile->Get("pulserwise");
  TTree *t1=(TTree*)rootfile->Get("snapshots");

  // Tree variables
  int numSamples;
  int mpsCoda;
  Float_t samples[1000];
  int laserState;
  int beamState; 
  int peakRAU;  // Peak in raw ADC units
  // Float_t pedestal=3849.;
  Float_t pedestal = 3825; // For 5-pass
  //Float_t pedestal=3850.; // Was 2370 to work for Spring 2015 run, making it larger (3850) for Spring 2016 - use pedestal found from .x plotSnaps.C 
  // Float_t pedestal=2370.; // Was 2370 to work for Spring 2015 run, making it larger (3850) for Spring 2016 - use pedestal found from .x plotSnaps.C 
  Float_t maxAccum=28000.; // Was 10000 for Spring 2015 run, changed to 28000 for Spring 2016

  // Get run number:
  int runNumber;
  t1->SetBranchAddress("runNumber", &runNumber);
  TH1F *runHisto = new TH1F("runHisto", "runHisto", 100, 0, 1E6);
  Int_t entries = (Int_t)t1->GetEntries();
  for (Int_t i=0; i<entries; i++){
    t1->GetEntry(i);
    runHisto->Fill(runNumber);
  }
  int runNum = runHisto->GetMean();
  delete runHisto;


  // Read tree variables:
  t1->SetBranchAddress("numSamples",&numSamples);
  t1->SetBranchAddress("snap",&samples);
  t1->SetBranchAddress("mpsCoda",&mpsCoda);
  t1->SetBranchAddress("laserState",&laserState);
  t1->SetBranchAddress("beamState",&beamState); 

  TH1F* hPeak= new TH1F("hPeak","Counts vs Wave Form Peak (small RAU=big peak)", 1024, 0., 4096.);
  TH1F* hPed0= new TH1F("hPed0","WF chan 0", 4096, 0., 4096.);
  TH1F* hPed1= new TH1F("hPed1","WF chan 1", 4096, 0., 4096.);
  TH1F* hPed2= new TH1F("hPed2","WF chan 2", 4096, 0., 4096.);
  TH1F* hsum=  new TH1F("hsum","ped-summedWF", 1000, -5000., maxAccum);
  TH2F* hpeakVsum= new TH2F("hpeakVsum","Peak vs Sum", 500, -5000., maxAccum, 512, 0, 4096);

  Int_t ntriggers=(Int_t)t1->GetEntries();
  printf("Number of waveform samples in file : %d\n", ntriggers);
  Float_t sum;

  for (int i=2; i<ntriggers; i++){ // Ignore first two events
    t1->GetEntry(i);
    if (beamState==0) continue; // Beam off: skip this entry
    if (laserState!=1) continue; // Laser off: skip this entry
    peakRAU = 1.E8; // Where this from?
    sum = 0;
    for(int j=0; j<numSamples; j++){
      sum+=pedestal-samples[j];  // Invert sums for sanity
      if (samples[j]<peakRAU) peakRAU=samples[j];
    }
    hPeak->Fill(peakRAU);
    hPed0->Fill(samples[0]);
    hPed1->Fill(samples[1]);
    hPed2->Fill(samples[2]);
    if(peakRAU>0){
      hsum->Fill(sum);
      hpeakVsum->Fill(sum,peakRAU);
    }
  }


  TString c1Title;
  c1Title.Form("Run %d: Waveform peaks", runNum);
  TCanvas* c1=new TCanvas("c1", c1Title);
  hPeak->GetXaxis()->SetTitle("Peak");
  hPeak->GetYaxis()->SetTitle("Counts");
  hPeak->Draw();
  int peakChannel = hPed1->GetMaximumBin();

  Double_t xmax = hPed1->GetBinCenter(peakChannel);
  printf("Channel 1 max counts at x=%d\n", xmax);  
  hPed0->SetAxisRange(xmax-10.,xmax+10.);
  hPed1->SetAxisRange(xmax-10.,xmax+10.);
  hPed1->SetLineColor(2);
  hPed2->SetAxisRange(xmax-10.,xmax+10.);
  hPed2->SetLineColor(8);
  

  TCanvas* c2 = new TCanvas("c2","FADC Pedestals", 20, 20, 800, 300);
  hPed0->Draw();
  hPed1->Draw("sames");
  hPed2->Draw("sames");
  hPed0->SetTitle("Waveforms in Channels 0-2");
  gPad->Update();

  TPaveStats *st1 = (TPaveStats*)hPed0->FindObject("stats");
  st1->SetX1NDC(0.8);
  st1->SetX2NDC(0.95);
  st1->SetY1NDC(0.8);
  st1->SetY2NDC(0.95);
  TPaveStats *st2 = (TPaveStats*)hPed1->FindObject("stats");
  st2->SetX1NDC(0.8);
  st2->SetX2NDC(0.95);
  st2->SetY1NDC(0.6);
  st2->SetY2NDC(0.77);
  TPaveStats *st3 = (TPaveStats*)hPed2->FindObject("stats");
  st3->SetX1NDC(0.8);
  st3->SetX2NDC(0.95);
  st3->SetY1NDC(0.4);
  st3->SetY2NDC(0.57);

  legend = new TLegend(0.6,0.73,0.79,0.85);
  legend->SetHeader("Legend");
  legend->AddEntry(hPed0, "hPed0", "l");
  legend->AddEntry(hPed1, "hPed1", "l");
  legend->AddEntry(hPed2, "hPed2", "l");
  legend->Draw("same");


  TCanvas* c3 = new TCanvas("c3","Peak and Pulse Sums", 30, 30, 800, 600);
  c3->Divide(1,2);

  c3->cd(1);
  hsum->Draw();
  hsum->SetTitle("Pedestal summed waveform");
  hsum->GetXaxis()->SetTitle("Sum");
  hsum->GetYaxis()->SetTitle("Peak [sRAU]");

  c3->cd(2);
  hpeakVsum->GetXaxis()->SetTitle("Sum");
  hpeakVsum->GetYaxis()->SetTitle("Peak [sRAU]");
  hpeakVsum->Draw();

  printf("Done.\n");

}

