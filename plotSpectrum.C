void plotSpectrum(){
  //
  //  float xMax=40000.; //set x-max for compton histo
  float xMax=8000.; //set x-max for compton histo

  TCanvas* c1=new TCanvas("c1","Triggered PMT Pulse  Data");
  TFile* rootfile= gROOT->GetFile();

  rootfile->cd("mpsHistos");

  int OnMPS=0;
  int OffMPS=0;
  int OnAcceptedTriggers=0;
  int OffAcceptedTriggers=0;
  int OnTriggers=0;
  int OffTriggers=0;

  //sum up laser-on counts
  //  TH1F* h_Trigger=hM_Trigs_Scaler;  //use if IP scaler trigger info availbled
  TH1F* h_Trigger=hM_Trigs_Prescaled; //use if only latching scaler info available

  for(int i=1;i<=4;i++){
    OnMPS+ =hM_SpinState_BeamOn->GetBinContent(i);
    OnAcceptedTriggers+ =hM_Trig_Accepted->GetBinContent(i);
    OnTriggers+ = h_Trigger->GetBinContent(i);
  }
  //sum up laser-off counts
  for(int i=5;i<=8;i++){
    OffMPS+ =hM_SpinState_BeamOn->GetBinContent(i);
    OffAcceptedTriggers+ =hM_Trig_Accepted->GetBinContent(i);
    OffTriggers+ = h_Trigger->GetBinContent(i);
  }

  printf("Laser-On  MPS count       %10d\n",OnMPS);
  printf("          Accepted Triggs %10d\n",OnAcceptedTriggers);
  printf("          Triggers        %10d\n",OnTriggers);

  printf("Laser-Off MPS count       %10d\n",OffMPS);
  printf("          AcceptedTriggs  %10d\n",OffAcceptedTriggers);
  printf("          Triggers        %10d\n",OffTriggers);

  c1->Divide(1,2);
  c1->cd(1);
  //  rootfile->cd("mpsHistos");
  hM_SpinState_BeamOn->Draw();



  rootfile->cd("triggeredHistos");

  TH1F *h_laserOn =hTrig_sums_laserOn->Clone();
  h_laserOn->Sumw2();
  h_laserOn->SetTitle("Prescaled Counts per MPS (Laser On and Laser Off)");
  TH1F *h_laserOff =hTrig_sums_laserOff->Clone();
  h_laserOff->Sumw2();

  c1->cd(2);
  if(OnMPS>0&& OnAcceptedTriggers>0){
    float scale=OnTriggers;
    scale=float(OnTriggers)/(float(OnMPS)*float(OnAcceptedTriggers));
    printf("Laser-on scale factor %f\n",scale);
    h_laserOn->Scale(scale);
    h_laserOn->GetXaxis()->SetRangeUser(-1000.,xMax);
    h_laserOn->Draw();
  }
  if(OffMPS>0&& OffAcceptedTriggers>0){
    float scale=OnTriggers;
    scale=float(OffTriggers)/(float(OffMPS)*float(OffAcceptedTriggers));
    h_laserOff->Scale(scale);
    printf("Laser-off scale factor %f\n",scale);
    h_laserOff->GetXaxis()->SetRangeUser(-1000.,xMax);
    h_laserOff->SetLineColor(kRed);
    h_laserOff->Draw("same");
  }
  TH1F* h_signal=h_laserOn->Clone();
  h_signal->Add(h_laserOff,-1.0);
  h_signal->SetLineColor(kGreen);
  h_signal->Draw("same");


  //----- Legend -----

  //legend = new TLegend(0.1, 0.7, 0.3, 0.9); // Top left corner
  legend = new TLegend(0.8, 0.5, 0.95, 0.7); // Below stats box
  legend->SetHeader("Legend");
  legend->AddEntry(h_signal, "h_signal", "l");
  legend->AddEntry(h_laserOff, "h_laserOff", "l");
  legend->AddEntry(h_laserOn, "h_laserOn", "l");
  legend->Draw("same");

  return;
}

