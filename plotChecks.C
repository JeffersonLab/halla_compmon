void plotChecks(){

  TFile* rootfile= gROOT->GetFile();
  rootfile->cd("mpsHistos");
  TCanvas* c1=new TCanvas("c1","Beam And Laser Cuts");
  c1->Divide(1,3);
  //Beam Current from BCM
  float xLow=-10.;
  float xHi=30.;
  float ymax=1.1*hM_BCM->GetMaximum();
  c1->cd(1);
  hM_BCM->SetMaximum(ymax);
  hM_BCM->GetXaxis()->SetRangeUser(xLow,xHi);
  hM_BCM->Draw();
  hM_BCM_BeamOn->SetLineColor(kBlue);
  hM_BCM_BeamOn->Draw("same");
  hM_BCM_BeamOff->SetLineColor(kRed);
  hM_BCM_BeamOff->Draw("same");
  //Beam Current from BPMSum
  xLow=-10.;
  xHi=30.;
  ymax=1.1*hM_BPMSum->GetMaximum();
  c1->cd(2);
  hM_BPMSum->SetMaximum(ymax);
  hM_BPMSum->GetXaxis()->SetRangeUser(xLow,xHi);
  hM_BPMSum->Draw();
  hM_BPMSum_BeamOn->SetLineColor(kBlue);
  hM_BPMSum_BeamOn->Draw("same");
  hM_BPMSum_BeamOff->SetLineColor(kRed);
  hM_BPMSum_BeamOff->Draw("same");
  //Laser Power
  xLow=-500.;
  xHi=3000.;
  ymax=1.1*hM_Cav_Power->GetMaximum();
  c1->cd(3);
  hM_Cav_Power->SetMaximum(ymax);
  hM_Cav_Power->GetXaxis()->SetRangeUser(xLow,xHi);
  hM_Cav_Power->Draw();
  hM_Cav_Power_LaserOn->SetLineColor(kBlue);
  hM_Cav_Power_LaserOn->Draw("same");
  hM_Cav_Power_LaserOff->SetLineColor(kRed);
  hM_Cav_Power_LaserOff->Draw("same");


  // CANVAS 2
  TCanvas* c2=new TCanvas("c2","MPS Counts",3);
  c2->Divide(1,2);
  hM_SpinState->Draw();
  hM_SpinState_BeamOn->SetLineColor(kBlue);
  hM_SpinState_BeamOn->Draw("same");
}
