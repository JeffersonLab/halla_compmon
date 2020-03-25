void exampleMacro1_old(){
  //cout<<"BLAH"<<endl;
  TPad *myPad = (TPad*)gPad;
  myPad->Divide(2,1);
  myPad->cd(1);
  TH1F *h1 = new TH1F("h1", "Histo", 200, -10, 10);
  h1->FillRandom("gaus", 10000);
  h1->SetLineColor(1);
  h1->Draw();
  TH1F *h2 = new TH1F("h2", "Histo", 200, -10, 10);
  h2->FillRandom("gaus", 5000);
  h2->SetLineColor(2);
  h2->Draw("same");
  myPad->cd(2);
  TH1F *h3 = new TH1F("h3", "Not Histo", 200, -10, 10);
  h3->FillRandom("gaus", 10000);
  h3->SetLineColor(3);
  h3->Draw();
  TH1F *h4 = new TH1F("h4", "Not Histo", 200, -10, 10);
  h4->FillRandom("gaus", 2500);
  h4->SetLineColor(4);
  h4->Draw("same");
}

void exampleMacro1(){
  TTree *snapshots = (TTree *)gDirectory->Get("snapshots");
  TCanvas *c1 = new TCanvas("c1", "Eat my ass", 1400, 800);
  TPad *myPad = new TPad("snapshots_pad", "Snapshots Pad", 0, 0, 1, 1); 
  myPad->cd(); myPad->Draw();
  float snapshot[220];
  float bcm;
  int randomTime, mpsCoda, numSamples, runNumber;
  int run_num = 0;

  snapshots->SetBranchAddress("randomTime", &randomTime);
  snapshots->SetBranchAddress("snap", &snapshot);
  snapshots->SetBranchAddress("mpsCoda", &mpsCoda);
  snapshots->SetBranchAddress("numSamples", &numSamples);
  snapshots->SetBranchAddress("runNumber", &runNumber);
  snapshots->SetBranchAddress("bcm", &bcm);

  TGraph *g_heightSum = new TGraph();
  int accepted = 0;

  for(int i = 0; i < snapshots->GetEntries(); i++){
    snapshots->GetEntry(i); //if(i == 0){run_num = runNumber;}
    g_heightSum->SetPoint(i, i, 50);
  }

  TH1F *h1 = new TH1F("h1", "Histo", 200, -10, 10);
  h1->FillRandom("gaus", 10000);
  h1->SetLineColor(1);
  h1->Draw();
}
