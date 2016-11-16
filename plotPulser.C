// Edited 2016-02-20 Larisa Thorne
// Want gr to be global
// ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch02s06.html

TGraphErrors *gr;

// -------- Main program ----------------------

void plotPulser(){
  
  // TFile *rootfile=new TFile("lnkCompMon.output");
  TFile *rootfile = gROOT->GetFile();
  TTree *t1 = (TTree*)rootfile->Get("pulserwise");

  // Tree variables
  Float_t p1; // In pulserwise
  Float_t p2; // In pulserwise
  Float_t p3; // In pulserwise
  Float_t p4; // In pulserwise
  int varIndex; // In pulserwise. Don't use varDAC
  Float_t delta_contribution;
  Float_t all_delta;
  const int nsteps = 200;
  Int_t Wanted_Dac_Index=13;
  char title[280];
  Float_t meanDiff[nsteps];
  Float_t meanVar[nsteps];
  Float_t entries[nsteps];
  Float_t rms[nsteps];
  Float_t sqrtentries[nsteps];
  Float_t error[nsteps];
  Float_t errorX[nsteps];
  sprintf(title,"Pulserwise Data for DAC Index %d",Wanted_Dac_Index);

  /* // Get run number:
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

  
  TCanvas* c1 = new TCanvas("c1",title);
  c1->Divide(2,3);
  TCanvas* c2 = new TCanvas("c2","All DAC Settings");
  c2->Divide(5,5);

  t1->SetBranchAddress("p1",&p1);
  t1->SetBranchAddress("p2",&p2);
  t1->SetBranchAddress("p3",&p3);
  t1->SetBranchAddress("p4",&p4);
  t1->SetBranchAddress("varIndex",&varIndex);
  
  
  char his_title[128];
  char his_label[128];

  TH1F* h_delta[nsteps];   // Delta minus background
  TH1F* h_var[nsteps];     // Var minus bckground
  TH1F* h_diff[nsteps];    // Both - var
  for(Int_t i=0;i<nsteps;i++){
    sprintf(his_title,"DAC Setting_%03d",i);
    //printf(his_title,"DAC Setting_%03d",i);
    sprintf(his_label,"h_delta_%03d",i);
    h_delta[i]=new TH1F(his_label,his_title,1000,-10000,10000);
    sprintf(his_label,"h_var_%03d",i);
    h_var[i]=new TH1F(his_label,his_title,6000,-10000,50000);
    sprintf(his_label,"h_diff_%03d",i);
    h_diff[i]=new TH1F(his_label,his_title,60000,1000,50000);
  }

  float xMax=30000;
  float xMin=-4000;
  TH1F*h1 = new TH1F("Var + Delta","Var + Delta",10000,xMin,xMax);
  TH1F*h2 = new TH1F("Delta","Delta",10000,xMin,xMax);
  TH1F*h3 = new TH1F("Var","Var",10000,xMin,xMax);
  TH1F*h4 = new TH1F("BG","BG",10000,xMin,xMax);  
  TH1F*h5 = new TH1F("Delta Contribution","Delta Contribution",1000,1000,10000);
  TH1F*h6 = new TH1F("All DAC Setting Delta Contrib.","All DAC Setting Delta Contrib.",1000,-10000,10000);

  Int_t nevents=(Int_t)t1->GetEntries(); // Found form in sample 
  printf("nevents: %d\n", nevents); // Must have nonzero value to proceed!

  Int_t buffer[nsteps];
  for(Int_t i=0;i<nsteps;i++){
    buffer[i]=i;
  }

  float pVar,pDelta,pBoth,pNothing;

  for(Int_t i=0;i<nevents;i++){
    t1->GetEntry(i);
    //if(i==0||i<100){
    //  cout<<"Dac Index is    "<<varIndex<<endl;
    // }
    pVar=p1;
    pDelta=p2;
    pBoth=p3;
    pNothing=p4;

    if(varIndex==Wanted_Dac_Index){
      h1->Fill(pBoth);
      h2->Fill(pDelta);
      h3->Fill(pVar);
      h4->Fill(pNothing);
      delta_contribution=pBoth-pVar;
      h5->Fill(delta_contribution);
    }
    all_delta=pBoth-pVar;
    h6->Fill(all_delta);
    if(varIndex<nsteps){
      h_delta[varIndex]->Fill(pDelta-pNothing);
      h_var[varIndex]->Fill(pVar-pNothing);
      h_diff[varIndex]->Fill(pBoth-pVar);
    }  
  }
  c1->cd(1);
  h1->Draw();
  c1->cd(2);
  h2->Draw();
  c1->cd(3);
  h3->Draw();
  c1->cd(4);
  h4->Draw();
  c1->cd(5);
  h5->Draw();
  c1->cd(6);
  h6->Draw();

  int nonzeroSteps=0;
  for(Int_t i=0;i<nsteps;i++){
    entries[i]=h_diff[i]->GetEntries();
    if(entries[i]==0) break;
    nonzeroSteps++;
    c2->cd(i+1);
    h_diff[i]->Draw();
    meanDiff[i]=h_diff[i]->GetMean();
    meanVar[i]=h_var[i]->GetMean();  
    rms[i]=h_diff[i]->GetRMS();
    sqrtentries[i]=sqrt(entries[i]);
    error[i]=rms[i]/sqrtentries[i];
    errorX[i]=0.;   //just to keep TGraph happy
    printf(" %4d Var %12.1f  Delta %10.1f  Diff %10.3f p/m %10.2f \n",
 	   i, h_var[i]->GetMean(),h_delta[i]->GetMean(), meanDiff[i],error[i]);
  }

  TCanvas* c3 = new TCanvas("c3","Linearity");

  //TGraph *gr = new TGraphErrors(nonzeroSteps, meanVar, meanDiff, errorX, error);
  gr = new TGraphErrors(nonzeroSteps, meanVar, meanDiff, errorX, error);

  gr->SetTitle("MM Pulser Linearity Plot");
  gr->GetXaxis()->SetTitle("VAR LED pulses (srau)");
  gr->GetYaxis()->SetTitle("(VAR+Delta)-VAR (srau)");
  gr->Draw("AP");

  /* You can store any named object in gROOT->GetListOfSpecials() with 
gROOT->GetListOfSpecials()->Add(myobject); 
and retrieve it later with gROOT->FindObject. */
  //gROOT->GetListOfSpecials()->Add(gr);
  //printf(gROOT->FindObject("gr"));

  printf("Done.\n\n");

}
