// Last edit: 20190301 Juan Carlos Cornejo
//   Make it work with ROOT 6 and push hard coded values to the top
//   Also fix correct sequence (Delta,Both,Nothing,Var,....)
// Edited 2016-02-20 Larisa Thorne
// Want gr to be global
// ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch02s06.html

// Set the plot limits here
const Float_t kDeltaLow  = -10000*0-2000;
const Float_t kDeltaHigh =  10000*0+5000*0+20000;
const Float_t kVarLow    = -10000*0-3000;
const Float_t kVarHigh   =  50000*0+55000*0 + 90e3;
const Float_t kDiffLow   =   1000*0-8000;
const Float_t kDiffHigh  =  50000*0+8000*0+24000;
const Float_t kBothLow   = -10000;
const Float_t kBothHigh  = 1.30*kVarHigh;
const Float_t kNothingLow = -4000;
const Float_t kNothingHigh = 30000;

const Int_t   kDeltaBins =  1000*0+500;
const Int_t   kVarBins   =  6000;
const Int_t   kDiffBins  = 60000*0+500;
const Int_t   kBothBins  = 10000*0+500;
const Int_t  kNothingBins = 10000*0+500;

TGraphErrors *gr;

// -------- Main program ----------------------

//void plotPulser(Int_t Wanted_Dac_Index=13){
void plotPulser(Int_t runnum,Int_t Wanted_Dac_Index=13){
  
  // TFile *rootfile=new TFile("lnkCompMon.output");
  //TFile *rootfile = gROOT->GetFile();
  //TTree *t1 = (TTree*)rootfile->Get("pulserwise");
  TChain *t1 = new TChain("pulserwise");
  TString filesPre = Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),runnum);
  Int_t nfiles = t1->Add(filesPre+".root");
  nfiles += t1->Add(filesPre+"_*.root");
  if(nfiles<=0) {
    std::cerr << "Looked for files under: " << filesPre+".root" << std::endl;
    std::cerr << "Found no files to plot!" << std::endl;
    return;
  }
  std::cout << "Found " << nfiles << (nfiles==1?" file":" files") << " for run " << runnum << std::endl;

  // Tree variables
  Float_t p1; // In pulserwise
  Float_t p2; // In pulserwise
  Float_t p3; // In pulserwise
  Float_t p4; // In pulserwise
  int varIndex; // In pulserwise. Don't use varDAC
  int varDAC; // We also want to know what varDAC was
  int helicityState; // Just a test
  int helicityStateReported; // Just a test
  Float_t delta_contribution;
  Float_t all_delta;
  const int nsteps = 200;
  //Int_t Wanted_Dac_Index=13;
  char title[280];
  Float_t meanDiff[nsteps];
  Float_t meanVar[nsteps];
  Float_t entries[nsteps];
  Float_t rms[nsteps];
  Float_t sqrtentries[nsteps];
  Float_t error[nsteps];
  Float_t errorX[nsteps];
  Int_t settingDAC[nsteps];
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
  t1->SetBranchAddress("varDAC",&varDAC);
  //t1->SetBranchAddress("helicityState",&helicityState);
  //t1->SetBranchAddress("helicityStateReported",&helicityStateReported);
  
  
  char his_title[128];
  char his_label[128];

  TH1F* h_delta[nsteps];   // Delta minus background
  TH1F* h_var[nsteps];     // Var minus bckground
  TH1F* h_diff[nsteps];    // Both - var
  for(Int_t i=0;i<nsteps;i++){
    sprintf(his_title,"DAC Setting_%03d",i);
    //printf(his_title,"DAC Setting_%03d",i);
    sprintf(his_label,"h_delta_%03d",i);
    h_delta[i]=new TH1F(his_label,his_title,kDeltaBins,kDeltaLow,kDeltaHigh);
    sprintf(his_label,"h_var_%03d",i);
    h_var[i]=new TH1F(his_label,his_title,kVarBins,kVarLow,kVarHigh);
    sprintf(his_label,"h_diff_%03d",i);
    h_diff[i]=new TH1F(his_label,his_title,kDiffBins,kDiffLow,kDiffHigh);
  }

  // (20190305 Cornejo): New, make these use the same limits as the new limits
  // defined at the top of this script
  //float xMax=30000;
  //float xMin=-4000;
  //TH1F*h1 = new TH1F("Var + Delta","Var + Delta",10000,xMin,xMax);
  //TH1F*h2 = new TH1F("Delta","Delta",10000,xMin,xMax);
  //TH1F*h3 = new TH1F("Var","Var",10000,xMin,xMax);
  //TH1F*h4 = new TH1F("BG","BG",10000,xMin,xMax);  
  //TH1F*h5 = new TH1F("Delta Contribution","Delta Contribution",1000,1000*0,10000*0);
  //TH1F*h6 = new TH1F("All DAC Setting Delta Contrib.","All DAC Setting Delta Contrib.",1000,-10000,10000);
  TH1F*h1 = new TH1F("Var + Delta","Var + Delta",kBothBins,kBothLow,kBothHigh);
  TH1F*h2 = new TH1F("Delta","Delta",kDeltaBins,kDeltaLow,kDeltaHigh);
  TH1F*h3 = new TH1F("Var","Var",kVarBins,kVarLow,kVarHigh);
  TH1F*h4 = new TH1F("BG","BG",kNothingBins,kNothingLow,kNothingHigh);
  TH1F*h5 = new TH1F("Delta Contribution","Delta Contribution",kDiffBins,kDiffLow,kDiffHigh);
  TH1F*h6 = new TH1F("All DAC Setting Delta Contrib.","All DAC Setting Delta Contrib.",kDiffBins,kDiffLow,kDiffHigh);


  Int_t nevents=(Int_t)t1->GetEntries(); // Found form in sample 
  printf("nevents: %d\n", nevents); // Must have nonzero value to proceed!

  Int_t buffer[nsteps];
  for(Int_t i=0;i<nsteps;i++){
    buffer[i]=i;
    settingDAC[i]=-1;
  }

  float pVar,pDelta,pBoth,pNothing;

  for(Int_t i=20;i<nevents;i++){
    t1->GetEntry(i);

    //if(i==0||i<100){
    //  cout<<"Dac Index is    "<<varIndex<<endl;
    // }
    // (cornejo 20190305: these seem to be different than what was in the notes,
    // that is, seems the sync now comes at the end with Variable (which
    // Brian somewhat recalls being the case now for the JLab setup)
    //pVar=p1;
    //pDelta=p2;
    //pBoth=p3;
    //pNothing=p4;
    pDelta=p1;
    pBoth=p2;
    pNothing=p3;
    pVar=p4;
    // 2019-08-09 (cornejo): How did the pattern change again? How is this possible?
    // But it did! Using the above bg all of a sudden becomes Variable, sigh..
    // I'll just fix it now...sigh...
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
      if(settingDAC[varIndex]==-1) {
        settingDAC[varIndex]=varDAC;
      }
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
    printf(" %4d Var %12.1f  Delta %10.1f  Diff %10.3f p/m %10.2f [DAC=%5d]\n",
 	   i, h_var[i]->GetMean(),h_delta[i]->GetMean(), meanDiff[i],error[i],settingDAC[i]);
  }

  TCanvas* c3 = new TCanvas("canvLinearity","Linearity");

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

  float var0,var1,dvar,dvardac;
  int dac0,dac1,ddac;
  for(Int_t i=1;i<nsteps*0;i++){
    entries[i]=h_diff[i]->GetEntries();
    if(entries[i]==0) break;
    var1=h_var[i-1]->GetMean();
    var0=h_var[i]->GetMean();
    dac1=settingDAC[i-1];
    dac0=settingDAC[i];
    dvar=var0-var1;
    ddac=dac0-dac1;
    dvardac = dvar/float(ddac);
    printf(" %4d Var: %9.1f - %9.1f = %9.1f, DACDiff[%5d - %5d] = %3d. Var/DAC=%9.1f\n",i,
     var1,var0,dvar,dac1,dac0,ddac,dvardac);
  }


}
