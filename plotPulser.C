#include "../LED/MiniMegan.h"
// Last edit: 20190301 Juan Carlos Cornejo
//   Make it work with ROOT 6 and push hard coded values to the top
//   Also fix correct sequence (Delta,Both,Nothing,Var,....)
// Edited 2016-02-20 Larisa Thorne
// Want gr to be global
// ftp://root.cern.ch/root/doc/ROOTUsersGuideHTML/ch02s06.html

// Set the plot limits here
//const Int_t kVarDACCE = 32661; // CRex setting
const Bool_t kUseStateCuts = 0;
const Int_t kBeamState = 1;
const Int_t kLaserState = 1;
const Int_t kVarDACCE = 33355; // PRex setting
const Float_t kDeltaLow  = -10000*0-2000;
const Float_t kDeltaHigh =  10000*0+5000*0+11000;
const Float_t kVarLow    = -10000*0-3000;
const Float_t kVarHigh   =  50000*0+50000;
const Float_t kDiffLow   =   1000*0-8000;
const Float_t kDiffHigh  =  50000*0+8000*0+24000;
const Float_t kBothLow   = -10000;
const Float_t kBothHigh  = 1.30*kVarHigh;
const Float_t kNothingLow = -4000;
const Float_t kNothingHigh = 30000;

const Float_t kCancelLow = -10e3;
const Float_t kCancelHigh = 10e3;


const Float_t kBinRes = 10*0+50; // srau

//const Int_t   kDeltaBins =  1000*0+500;
//const Int_t   kVarBins   =  6000;
//const Int_t   kDiffBins  = 60000*0+500;
//const Int_t   kBothBins  = kVarBins*1.3;
//const Int_t  kNothingBins = 10000*0+500;
const Int_t   kDeltaBins =  (kDeltaHigh-kDeltaLow)/kBinRes;
const Int_t   kVarBins =  (kVarHigh-kVarLow)/kBinRes;
const Int_t   kDiffBins =  (kDiffHigh-kDiffLow)/kBinRes;
const Int_t   kBothBins =  (kBothHigh-kBothLow)/kBinRes;
const Int_t   kNothingBins =  (kNothingHigh-kNothingLow)/kBinRes;
const Int_t   kCancelBins =  (kCancelHigh-kCancelLow)/kBinRes;


const Int_t kSkipAfterTrip = 5*120; //


TGraphErrors *gr;

TH1F *hVarCE;

const int nsteps=200;
MM::PulseHistos histo[nsteps];

// -------- Main program ----------------------

//void plotPulser(Int_t Wanted_Dac_Index=13){
void plotPulser(Int_t runnum,Int_t Wanted_Dac_Index=13){
  gStyle->SetOptFit(1);
  MM::init(runnum);
  // TFile *rootfile=new TFile("lnkCompMon.output");
  //TFile *rootfile = gROOT->GetFile();
  //TTree *t1 = (TTree*)rootfile->Get("pulserwise");
  TChain *t1 = new TChain("pulserwise");
  Int_t nfiles;
  TString filesPre;
  if(runnum>70e3) {
    ifstream in;
    in.open(Form("%s/rungroup_%d.txt",getenv("COMP_RUNGROUPS"),runnum),std::ios::in);
    if(!in.is_open())
      return;
    int rn;
    while(in >> rn && !in.eof()) {
      int new_files = 0;
      filesPre=Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),rn);
      new_files += t1->Add(filesPre+".root");
      new_files += t1->Add(filesPre+"_*.root");
      if(new_files<=0) {
        std::cerr << "Looked for files under: " << filesPre+".root" << std::endl;
        std::cerr << "Found no files to plot!" << std::endl;
      }
    }
  } else {
    filesPre = Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),runnum);
    nfiles = t1->Add(filesPre+".root");
    nfiles += t1->Add(filesPre+"_*.root");
    if(nfiles<=0) {
      std::cerr << "Looked for files under: " << filesPre+".root" << std::endl;
      std::cerr << "Found no files to plot!" << std::endl;
      return;
    }
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
  int beamState; //
  int laserState; //
  Float_t delta_contribution;
  Float_t all_delta;
  //const int nsteps = 200;
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
  if(kUseStateCuts) {
    t1->SetBranchAddress("beamState",&beamState);
    t1->SetBranchAddress("laserState",&laserState);
  }

  char his_title[128];
  char his_label[128];

  //TH1F* h_delta[nsteps];   // Delta minus background
  //TH1F* h_var[nsteps];     // Var minus bckground
  //TH1F* h_diff[nsteps];    // Both - var
  //TH1F *h_nothing[nsteps];  // The pedestal (LED off)
  //TH1F *h_cancel[nsteps];  // Check that things cancel out
  std::cout << "Bins::"
      << "  kDeltaBins: " << kDeltaBins
      << "  kVarBins: " << kVarBins
      << "  kDiffBins: " << kDiffBins
      << "  kBothBins: " << kBothBins
      << "  kNothingBins: " << kNothingBins
      << std::endl;
  for(Int_t i=0;i<nsteps;i++){
    sprintf(his_title,"DAC Setting_%03d",i);
    //printf(his_title,"DAC Setting_%03d",i);
    sprintf(his_label,"h_delta_%03d",i);
    histo[i].delta=new TH1F(his_label,his_title,kDeltaBins,kDeltaLow,kDeltaHigh);
    sprintf(his_label,"h_var_%03d",i);
    histo[i].var=new TH1F(his_label,his_title,kVarBins,kVarLow,kVarHigh);
    sprintf(his_label,"h_diff_%03d",i);
    histo[i].diff=new TH1F(his_label,his_title,kDiffBins,kDiffLow,kDiffHigh);
    sprintf(his_label,"h_nothing_%03d",i);
    histo[i].nothing=new TH1F(his_label,his_title,kNothingBins,kNothingLow,kNothingHigh);
    sprintf(his_label,"h_cancel_%03d",i);
    histo[i].cancel=new TH1F(his_label,his_title,kCancelBins,kCancelLow,kCancelHigh);
    sprintf(his_label,"h_both_%03d",i);
    histo[i].both=new TH1F(his_label,his_title,kBothBins,kBothLow,kBothHigh);

  // (20190305 Cornejo): New, make these use the same limits as the new limits
  // defined at the top of this script
  //float xMax=30000;
  //float xMin=-4000;
  //TH1F*h1 = new TH1F("Var + Delta","Var + Delta",10000,xMin,xMax);
  //TH1F*h2 = new TH1F("Delta","Delta",10000,xMin,xMax);
  //TH1F*h3 = new TH1F("Var","Var",10000,xMin,xMax);
  //TH1F*h4 = new TH1F("BG","BG",10000,xMin,xMax);  
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
  //hVarCE = new TH1F("hVarCE","Variable @ CE", kVarBins,kVarLow,kVarHigh);


  Int_t nevents=(Int_t)t1->GetEntries(); // Found form in sample 
  printf("nevents: %d\n", nevents); // Must have nonzero value to proceed!

  Int_t buffer[nsteps];
  for(Int_t i=0;i<nsteps;i++){
    buffer[i]=i;
    settingDAC[i]=-1;
  }

  //float pVar,pDelta,pBoth,pNothing;
  float &pVar     = p1;
  float &pDelta   = p2;
  float &pBoth    = p3;
  float &pNothing = p4;

  Int_t start_entry = 20;
  Int_t nevents_avail = nevents-start_entry;
  Int_t nevents_used = 0;
  //for(Int_t i=20*0+10e3;i<nevents;i++){
  int last_beam_trip = -1e4;
  for(Int_t i=start_entry;i<nevents;i++){
    t1->GetEntry(i);

    if(kUseStateCuts) {
      if((kBeamState!=-1) && (beamState!=kBeamState)) {
        last_beam_trip = i;
        continue;
      }

      if(i-last_beam_trip<kSkipAfterTrip)
        continue;

      if((kLaserState!=-1) && (laserState!=kLaserState))
        continue;
    }
    nevents_used++;

    //if(i==0||i<100){
    //  cout<<"Dac Index is    "<<varIndex<<endl;
    // }
    /*
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
    */

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
      histo[varIndex].delta->Fill(pDelta-pNothing);
      histo[varIndex].var->Fill(pVar-pNothing);
      histo[varIndex].diff->Fill(pBoth-pVar);
      histo[varIndex].nothing->Fill(pNothing);
      histo[varIndex].cancel->Fill((pBoth-pVar)-(pDelta-pNothing) - MM::crossTalk(pVar-pNothing));
      histo[varIndex].both->Fill(pBoth-pNothing);
      if(settingDAC[varIndex]==-1) {
        settingDAC[varIndex]=varDAC;
      }
      //if(varDAC==kVarDACCE) {
        //hVarCE->Fill(pVar-pNothing);
      //}
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
    entries[i]=histo[i].diff->GetEntries();
    if(entries[i]==0) break;
    nonzeroSteps++;
    c2->cd(i+1);
    histo[i].diff->Draw();
    // Just get mean from histogram
    meanDiff[i]=histo[i].diff->GetMean();
    meanVar[i]=histo[i].var->GetMean();  
        rms[i]=histo[i].diff->GetRMS();
    sqrtentries[i]=sqrt(entries[i]);
    error[i]=rms[i]/sqrtentries[i];
    errorX[i]=0.;   //just to keep TGraph happy


    // Alternative method, fit with a gausian
    // Fit once to setup the fit
    double tmp_min,tmp_max;
    tmp_min = meanDiff[i]-histo[i].diff->GetRMS()*3.;
    tmp_max = meanDiff[i]+histo[i].diff->GetRMS()*3.;
    TFitResultPtr fr = histo[i].diff->Fit("gaus","SQ","",tmp_min,tmp_max);
    // Fit a second time for better results and to record fit results
    if(fr==0) {
      tmp_min = fr->Parameter(1)-fr->Parameter(2)*5.;
      tmp_max = fr->Parameter(1)+fr->Parameter(2)*5.;
    }
    fr = histo[i].diff->Fit("gaus","SQ","",tmp_min,tmp_max);
    if(fr==0) {
      meanDiff[i] = fr->Parameter(1);
      error[i] = fr->ParError(1);
    }

    printf(" %4d Var %12.1f  Delta %10.1f  Diff %10.3f p/m %10.2f [DAC=%5d]\n",
 	   i, histo[i].var->GetMean(),histo[i].delta->GetMean(), meanDiff[i],error[i],settingDAC[i]);
  }
  if(nonzeroSteps<=0)
    return;
  int draw_count = 0;
  float draw_every=float(nonzeroSteps)/24.;
  for(int i = 0; i < 24; i++ ) {
    int j = i*draw_every;
    c2->cd(i+1);
    histo[j].diff->Draw();
  }
  c2->cd(25);
  histo[nonzeroSteps-1].diff->Draw();

  TCanvas* c3 = new TCanvas("canvLinearity","Linearity");

  //TGraph *gr = new TGraphErrors(nonzeroSteps, meanVar, meanDiff, errorX, error);
  gr = new TGraphErrors(nonzeroSteps, meanVar, meanDiff, errorX, error);

  if(kUseStateCuts) {
    gr->SetTitle(Form("MM Pulser Linearity Plot (B=%s,L=%s)",
          kBeamState==1?"ON":"OFF",kLaserState==1?"ON":"OFF"));
  } else {
    gr->SetTitle("MM Pulser Linearity Plot");
  }
  gr->GetXaxis()->SetTitle("VAR LED pulses (srau)");
  gr->GetYaxis()->SetTitle("(VAR+Delta)-VAR (srau)");
  gr->Draw("AP");
  if(kUseStateCuts) {
    TString oname = Form("lin_beam/lin_%d_B%dB%d",runnum,
          kBeamState,kLaserState==1?1:0);
    c3->SaveAs(oname+".C");
    c3->SaveAs(oname+".png");
  }

  /* You can store any named object in gROOT->GetListOfSpecials() with 
gROOT->GetListOfSpecials()->Add(myobject); 
and retrieve it later with gROOT->FindObject. */
  //gROOT->GetListOfSpecials()->Add(gr);
  //printf(gROOT->FindObject("gr"));

  printf("Done. nevents=%10d, used=%10d (%3.2f%%)\n\n",nevents,nevents_used,100.*Float_t(nevents_used)/nevents_avail);

  float var0,var1,dvar,dvardac;
  int dac0,dac1,ddac;
  for(Int_t i=1;i<nsteps*0;i++){
    entries[i]=histo[i].diff->GetEntries();
    if(entries[i]==0) break;
    var1=histo[i-1].var->GetMean();
    var0=histo[i].var->GetMean();
    dac1=settingDAC[i-1];
    dac0=settingDAC[i];
    dvar=var0-var1;
    ddac=dac0-dac1;
    dvardac = dvar/float(ddac);
    printf(" %4d Var: %9.1f - %9.1f = %9.1f, DACDiff[%5d - %5d] = %3d. Var/DAC=%9.1f\n",i,
     var1,var0,dvar,dac1,dac0,ddac,dvardac);
  }
}
