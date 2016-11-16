// Last edit: Larisa Thorne, 2016-02-22
// Description: 

enum combinedSpinFlag_t{SPIN_LN,SPIN_LP,SPIN_RN,SPIN_RP,
	  SPIN_LNOFF,SPIN_LPOFF,SPIN_RNOFF,SPIN_RPOFF,SPIN_UNKNOWN};

TString getSpinLabel(combinedSpinFlag_t spinState){

  switch(spinState){
  case SPIN_LN:
    return "LN";
    break;
  case SPIN_LP:
    return "LP";
    break;
   case SPIN_RN:
    return "RN";
    break;
  case SPIN_RP:
    return "RP";
    break;
  case SPIN_LNOFF:
    return "LNOFF";
    break;
  case SPIN_LPOFF:
    return "LPOFF";
    break;
   case SPIN_RNOFF:
    return "RNOFF";
    break;
  case SPIN_RPOFF:
    return "RPOFF";
    break;
  }

  return "UNKNOWN";
}


// --- Main program -----------------------------------------------

void plotAcc0(){
  // Reads mpswise ROOT tree
  // Set up to read ROOT file

  TFile *rootfile = new TFile("lnkCompMon.output");
  TTree *t1 = (TTree*)rootfile->Get("mpswise");

  // Tree variables
  int mpsCoda;
  int laserState;
  int helicityState;
  combinedSpinFlag_t combinedSpinState;
  int beamState;
  Double_t Acc0;

  TH1F* hAcc0[9]; // Creates array of hAcc0 histograms, one for each CSS
  TString hName;
  TString hTitle; 
  TH1F* hMPS = new TH1F("hMPS", "MPS vs Combined Spin State", 9, 0, 9);
  TH1F* hAcc0Sum = new TH1F("hAcc0Sum", "Acc0 Combined Spin State", 9, 0, 9);

  // Create histogram for each spin/helicity/laser state combination:
  for (Int_t i=0; i<9; i++){
    hName = "hAcc0_" + getSpinLabel(i);
    hTitle = "Acc0 for Spin State " + getSpinLabel(i);
    hAcc0[i] = new TH1F(hName, hTitle, 1000, 0., 2E7);
  }

  t1->SetBranchAddress("mpsCoda", &mpsCoda);
  t1->SetBranchAddress("laserState", &laserState);
  t1->SetBranchAddress("helicityState", &helicityState);
  t1->SetBranchAddress("combinedSpinState", &combinedSpinState);
  t1->SetBranchAddress("beamState", &beamState);
  t1->SetBranchAddress("Acc0", &Acc0);

  Int_t ntriggers = (Int_t)t1->GetEntries();
  printf("Number of triggers in file = %d\n", ntriggers);

  // Fill Acc0 histograms:
  for (int i=0; i<ntriggers; i++){
    t1->GetEntry(i);
    hMPS->Fill(combinedSpinState);
    hAcc0Sum->Fill(combinedSpinState,Acc0);
    hAcc0[combinedSpinState]->Fill(Acc0);
  }
  TH1F *hAcc0Normalized = (TH1F*)hAcc0Sum->Clone("hAcc0Normalized");
  hAcc0Normalized->Divide(hMPS);


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
  delete runHisto;



  // Plot combined spin state event count information:
  TString c1Title;
  c1Title.Form("Run %d: Comparison of Combined Spin State", runNum);
  TCanvas *c1 = new TCanvas("c1", c1Title, 30, 30, 600, 500);
  c1 -> Divide(1,2);

  c1->cd(1);
  hMPS->Draw();

  c1->cd(2);
  hAcc0Normalized->Draw();
  

  // Plot hAcc0 for each combined spin state:
  // Note that some plots will probably be empty. Compare to empty bins in c1
  TCanvas *c2 = new TCanvas("c2", "hAcc0s for all combined spin states", 30, 30, 800, 600);
  c2->Divide(3,3);

  c2->cd(1);
  hAcc0[0]->Draw();
  c2->cd(2);
  hAcc0[1]->Draw();
  c2->cd(3);
  hAcc0[2]->Draw();
  c2->cd(4);
  hAcc0[3]->Draw();
  c2->cd(5);
  hAcc0[4]->Draw();
  c2->cd(6);
  hAcc0[5]->Draw();
  c2->cd(7);
  hAcc0[6]->Draw();
  c2->cd(8);
  hAcc0[7]->Draw();
  c2->cd(9);
  hAcc0[8]->Draw();
  

  printf("Done.\n");

}
