void plotSnaps(){
  // read snapshot info from tree and plot'
  //setyo ti read riit fuke
  //  TFile *rootfile=new TFile("lnkCompMon.output");
  TFile *rootfile=gROOT->GetFile();
  TTree *t1=(TTree*)rootfile->Get("snapshots");
  // Tree variables
  int numSamples;
  int mpsCoda;
  Float_t samples[1000];
  int laserState;
  int beamState;
  int randomTime;
  int triggerTypeRequest=2;  //0 for triggered, 1 for random, 2 for either 

  t1->SetBranchAddress("numSamples",&numSamples);
  t1->SetBranchAddress("randomTime",&randomTime);
  t1->SetBranchAddress("snap",&samples);
  t1->SetBranchAddress("mpsCoda",&mpsCoda);
  t1->SetBranchAddress("laserState",&laserState);
  t1->SetBranchAddress("beamState",&beamState);

  int gotoMPS=0;
  int skipMPS=0;
  int skipMPScounter=0;
  int oldMPS=0;

  //
  // setup for histogram
  int MaxSamples=1000;
  TString cresponse;

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

  TString c1Title;
  c1Title.Form("Run %d: Snapshots", runNum);
  TCanvas* c1=new TCanvas("c1", c1Title);
  TH1F* hsnap==NULL;
  Int_t ntriggers=(Int_t)t1->GetEntries();
  printf("Number of triggers in file = %d\n",ntriggers);
  for(int i=0; i<ntriggers; i++){
    t1->GetEntry(i);
    if(oldMPS!=mpsCoda){
      skipMPScounter++;
      oldMPS=mpsCoda;
    }
    if(skipMPScounter<skipMPS) continue;
    if(gotoMPS>0 && mpsCoda<gotoMPS) continue;
    skipMPScounter=0;
    gotoMPS=0;
    skipMPS=0;
    if(hsnap==NULL&&numSamples>0){
        hsnap=new TH1F("hsnap","snapshot",numSamples,0,numSamples);
    }
    if(hsnap!=NULL){
      if(randomTime==0){
 	if(triggerTypeRequest==0||triggerTypeRequest==2){
	  for(int j=0; j<numSamples; j++){
	    hsnap->Fill(j,samples[j]);
	  }
	  hsnap->SetTitle("Triggered Snapshot");
	  hsnap->Draw();
	  c1->Update();
	}else{
	  continue;
	}
      }else{
 	if(triggerTypeRequest==1||triggerTypeRequest==2){
	  for(int j=0; j<numSamples; j++){
	    hsnap->Fill(j,samples[j]);
	  }
	  hsnap->SetTitle("Random Time Snapshot");
	  hsnap->Draw();
	  c1->Update();
	}else{
	  continue;  //go to next event without plotting
	}
      }
    }
    bool wait=true;
    printf("RandomTime %d  MPS  %6d Laser %2d BeamOn %2d >",
	   randomTime, mpsCoda, laserState, beamState);
    while(wait){

      cin >> cresponse;
      cresponse.ToUpper();
      if(cresponse.Index("Q")==0) {
	return;
      }else if (cresponse.Index("N")==0) {
	wait=false;
	triggerTypeRequest=2;   //plot any waveform
      }else if (cresponse.Index("T")==0) {
	wait=false;    //plot triggered waveform
	triggerTypeRequest=0;
      }else if (cresponse.Index("R")==0) {
	wait=false;
	triggerTypeRequest=1;   //plot random time waveform
      }else if(cresponse.Index("G")==0) {
	cin >>cresponse;   //get number MPS to Go To
	gotoMPS=cresponse.Atoi();
	wait=false;
      }else if(cresponse.Index("S")==0) {
	cin >>cresponse;   //get number MPS to Skip
	skipMPS=cresponse.Atoi();
	wait=false;
      }else if (cresponse.Index("?")==0){
	printf("q    Quit \n");
	printf("n    Next sample\n");
	printf("t    Next triggered sample\n");
	printf("r    Next random time sample\n");
	printf("s x  Skip x MPSs\n");
	printf("g x  Goto MPS x\n");
	printf("?    Print commands\n");
	printf(" ?                           >");
      }else{
	printf(" ?                           >");
      }
    }
    hsnap->Reset();
  }
}
