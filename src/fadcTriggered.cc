//*****************************
//*  fadcTriggered.cc
//*  created 02/06/2015 G.B. Franklin
//*  Designed for processing Compton FADC triggered data
//*
//*****************************
#include <iostream>
using namespace std;
#include "bankstructure.h"
#include "fadcTriggered.h"

#define TIMEWINDOWlow 566
#define TIMEWINDOWhi 150

#define TIMEWINDOW2low 150
#define TIMEWINDOW2hi 500

#define SUMWINDOWlow 10000
#define SUMWINDOWhi 30000

#define SNAPLIMIT 1     //

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
fadcTriggered::fadcTriggered(){
  //constructor
  return;
}
fadcTriggered::fadcTriggered(textParams* theParamsIn, 
			     comptonStatus* theStatusIn){
  //constructor
  theParams=theParamsIn;  //pointer to parameter handeling class
  theStatus=theStatusIn;  //pointer to the comptonStatushandeling class
  return;
}
void fadcTriggered::newRun(){
  mpsLaserOnCount=0;
  mpsLaserOffCount=0;
  bcmLaserOnSum=0.;
  bcmLaserOffSum=0.;
  ped_value=theParams->getFloat("ped_value");
  channel_calorimeter_PMT=theParams->getFloat("channel_calorimeter_PMT");
  calculate_sum_pedestal=theParams->getInt("calculate_sum_pedestal");
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Define histos
int fadcTriggered::DefineHistos() {
   int fullscale=500000;
  int offset=-1000;
  //  int smallscale=5000;
  //  int smalloffset=-1000;

  TDirectory *topdir= gDirectory;
  TDirectory *triggeredHistos=topdir->mkdir("triggeredHistos");
  triggeredHistos->cd();
  hTrig_numSums=new TH1F("hTrig_numSums","Number Summed FADC Triggers in MPS",
			   100,0,100);
  hTrig_numSamples=new TH1F("hTrig_numSamples",
			    "Number of WF samples in MPS",100,0,100);
  //
  // histograms of integrated triggered pulses
  hTrig_sums_All=new TH1F("hTrig_sums_All","Sums of All Triggered Pulses",
			 5000,offset,fullscale+offset);

  hTrig_sums_laserOn=new TH1F("hTrig_sums_laserOn",
			      "Sums of Laser On Triggered Pulses: Laser On",
			      5000,offset,fullscale+offset);
  hTrig_sums_laserOff=new TH1F("hTrig_sums_laserOff",
			      "Sums of Laser On Triggered Pulses: Laser Off",
			      5000,offset,fullscale+offset);
  hTrig_sums[0]=new TH1F("hTrig_sums_L_N","Sums for (Laser Left)&(Helicty -1)",
			 5000,offset,fullscale+offset);

  hTrig_sums[1]=new TH1F("hTrig_sums_L_P","Sums for (Laser Left)&(Helicty +1)",
			 5000,offset,fullscale+offset);
  hTrig_sums[2]=new TH1F("hTrig_sums_R_N","Sums for (Laser Right)&(Helicty -1)",
			 5000,offset,fullscale+offset);
  hTrig_sums[3]=new TH1F("hTrig_sums_R_P","Sums for (Laser Right)&(Helicty +1)",
			 5000,offset,fullscale+offset);
  //
  // Histograms that count triggers sorted by Laser State
  // MODIFIED to hMPS stuff 6/18/15 gbf
  //
  //
  //normalized histos
  hNorm_Trigs_Scaler=new TH1F("hNorm_Trigs_Scaler",
				 "Num Triggers per MPS vs Spin State",10,0,10);
  labelSpinSortedHistos(hNorm_Trigs_Scaler);

  hNorm_sums_subtracted=new TH1F("hNorm_sums_subtracted",
			      "Background Subtracted Triggered Pulses",
			      5000,offset,fullscale+offset);
  hNorm_sums_asym=new TH1F("hNorm_sums_asym",
	      "Energy-Dependent Asymmetry",
			      5000,offset,fullscale+offset);
  //
  //waveform snapshop histos
  hTrig_wf=new TH1F("hTrig_wf","Sum of Sampled Waveforms",500,0,500);
  hTrig_wf_Ped=new TH1F("hTrig_wf_Ped","Early Channels in Snapshots",
			3000,1000,4000);
  topdir->cd();
  return 0;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcTriggered::DefineTriggeredTree(){
  // data output to tree for each Compton trigger (summed pulses)
  triggerWiseTree=new TTree("triggerwise","Pulse-wise triggered data");
  triggerWiseTree->Branch("sum",&sumVal,0);
  triggerWiseTree->Branch("sumClock",&sumClock,0);
  triggerWiseTree->Branch("sumPedestal",&sumPedestal,0);
  triggerWiseTree->Branch("sumIsRandom",&sumIsRandom,0);
  //now add on variables from comptonStatus 
  theStatus->DefineStatusBranches(triggerWiseTree);
  //
  // data for sampled waveforms
  snapshotsTree=new TTree("snapshots","sampled snapshops");
  snapshotsTree->Branch("randomTime",&randomTime,"randomTime/I");
  snapshotsTree->Branch("numSamples",&NumSamples,"numSamples/I");
  snapshotsTree->Branch("snapClock",&snapshotClock,"snapClock/I");
  snapshotsTree->Branch("snap",&snapshot,"snapshot[numSamples]/F");
  //now add on variables from comptonStatus 
  theStatus->DefineStatusBranches(snapshotsTree);
  //
  // now special tree for MiniMegan pulser data
  pulserWiseTree=new TTree("pulserwise","MiniMegan pulser-wise data");
  pulserWiseTree->Branch("p1",&MMpulse[0],"p1/F");
  pulserWiseTree->Branch("p2",&MMpulse[1],"p2/F");
  pulserWiseTree->Branch("p3",&MMpulse[2],"p3/F");
  pulserWiseTree->Branch("p4",&MMpulse[3],"p4/F");
  pulserWiseTree->Branch("varIndex",&MMVarDacIndex,"varIndex/I");
  pulserWiseTree->Branch("varDAC",&MMVarDacSetting,"varDAC/I");
  pulserWiseTree->Branch("synchIndex",&MMSynchIndex,"synchIndex/I");
  pulserWiseTree->Branch("synchIndexClock",&MMSynchIndexClock,"synchIndexClock/I");
  theStatus->DefineStatusBranches(pulserWiseTree);
  return 0;
}
void fadcTriggered::labelSpinSortedHistos(TH1F* histo){
  //helper function to DefineHisos
  //
  TAxis* pAxis= histo->GetXaxis();
  pAxis->SetBinLabel(1,"LN");
  pAxis->SetBinLabel(2,"LP");
  pAxis->SetBinLabel(3,"RN");
  pAxis->SetBinLabel(4,"RP");
  pAxis->SetBinLabel(5,"LNoff");
  pAxis->SetBinLabel(6,"LPoff");
  pAxis->SetBinLabel(7,"RNoff");
  pAxis->SetBinLabel(8,"RPoff");
  pAxis->SetBinLabel(9,"Unk");
  return;
}

int fadcTriggered::DoSummedPulses(vmeauxdata* theVMEauxdata,
				  fadcdata *theFADCdata){
  // histogram triggered data (pre-summed by CODA )
  int chan=channel_calorimeter_PMT;  
  if(theFADCdata->IsSumsValid(chan)){
    /* 10/01/2016 jc2: For Fall running, we take some channels right
     * before the pulse and use them to determine the pedestal. Use
     * that for this pedestal correction.
     */
    int numTriggersAccepted=theFADCdata->GetSumsNumberTriggersSummed(chan);
    int numRandomsAccepted=theFADCdata->GetSumsNumberRandomsSummed(chan);
    int numInSum=theFADCdata->GetNumberSamplesSummed(chan);
    //Float_t PedCorrection = 0;
    sumPedestal = 0;
    int numInPedSum = 0;

    int crlVersion = theFADCdata->GetCRLVersion();
    int enableNewWaveformReadout = theFADCdata->GetWaveformReadoutVersion();
    int sumSign = 1;
    bool calculatePed = false;
    /* (jc2) Early versions of the CRL (before version 3) corrected the
     * pedestal on-the-fly (meaning, it's in the CODA file). We want to undo
     * this for those versions. So the sign is positive in this case.
     * For CRL versions >= 3 the sign should be negative. */
    if(crlVersion>=3 && enableNewWaveformReadout) {
      if(calculate_sum_pedestal) {
        numInPedSum = theFADCdata->GetNumberPreSamplesSummed(chan);
        calculatePed = true;
      } else {
        sumPedestal = numInSum*ped_value;
      }
      sumSign = -1;
    } else {
      // We want to undo the pedestal correction that was done by the
      // old CRL. Instead, we want to use the one supplied in the
      // compmon.params file.
      sumPedestal = numInSum*ped_value -
        theFADCdata->GetSumsPedestalSubtracted(chan);
    }

    //int numTriggers=theVMEauxdata->GetTriggerScaler();
    mpsCount=theStatus->GetCountMPS();
    helicityState=theStatus->GetHelicityState();
    laserState=theStatus->GetLaserState();
    float bcm=theStatus->GetCalibratedBCM();
    bool beamOn= (theStatus->GetBeamState()==BEAM_ON);
    int sort=theStatus->GetCombinedSpinState();  //sortlaser polarization and beam helicity, etc.
    //prepare to ignore LASER_UNKNOWN data
    bool laserOn= (laserState==LASER_RIGHT||laserState==LASER_LEFT);
    bool laserOff= (laserState==LASER_RIGHTOFF||laserState==LASER_LEFTOFF);

    //following sorted fills are only for BEAM ON condiation
    // (note jc2 11/03/2016: Removed this beam on requirement)
    if(beamOn||true){
      if(laserOn){
        mpsLaserOnCount++;
        bcmLaserOnSum+=bcm;
      }else if (laserOff){
        mpsLaserOffCount++;
        bcmLaserOffSum+=bcm;
      }
      hTrig_numSums->Fill(numTriggersAccepted);
      sumIsRandom = false;
      for(int i=0; i<numTriggersAccepted; i++){
        if(calculatePed) {
          sumPedestal = numInSum*
            (theFADCdata->GetPreSums(chan,i)/double(numInPedSum));
        }
        sumVal = sumSign*theFADCdata->GetSums(chan,i)+sumPedestal;  //move to roottree slot
        sumClock = theFADCdata->GetSumsClock(chan,i);
        //std::cout << "sumPedestal: " << sumPedestal
        //  << ", sumVal: " << sumVal
        //  << ", numTriggersAccepted: " << numTriggersAccepted
        //  << std::endl;
        hTrig_sums_All->Fill(sumVal);
        if(laserOn) hTrig_sums_laserOn->Fill(sumVal);
        if(laserOff) hTrig_sums_laserOff->Fill(sumVal);
        if(sort>=0 && sort<4)  hTrig_sums[sort]->Fill(sumVal);
        triggerWiseTree->Fill();
      }
    }
    // sort pulses for pulserWiseTree if MiniMegan pulser running
    MMVarDacIndex=theFADCdata->GetSumsPulserIndex();//pulser setting index
    MMVarDacSetting=theVMEauxdata->GetDACSetting();//pulser DAC setting
    MMSynchIndex=-1;
    int synchIndex=-1;
    if(numTriggersAccepted>4){
      //Synch Index should be set for one of every 4.  Find the first one
      //ignore the first one in case it overlaps an integration period start
      int bitPattern=0;
      for(int i=1; i<5 ; i++){
        bitPattern= (bitPattern<<1) | theFADCdata->GetSumsPulserSynch(0,i);
      }
      if(bitPattern==0x8) {
        synchIndex=1;
      }else if (bitPattern==0x4){
        synchIndex=2;
      }else if (bitPattern==0x2){
        synchIndex=3;
      }else if (bitPattern==0x1){
        synchIndex=4;
      }
    }
    // if(synchIndex>=0){
    //   for(int i=0; i<numTriggersAccepted; i++){
    // 	MMpulse[ (i-synchIndex+4)%4 ]=
    // 	  theFADCdata->GetSums(chan,i)+PedCorrection;
    // 	if(i%4==3){
    // 	  if(mpsCount>2) pulserWiseTree->Fill(); //Ignore first two MPSs
    // 	}
    //   }
    // }
    int indexPulser;
    bool goodPulserSet=true;
    indexPulser=0;
    MMSynchIndex=synchIndex;
    if(synchIndex>=0){
      int tmpBit = 0;
      for(int i=synchIndex; i<numTriggersAccepted-1; i++){
        tmpBit = (tmpBit<<1)|(theFADCdata->GetSumsPulserSynch(0,i)&&0x1);
        if(indexPulser==0 &&(!theFADCdata->GetSumsPulserSynch(0,i))){
          printf("Error, Missing MM Pulser synch mps=%6d trigger=%d NumTrig %3d\n",
              //mpsCount,i,numTriggersAccepted);
              theStatus->GetMPSCoda(),i,numTriggersAccepted);
          std::cout << "i: " << i << ", tmpBit: ";
          for(int bitI = 31; bitI >= 0; bitI--) {
            if(bitI%8 == 7)
              std::cout << " ";
            std::cout << ((tmpBit>>bitI)&0x1);
          }
          std::cout << ", " << theFADCdata->GetSumsPulserSynch(0,i) << endl;
          goodPulserSet=false;
        }
        if(indexPulser!=0 &&( theFADCdata->GetSumsPulserSynch(0,i))){
          //printf("Error,   Out of place MM Pulser synch mps=%6d trigger=%d\n",mpsCount,i);
          printf("Error,   Out of place MM Pulser synch mps=%6d trigger=%d\n",theStatus->GetMPSCoda(),i);
          goodPulserSet=false;
        }

        indexPulser++;
        if(indexPulser>3)indexPulser=0;
      }
      if(goodPulserSet){
        //output good sets to Tree (ignore first mps data)
        indexPulser=0;
        MMSynchIndexClock = theFADCdata->GetSumsClock(chan,synchIndex);
        for(int i=synchIndex; i<numTriggersAccepted-1; i++){
          if(calculatePed) {
            sumPedestal = numInSum*
              (theFADCdata->GetPreSums(chan,i)/double(numInPedSum));
          }
          MMpulse[indexPulser ]=
            sumSign*theFADCdata->GetSums(chan,i)+sumPedestal;
          if(indexPulser==3 &&mpsCount>0){
            pulserWiseTree->Fill();
          }
          indexPulser++;
          if(indexPulser>3)indexPulser=0;
        }
      // } else {
      // 	//debug MPS data with bad syncs
      // 	int sum;
      // 	int synch;
      // 	int triggerClock;
      // 	int clockLast=0;
      // 	printf("%5d MPS Bad\n",mpsCount);
      // 	for(int i=synchIndex; i<numTriggersAccepted-1; i++){
      // 	  sum=theFADCdata->GetSums(chan,i)+PedCorrection;
      // 	  synch=theFADCdata->GetSumsPulserSynch(0,i);
      // 	  triggerClock=theFADCdata->GetSumsClock(0,i);
      // 	  printf("%3d  %8d  %2d   %8d  %d \n",
      // 		 i,sum,synch,triggerClock,triggerClock-clockLast);
      // 	  clockLast=triggerClock;
      // 	}
      }
    }
    // Now store the randoms
    sumIsRandom = true;
    for(int i=0; i<numRandomsAccepted; i++){
      sumPedestal = numInSum*(theFADCdata->GetRandomPreSums(chan,i)/double(numInPedSum));
      sumVal = sumSign*theFADCdata->GetRandomSums(chan,i)+sumPedestal;
      triggerWiseTree->Fill();
    }
  }else{
    printf("Invalid Summed Triggered Data\n");
  }
  return 0; 
}
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcTriggered::DoNormalizedHistos(){
  // calculates normalized histos (from triggered data)
  //really stupid way to copy histogram.  (Must be a better way)
  //  hNorm_Trigs_Scaler->Add(hTrig_Trigs_Scaler,hTrig_Trigs_Scaler,1.,0.);
  //normalizer per MPS for each laser state
  // hNorm_Trigs_Scaler->Divide(hM_SpinState_BeamOn);

  //
  if(bcmLaserOffSum>0.0){
    Double_t C1=1.0;
    //Scaling by bcm has to also be adjusted by scaler # triggers to #trigger accepts
    //We don't have that in test run, so we'll hardwire a quick scaling of data
    //past the compton edge
    int binlo=hTrig_sums_laserOn->FindBin(10000.);
    int binhi=hTrig_sums_laserOn->FindBin(20000.);
    Double_t sumOn=hTrig_sums_laserOn->Integral(binlo,binhi);
    Double_t sumOff=hTrig_sums_laserOff->Integral(binlo,binhi);
    Double_t C2= -sumOn/sumOff;
    //    Double_t C2= -bcmLaserOnSum/(bcmLaserOffSum);  //note negative sign
    //Not clear if bcm sums are good.  Lets try just counting MPSs
    //Double_t C2= -float(mpsLaserOnCount)/float(mpsLaserOffCount);
    hNorm_sums_subtracted->Add(hTrig_sums_laserOn,hTrig_sums_laserOff,C1,C2);
    // build asymmetry histogram
    //start by building denominator
    hNorm_sums_asym->Add(hTrig_sums[1],hTrig_sums[2],1.,1.);
    hNorm_sums_asym->Add(hTrig_sums[0],-1.);  //note negative signs
    hNorm_sums_asym->Add(hTrig_sums[3],-1.);
    //divide by background subtracted sums
    hNorm_sums_asym->Divide(hNorm_sums_subtracted);
  }
  return 0;
}

//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcTriggered::DoSampledWaveforms(fadcdata *theFADCdata){
  // histogram triggered data (presummed by CODA )
   if(theFADCdata->IsSamplesValid(0)){
    int chan=0;   //only FADC channel zero for now
 //number FADC channel read out
    int TotalSamples=theFADCdata->GetNumberSamples(chan);
 //# samples stored for each sample trigger
    randomTime=theFADCdata->IsRandomTimes(chan); 
    NumSamples=theFADCdata->GetSamplesPerEvent(chan);  //root tree variable
    //# of waveforms stored
    int NumEvents= theFADCdata->GetNumberEvents(chan);
    if(TotalSamples!= NumSamples*NumEvents){
      NumEvents=TotalSamples/NumSamples;   //fix if data overflow
      printf("Total Samples doesn't agree with number of events\n");
    }
    mpsCount=theStatus->GetCountMPS();
    helicityState=theStatus->GetHelicityState();
    laserState=theStatus->GetLaserState();
    hTrig_numSamples->Fill(NumEvents);
    pulsestructure Pulse;
    int* bits;
    int* data;
    for(int event=0; event<NumEvents; event++){
      Pulse=theFADCdata->GetPulse(chan,event);
      bits=Pulse.UserBits;
      snapshotClock=Pulse.Clock;  //3/14/2016  gbf
      data=Pulse.Data;
      for(int i=0; i<NumSamples&&i<1002; i++){
 	snapshot[i]=data[i];
	hTrig_wf->Fill(i,data[i]);  //use weighted fill
      }
      for(int i=0; i<50; i++){
	hTrig_wf_Ped->Fill(data[i]);  //histogram early channels (Pedestal info)
      }
      snapshotsTree->Fill();
    }
  }else{
    printf("Invalid Waveform Data\n");
  }
  return 0;   
}

