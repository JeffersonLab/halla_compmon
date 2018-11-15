//*****************************
//*  fadcAccums.cc
//*  created 02/06/2015 G.B. Franklin
//*  For processing Compton FADC accumuator data
//*
//*****************************
#include <iostream>
using namespace std;
#include "bankstructure.h"
#include "fadcAccums.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
fadcAccums::fadcAccums(textParams* theParamsIn, comptonStatus* theStatusIn){
  //constructor
  theParams=theParamsIn;  //pointer to parameter handeling class
  theStatus=theStatusIn;  //pointer to status handeling class
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void fadcAccums::newRun(){
  channel_calorimeter_PMT=theParams->getInt("channel_calorimeter_PMT");
  ped_value=theParams->getFloat("ped_value");
  lastSubQuad=-1;    //last subquad record not yet defined
  prevLaserStateEncountered=-1;  //last Laser State not yet defined
  laserStateBeingSummed=-1;
  laserStateThisQuad=LASER_UNKNOWN;  //used within 1 Quad 
  quartetValid=false;  //used to check 4 consecutive helicity windows 
  quartetStable=false;
  countQuadsLaserPeriod=0;
  pointerHistory=0;      //initialize pointer
  for(int i=0; i<NUM_ACCUM_TYPES; i++){
    accSumLastLaserOff[i]=-1.e9;  //keeps last laser-off data
  }
  for(int i=0; i<NUM_HISTORY_LASER_PERIODS; i++){
    countQuadsHistory[i]=0;  //zero out numer of quads summed into each history period
  }
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//Define histos
int fadcAccums::DefineHistos() {
  //  float accumFullScale=2E8;
  //  float scaledRange=accumFullScale; // No scaling for now
  float accumFullScale=theParams->getFloat("accum_histo_range");
  float accumDiffMax=theParams->getFloat("accum_histo_diff_max");
  float accumHistoScaledRange=theParams->getFloat("accum_histo_scaled_range");

  float sumMax   =2*accumHistoScaledRange;
  float sumMin   =-0.2*sumMax;
  float diffMax = accumDiffMax;

  TString label;
  TString title;
  TDirectory *topdir= gDirectory;
  TDirectory *mpsHistos=topdir->mkdir("mpsHistos");
  TDirectory *quartetHistos=topdir->mkdir("quartetHistos");
  TDirectory *laserCycleHistos=topdir->mkdir("laserCycleHistos");
  mpsHistos->cd();
  //
  // Bookkeeping histos
  hM_BCM=new TH1F("hM_BCM","BCM per MPS",1700.,-20.,150.);
  hM_BCM_BeamOn=new TH1F("hM_BCM_BeamOn",
      "BCM per MPS (Beam On)",1700.,-20.,150.);
  hM_BCM_BeamOff=new TH1F("hM_BCM_BeamOff",
      "BCM per MPS (Beam Off)",1700.,-20.,150.);

  hM_BPMSum=new TH1F("hM_BPMSum","BPM Sum per MPS",1700.,-20.,150.);
  hM_BPMSum_BeamOn=new TH1F("hM_BPMSum_BeamOn","BPM Sum per MPS (Beam On)",
			    1700.,-20.,150.);
  hM_BPMSum_BeamOff=new TH1F("hM_BPMSum_BeamOff","BPM Sum per MPS (Beam Off)",
			     1700.,-20.,150.);

  hM_CavPower=new TH1F("hM_Cav_Power","Laser Cavity Power",1050,-1000.,10000.);
  hM_CavPower_LaserOn=new TH1F("hM_Cav_Power_LaserOn",
			       "Laser Cavity Power (laser ON)",
			      1050,-1000.,10000.);
  hM_CavPower_LaserOff=new TH1F("hM_Cav_Power_LaserOff",
				"Laser Cavity Power (Laser Off)",
			      1050,-1000.,10000.);

  // Spin-State sorted
  hM_SpinState=new TH1F("hM_SpinState",
			       "Num MPS  vs Spin State",10,0,10);
  labelSpinSortedHistos(hM_SpinState);
  //
  hM_SpinState_BeamOn=new TH1F("hM_SpinState_BeamOn",
			       "Num MPS with Beam On  vs Spin State",10,0,10);
  labelSpinSortedHistos(hM_SpinState_BeamOn); //set labels of x-axis
  //
  hM_BCM_SummedBySpinState=new TH1F("hM_BCM_SummedBySpinState",
     "BCM summed for each Spin State (Beam On)",10,0,10);
  labelSpinSortedHistos(hM_BCM_SummedBySpinState);  //set labels of x-axis

  hM_Trigs_Accepted=new TH1F("hM_Trig_Accepted",
		        "Num Trig Evnts Accepted vs Spin State (Beam On)",10,0,10);
 labelSpinSortedHistos(hM_Trigs_Accepted);
  //
  hM_Trigs_Scaler=new TH1F("hM_Trigs_Scaler",
				 "Num Triggers vs Spin State (Beam oN)",10,0,10);
 labelSpinSortedHistos(hM_Trigs_Scaler);
 //
  hM_Trigs_Prescaled=new TH1F("hM_Trigs_Prescaled",
				 "Num Prescaled Triggers vs Spin State (Beam On)",10,0,10);
 labelSpinSortedHistos(hM_Trigs_Scaler);

  hM_numSums=new TH1F("hM_numSums","Number Summed FADC Triggers in MPS (Beam On)",
			   100,0,100);

  //
  //Accumulator Histoes
  // catch all acc0 zero histogram
  hM_acc0_everything=new TH1F("hM_acc0_everything",
			      "Accum 0 wide all conditions",
			      500,-100.*accumFullScale,100.*accumFullScale);
  //
  // Beam Off Histos
  hM_acc0_beamOff_LaserOn=new TH1F("hM_acc0_beamOff_LaserOn","Accum 0 Beam Off Laser On",
			   25000,-accumFullScale,accumFullScale);
  hM_acc0_beamOff_laserOff=new TH1F("hM_acc0_beamOff_LaserOff",
				    "Accum 0 Beam Off Laser Off",
				    25000,-accumFullScale,accumFullScale);


  hM_acc4_beamOff_LaserOn=new TH1F("hM_acc4_beamOff_LaserOn","Accum 4 Beam Off Laser On",
			   25000,-accumFullScale,accumFullScale);
  hM_acc4_beamOff_laserOff=new TH1F("hM_acc4_beamOff_LaserOff",
				    "Accum 4 Beam Off Laser Off",
				    25000,-accumFullScale,accumFullScale);
  //MPS-wise Beam On Histos
  for(int accum=0; accum<8; accum++){
    label=Form("hM_acc%d_All",accum);
    title=Form("Accum %d all MPS",accum);
    // include large negative range incase pedestal subtraction not correct
    hM_acc_All[accum]=new TH1F(label,title,
			   25000,-accumFullScale,accumFullScale);

    label=Form("hM_acc%d_LaserOn",accum);
    title=Form("Accum %d  (Laser On) ",accum);
    hM_acc_laserOn[accum]=new TH1F(label,title,
			   25000,sumMin,sumMax);
    label=Form("hM_acc%d_LaserOff",accum);
    title=Form("Accum %d  (Laser Off) ",accum);
    hM_acc_laserOff[accum]=new TH1F(label,title,
				    25000,sumMin,sumMax);
  }
  //quartet  histos
  //Acc 0 Beam-off histos
  quartetHistos->cd();
  hQ_aligned0_beamOff=new TH1F("hQ_aligned0_beamOff",
			       "Acc0 Quad Spin Aligned Beam Off",
			       1000,sumMin,sumMax);
  hQ_anti0_beamOff=new TH1F("hQ_anti0_beamOff",
			       "Acc0 Quad Spin Anti-aligned Beam Off",
			       1000,sumMin,sumMax);

  hQ_sum0_beamOff=new TH1F("hQ_sum0_beamOff","Acc0 Quad Sum Beam Off",
		       1000,sumMin,sumMax);
  hQ_diff0_beamOff=new TH1F("hQ_diff0_beamOff","Acc0 Quad Diff Beam Off",
			    1000,-diffMax,diffMax);
  hQ_Asym0_raw_beamOff=new TH1F("hQ_Asym0_raw_beamOff",
				     "Acc0 Raw Asym BeamOff",
				1000,-0.3,0.3);
  //quartet  histos (Beam On)
  //accumulator 0
  hQ_aligned0_laserOff=new TH1F("hQ_aligned0_laserOff",
			       "Acc0 Quad Spin Aligned Laser Off",
			       1000,sumMin,sumMax);
  hQ_anti0_laserOff=new TH1F("hQ_anti0_laserOff",
			       "Acc0 Quad Spin Anti Aligned Laser Off",
			       1000,sumMin,sumMax);

  hQ_aligned0_laserOn=new TH1F("hQ_aligned0_laserOn",
			       "Acc0 Quad Spin Anti-aligned laser On",
			       1000,sumMin,sumMax);
  hQ_anti0_laserOn=new TH1F("hQ_anti0_laserOn",
			       "Acc0 Quad Spin Anti-aligned laser On",
			       1000,sumMin,sumMax);

  hQ_sum0_laserLeft=new TH1F("hQ_sum0_laserLeft","Acc0 Quad Sum laser left",
		       1000,sumMin,sumMax);
  hQ_sum0_laserOff=new TH1F("hQ_sum0_laserOff","Acc0 Quad Sum laser Off",
		       1000,sumMin,sumMax);
  hQ_diff0_laserLeft=new TH1F("hQ_diff0_laserLeft","Acc0 Quad Diff laser Left",
		       1000,-diffMax,diffMax);
  hQ_diff0_laserOff=new TH1F("hQ_diff0_laserOff","Acc0 Quad Diff Laser Off",
		       1000,-diffMax,diffMax);
  hQ_Asym0_raw_laserLeft=new TH1F("hQ_Asym0_raw_laserLeft",
				     "Acc0 Raw Asym Laser Left",
				   1000,-0.3,0.3);
  hQ_Asym0_raw_laserOff=new TH1F("hQ_Asym0_raw_laserOff",
				    "Acc0 Raw Asym Laser Off",
				 1000,-0.3,0.3);
  hQ_Signal0=new TH1F("hQ_Signal0","Acc0 Background Subtracted",
		      1000,2*sumMin,sumMax);
  hQ_Asym0=new TH1F("hQ_Asym0","Acc0 Asym (background subtracted)",
		    1000,-0.3,0.3);

  //accumulator 4
  float diffMax4=0.2*diffMax;
  float sumMin4=0.01*sumMin;
  float sumMax4=0.01*sumMax;
  hQ_sum4_beamOff=new TH1F("hQ_sum4_beamOff","Acc4 Quad Sum laser left",
			   1000,sumMin4,sumMax4);
  hQ_diff4_beamOff=new TH1F("hQ_diff4_beamOff","Acc4 Quad Diff laser Left",
			    1000,-diffMax4,diffMax4);
  hQ_sum4_laserLeft=new TH1F("hQ_sum4_laserLeft","Acc4 Quad Sum laser left",
			     1000,sumMin4,sumMax4);
  hQ_sum4_laserOff=new TH1F("hQ_sum4_laserOff","Acc4 Quad Sum laser Off",
		       1000,sumMin4,sumMax4);
  hQ_diff4_laserLeft=new TH1F("hQ_diff4_laserLeft","Acc4 Quad Diff laser Left",
		       1000,-diffMax4,diffMax4);
  hQ_diff4_laserOff=new TH1F("hQ_diff4_laserOff","Acc4 Quad Diff Laser Off",
			     1000, -diffMax4,diffMax4);
  hQ_Asym4_raw_laserLeft=new TH1F("hQ_Asym4_raw_laserLeft",
				     "Acc4 Raw Asym Laser Left",
				   1000,-0.5,0.5);
  hQ_Asym4_raw_laserOff=new TH1F("hQ_Asym4_raw_laserOff",
				    "Acc4 Raw Asym Laser Off",
				   1000,-0.5,0.5);
  hQ_Signal4=new TH1F("hQ_Signal4","Acc4 Background Subtracted",
		      1000,sumMin,sumMax);
  hQ_Asym4=new TH1F("hQ_Asym4","Acc4 Asym (background subtracted)",
				   1000,-0.5,0.5);
  //beam
  hQ_BCM_laserOff_P=new TH1F("hQ_BCM_laserOff_P",
			    "Beam Charge Positive Hel. LaserOff",
			    6000,-10.,50.);
  hQ_BCM_laserOff_N=new TH1F("hQ_BCM_laserOff_N",
			    "Beam Charge Negative Hel. Laser Off",
			    6000,-10.,50.);

  hQ_BCM_laserOn_P=new TH1F("hQ_BCM_laserOn_P","Beam Charge Positive Hel. Laser On",
			     6000,-10.,50.);
  hQ_BCM_laserOn_N=new TH1F("hQ_BCM_laserOn_N","Beam Charge Negative Hel. Laser On",
			     6000,-10.,50.);
  hQ_BCM_Asym=new TH1F("hQ_BCM_Asym","Beam Asym",
				   1000,-0.05,0.05);
  laserCycleHistos->cd();
  //laser-wise histos
  hL_sum0_laserLeft=new TH1F("hL_sum0_laserLeft",
			     "Acc0 LaserPeriod Sum per Quad- Laser left",
			     1000,sumMin,sumMax);
  hL_diff0_laserLeft=new TH1F("hL_diff0_laserLeft",
			     "Acc0 LaserPeriod Diff per Quad- Laser left",
			     1000,-0.1*diffMax,0.1*diffMax);

  hL_sum0_laserOff=new TH1F("hL_sum0_laserOff",
			    "Acc0 LaserPeriod Sum per Quad- Laser Off",
			    1000,sumMin,sumMax);
  hL_diff0_laserOff=new TH1F("hL_diff0_laserOff",
			     "Acc0 LaserPeriod Diff  per Quad- Laser Off",
			     1000,-0.1*diffMax,0.1*diffMax);


  hL_Asym0_raw_laserLeft=new TH1F("hL_Asym0_raw_laserLeft",
				     "Acc0 LaserPeriod Raw Asym Laser Left",
				   1000,-0.1,0.1);
  hL_Asym0_raw_laserOff=new TH1F("hL_Asym0_raw_laserOff",
				     "Acc0 LaserPeriod Raw Asym Laser Off",
				   1000,-0.1,0.1);
  //
  //Strip Chart style
  mpsHistos->cd();
  hS_Asym0_raw=new TH2F("hS_Asym0_raw","Acc 0 Quad Asym vs MPS",
			1000,0,40000,500,-.1,.1);
  topdir->cd();
  return 0;
}
void fadcAccums::labelSpinSortedHistos(TH1F* histo){
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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int fadcAccums::DefineTree(){
  // data output to tree for each Compton trigger (summed pulses)
  mpsWiseTree=new TTree("mpswise",
			"Accumulator data organized by helicity quartets");
  mpsWiseTree->Branch("Acc0",&accsig[0]);
  mpsWiseTree->Branch("Acc1",&accsig[1]);
  mpsWiseTree->Branch("Acc2",&accsig[2]);
  mpsWiseTree->Branch("Acc3",&accsig[3]);
  mpsWiseTree->Branch("Acc4",&accsig[4]);
  mpsWiseTree->Branch("Acc5",&accsig[5]);
  mpsWiseTree->Branch("Acc6",&accsig[6]);
  mpsWiseTree->Branch("Acc7",&accsig[7]);
  mpsWiseTree->Branch("NAcc0",&naccsig[0]);
  mpsWiseTree->Branch("NAcc1",&naccsig[1]);
  mpsWiseTree->Branch("NAcc2",&naccsig[2]);
  mpsWiseTree->Branch("NAcc3",&naccsig[3]);
  mpsWiseTree->Branch("NAcc4",&naccsig[4]);
  mpsWiseTree->Branch("NAcc5",&naccsig[5]);
  mpsWiseTree->Branch("mpsPedestal",&mpsPedestal);
  mpsWiseTree->Branch("mpsRandomPedestal",&mpsRandomPedestal);
  mpsWiseTree->Branch("mpsTriggerPedestal",&mpsTriggerPedestal);

  //now add on variables from comptonStatus 
  theStatus->DefineStatusBranches(mpsWiseTree);
  //quarter-wise tree
  quartetWiseTree=new TTree("quartetwise",
       "Quartet-summed accumualtor and beam info");
  quartetWiseTree->Branch("PosHelAcc0",&accPosQuad[0]);
  quartetWiseTree->Branch("PosHelAcc1",&accPosQuad[1]);
  quartetWiseTree->Branch("PosHelAcc2",&accPosQuad[2]);
  quartetWiseTree->Branch("PosHelAcc3",&accPosQuad[3]);
  quartetWiseTree->Branch("PosHelAcc4",&accPosQuad[4]);
  quartetWiseTree->Branch("PosHelAcc5",&accPosQuad[5]);
  quartetWiseTree->Branch("PosHelAcc6",&accPosQuad[6]);
  quartetWiseTree->Branch("PosHelAcc7",&accPosQuad[7]);
  //
  quartetWiseTree->Branch("PosHelNSamples0",&naccPosQuad[0]);
  quartetWiseTree->Branch("PosHelNSamples1",&naccPosQuad[1]);
  quartetWiseTree->Branch("PosHelNSamples2",&naccPosQuad[2]);
  quartetWiseTree->Branch("PosHelNSamples3",&naccPosQuad[3]);
  quartetWiseTree->Branch("PosHelNSamples4",&naccPosQuad[4]);
  quartetWiseTree->Branch("PosHelNSamples5",&naccPosQuad[5]);
  quartetWiseTree->Branch("PosHelNSamples6",&naccPosQuad[6]);
  quartetWiseTree->Branch("PosHelNSamples7",&naccPosQuad[7]);
  //

  quartetWiseTree->Branch("NegHelAcc0",&accNegQuad[0]);
  quartetWiseTree->Branch("NegHelAcc1",&accNegQuad[1]);
  quartetWiseTree->Branch("NegHelAcc2",&accNegQuad[2]);
  quartetWiseTree->Branch("NegHelAcc3",&accNegQuad[3]);
  quartetWiseTree->Branch("NegHelAcc4",&accNegQuad[4]);
  quartetWiseTree->Branch("NegHelAcc5",&accNegQuad[5]);
  quartetWiseTree->Branch("NegHelAcc6",&accNegQuad[6]);
  quartetWiseTree->Branch("NegHelAcc7",&accNegQuad[7]);
  //
  quartetWiseTree->Branch("NegHelNSamples0",&naccNegQuad[0]);
  quartetWiseTree->Branch("NegHelNSamples1",&naccNegQuad[1]);
  quartetWiseTree->Branch("NegHelNSamples2",&naccNegQuad[2]);
  quartetWiseTree->Branch("NegHelNSamples3",&naccNegQuad[3]);
  quartetWiseTree->Branch("NegHelNSamples4",&naccNegQuad[4]);
  quartetWiseTree->Branch("NegHelNSamples5",&naccNegQuad[5]);
  quartetWiseTree->Branch("NegHelNSamples6",&naccNegQuad[6]);
  quartetWiseTree->Branch("NegHelNSamples7",&naccNegQuad[7]);
  //

  quartetWiseTree->Branch("quartetHelicityPattern",&quartetHelicity,
			  "quartetHelicityPattern/I");
  quartetWiseTree->Branch("firstMPSnumber",&firstMPS,"firstMPSnumber/I");
  quartetWiseTree->Branch("epics_quartetBCM",&epics_quartetBCM);
  //now add on variables from comptonStatus 
  theStatus->DefineStatusBranches(quartetWiseTree);

  // Add quartet wise bcm info
  quartetWiseTree->Branch("PosHelBCM",&beamPosQuad);
  quartetWiseTree->Branch("NegHelBCM",&beamNegQuad);

  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
int fadcAccums::DoAccums(vmeauxdata* theVMEauxdata,fadcdata *theFADCdata){
  // Should really be called DoMPSstuff
  //
  //  Turn raw Accumulator data into signal (subtract from pedestal, etc.)
  // Do MPS-Wise stuff
  //
  int chan=channel_calorimeter_PMT;
  int laserState=theStatus->GetLaserState();  //sortlaser polarization
  bool laserOn=(laserState==LASER_LEFT || laserState==LASER_RIGHT);
  bool laserOff=(laserState==LASER_LEFTOFF || laserState==LASER_RIGHTOFF);

  float bcm=theStatus->GetCalibratedBCM();
  bool beamOn= (theStatus->GetBeamState()==BEAM_ON);
  bool beamOff= (theStatus->GetBeamState()==BEAM_OFF);

  float bpmsum=theStatus->GetBPMSummed();  //Sum of BPM wires
  
  int sort=theStatus->GetCombinedSpinState();
  int numTriggersPrescaled=theFADCdata->GetSumsNumberTriggers(chan); //number prescaled trigger latches
  int numTriggersAccepted=theFADCdata->GetSumsNumberTriggersSummed(chan);
  int numTriggersScaler=theVMEauxdata->GetTriggerScaler(); //number PMT triggers via IP scaler

  float cavPower=theStatus->GetCalibratedCavityPower();
  hM_SpinState->Fill(sort);    //MPS counts sorted by laser cycle
  hM_BCM->Fill(bcm);     //calibrated bcm
  hM_CavPower->Fill(cavPower);
  hM_BPMSum->Fill(bpmsum);


  if(laserOn){
    hM_CavPower_LaserOn->Fill(cavPower);
  }else if (laserOff){
    hM_CavPower_LaserOff->Fill(cavPower);
  }

  if(beamOn){
    hM_BCM_BeamOn->Fill(bcm);
    hM_BPMSum_BeamOn->Fill(bpmsum);
    hM_BCM_SummedBySpinState->Fill(sort,bcm); //fill weighted by bcm value
    hM_SpinState_BeamOn->Fill(sort);
    hM_Trigs_Scaler->Fill( sort,numTriggersScaler);
    hM_Trigs_Prescaled->Fill( sort,numTriggersPrescaled);
    hM_Trigs_Accepted->Fill(sort,numTriggersAccepted);
    hM_numSums->Fill(numTriggersAccepted);
  }else if(beamOff){
    hM_BCM_BeamOff->Fill(bcm);
    hM_BPMSum_BeamOff->Fill(bpmsum);
  }
  if(theFADCdata->IsAccumValid(chan)){
    for(int accum=0; accum<6; accum++){
      accraw[accum] = theFADCdata->GetAccumValue(chan, accum); //64 bit
      nacc[accum] = theFADCdata->GetAccumNumSamples(chan,accum);		
    }
    // Now handle accumulator combinations (Accums 6 and 7)
    accraw[6] = accraw[1] + accraw[2];
    accraw[7] = accraw[2] + accraw[3];
    nacc[6] = nacc[1] + nacc[2];
    nacc[7] = nacc[2] + nacc[3];
    
    int64_t IntegratedPed;
    double ped=ped_value;  //from textParams class
    // get pedestal from random triggers
    //ped = theFADCdata->GetMPSPedestal();
    mpsPedestal = theFADCdata->GetMPSPedestal();
    mpsRandomPedestal = theFADCdata->GetMPSRandomPedestal();
    mpsTriggerPedestal = theFADCdata->GetMPSTriggerPedestal();
    for(int accum=0; accum<8; accum++){
      IntegratedPed =ped*nacc[accum];
      accsig[accum]=IntegratedPed -accraw[accum];
      naccsig[accum]=nacc[accum];
    }

    float tmp;
    hM_acc0_everything->Fill(accsig[0]);
    //Fill most histograms only if "BeamOn"
    if(beamOn){
      for(int accum=0; accum<8; accum++){
	hM_acc_All[accum]->Fill(accsig[accum]);
	tmp=accsig[accum];   //Pedestal correct accum divided byh bcm
	if(laserState==LASER_RIGHT || laserState==LASER_LEFT){
	  hM_acc_laserOn[accum]->Fill(tmp);
	}else if(laserState==LASER_RIGHTOFF || laserState==LASER_LEFTOFF){
	  hM_acc_laserOff[accum]->Fill(tmp);
	}
      }
    }
    if(beamOn){
      if(laserState==LASER_RIGHT || laserState==LASER_LEFT){
	hM_acc0_beamOff_LaserOn->Fill(accsig[0]);
	hM_acc4_beamOff_LaserOn->Fill(accsig[4]);
      }else if(laserState==LASER_RIGHTOFF || laserState==LASER_LEFTOFF){
	hM_acc0_beamOff_laserOff->Fill(accsig[0]);
	hM_acc4_beamOff_laserOff->Fill(accsig[4]);
      }  
    }
    //if(theStatus->GetBeamState()==BEAM_ON)mpsWiseTree->Fill();
    //DEBUG turn off fr ow
    mpsWiseTree->Fill();
  }
  return 0;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
int fadcAccums::BuildQuad(){
  //called for each MPS (subquad)
  // wrapup==0 for handeling of normal MPS event
  // wrapup==1 to force booking of last laser-period.
  //it builds up the data for one quad
  // valid quads must have stable beam  (all subquads on or all subquads off,etc.)
  // after processing 4th subquad, data is LaserPeriod summation variables
  //  LaserPeriod summation variables are stored when a laser transition is detected
  // by the member function DoLaserTransitin and summation variables are cleared
  //
  if(   !(theStatus->IsHelicityValid() &&theStatus->IsHelicitySynchValid()) ){
    quartetValid=false;
    quartetStable=false;
    return -1;
  }
  int subquad=theStatus->GetIndexQuartet();
  int helicity=theStatus->GetHelicityState();
  int laserState=theStatus->GetLaserState();
  int beamState=theStatus->GetBeamState();
  float bcm=theStatus->GetCalibratedBCM();
  int MPS=theStatus->GetCountMPS();
  /*
  printf("DEBUG BuildQuad  subquad=%8d\n",subquad);
  printf("                helicity=%8d\n",helicity);
  printf("              laserState=%8d\n",laserState);
  printf("               beamState=%8d\n",beamState);
  */
  //printf("DEBUG MPS:index  helicity currentBit %5d %d %d  %d\n",
  //	 MPS,subquad, helicity, theStatus->GetCurrentHelicityBit());
  if(0==subquad){
    //quartetValid and quartetStable used to verify legal quartet pattern
    // and all 4 MPSs are same beam state and laser state
    quartetValid=true;
    quartetStable=true;
    quartetHelicity=helicity;  //set lower bit 
    firstMPS=MPS;
    lastSubQuad=-1;
    laserStateThisQuad=laserState;
    beamStateThisQuad=beamState;
    beamPosQuad=0.;    //start by clearing summing variables if first subquad
    beamNegQuad=0.;
    epics_quartetBCM=0.;
    countPosMinNeg=0;  //count number # helicity windows minus # negative
    for(int i=0; i<NUM_ACCUM_TYPES; i++){
      accPosQuad[i]=0.;   //summing signals
      accNegQuad[i]=0.;
      naccPosQuad[i]=0.;  //summing number of samples in FADC accumuator
      naccNegQuad[i]=0.;
    }
  } else {
    quartetHelicity=quartetHelicity<<4;
    quartetHelicity+=helicity;
  }
  //  if(subquad==3) printf("MPS %8d Helicity pattern %08x \n",firstMPS,quadHelicity);
  if(helicity!=0 && helicity!=1){
    quartetValid=false;
  }
  //
  // check each subquad to assure valid quartet
  if(subquad!=lastSubQuad+1) {
    quartetValid=false;
  }
  if(laserState!=laserStateThisQuad){  //all MPSs must be same laser state
    quartetStable=false;
  }
  if(beamState!=beamStateThisQuad){
    quartetStable=false;
  }
  //
  //Now start summing the subquartets
  // BCM is divided by 4 to get average BCM
  epics_quartetBCM+= 0.25*(theStatus->GetEpicsBCMAverage());
  if(helicity==1){
    beamPosQuad+=bcm;
    countPosMinNeg++;
    for(int i=0; i<NUM_ACCUM_TYPES; i++){
      accPosQuad[i]+=accsig[i];  //sum up the 2 positive helicity subquads
      naccPosQuad[i]+=nacc[i];  //sum up corresponding # samples
    }
  }else{
    beamNegQuad+=bcm;
    countPosMinNeg+=-1;
    for(int i=0; i<NUM_ACCUM_TYPES; i++){
      accNegQuad[i]+=accsig[i];
      naccNegQuad[i]+=nacc[i];
    }
  }
  //
  // If we've summed an entire valid quad, sum it into the laser period stuff
  //
  if(3==subquad){
    if(countPosMinNeg!=0){
      printf("Helicity Quartet Error-Unbalanced helicity sets= %d  MPS= %d\n",
	     countPosMinNeg,MPS);
      quartetValid=false;
    }
  }
  if( !quartetValid){
    printf("Invalid Quartet MPS= %d\n",MPS);
  }
  if(3==subquad && quartetValid ) {
    //valid quads get to here
    if(laserStateThisQuad==laserStateBeingSummed &&beamState==beamStateBeingSummed){
      //quads we want to include in summation get to here
      for(int i=0; i<NUM_ACCUM_TYPES; i++){
	accPosLaserPeriod[i]+=accPosQuad[i];
	accNegLaserPeriod[i]+=accNegQuad[i];
      }
      countQuadsLaserPeriod++;
      beamPosLaserPeriod+=beamPosQuad;
      beamNegLaserPeriod+=beamNegQuad;
    }
  }
  lastSubQuad=subquad;
  //
  //now fill quad-wise histograms
  //
  double diff,sum,background,signal;
  if(3==subquad && quartetValid && quartetStable &&beamState==BEAM_ON) {
    //Accumulator 0 histograms
      if(laserState==LASER_LEFT){
	hQ_aligned0_laserOn->Fill(accPosQuad[0]);
	hQ_anti0_laserOn->Fill(accNegQuad[0]);
      }else if(laserState==LASER_RIGHT){
	hQ_aligned0_laserOn->Fill(accNegQuad[0]);
	hQ_anti0_laserOn->Fill(accPosQuad[0]);
      }else if(laserState==LASER_LEFTOFF){
	hQ_aligned0_laserOff->Fill(accPosQuad[0]);
	hQ_anti0_laserOff->Fill(accNegQuad[0]);
      }else if(laserState==LASER_RIGHTOFF){
	hQ_aligned0_laserOff->Fill(accNegQuad[0]);
	hQ_anti0_laserOff->Fill(accPosQuad[0]);
      }
      diff=accPosQuad[0]-accNegQuad[0];
      sum=accPosQuad[0]+accNegQuad[0];
      background=accSumLastLaserOff[0]*(beamPosQuad+beamNegQuad);
      signal=sum-background;
      if(laserState==LASER_LEFT){
	//Laser On Left Polarition Histograms
	hQ_sum0_laserLeft->Fill(sum);
	hQ_diff0_laserLeft->Fill(diff);
	if(sum>1.0E-5) 
	  hQ_Asym0_raw_laserLeft->Fill(diff/sum);
	if(prevLaserStateEncountered==1){
	  if(abs(signal)>1.e-5){
	    hQ_Signal0->Fill(signal);
	    hQ_Asym0->Fill(diff/signal);
	  }
	}
      } else {
	hQ_sum0_laserOff->Fill(sum);
	hQ_diff0_laserOff->Fill(diff);
	if(sum>1.0E-5) {
	  hQ_Asym0_raw_laserOff->Fill(diff/sum);
	  hS_Asym0_raw->Fill(firstMPS,diff/sum);
	}
      }
    //Accumulator 4 histograms
      diff=accPosQuad[4]-accNegQuad[4];
      sum=accPosQuad[4]+accNegQuad[4];
      background=accSumLastLaserOff[4]*(beamPosQuad+beamNegQuad);
      signal=sum-background;
      if(laserState==LASER_LEFT){
	hQ_sum4_laserLeft->Fill(sum);
	hQ_diff4_laserLeft->Fill(diff);
	if(sum>1.0E-5) 
	  hQ_Asym4_raw_laserLeft->Fill(diff/sum);
	if(prevLaserStateEncountered==1){
	  if(abs(signal)>1.e-5){
	    hQ_Signal4->Fill(signal);
	    hQ_Asym4->Fill(diff/signal);
	  }
	}
      } else {
	//Laser Off Histograms
	hQ_sum4_laserOff->Fill(sum);
	hQ_diff4_laserOff->Fill(diff);
	if(sum>1.0E-5) {
	  hQ_Asym4_raw_laserOff->Fill(diff/sum);

      }
    }
  }
  //Beam charge histos
  if(3==subquad && quartetValid && quartetStable && beamState==BEAM_ON){
    if(laserState==LASER_RIGHT || laserState==LASER_LEFT) {
      diff=beamPosQuad-beamNegQuad;
      sum=beamPosQuad+beamNegQuad;;
      if(sum>1.0E-5) hQ_BCM_Asym->Fill(diff/sum);
      hQ_BCM_laserOn_P->Fill(beamPosQuad);
      hQ_BCM_laserOn_N->Fill(beamNegQuad);
    }else if(laserState==LASER_RIGHTOFF || laserState==LASER_LEFTOFF){
      hQ_BCM_laserOff_P->Fill(beamPosQuad);
      hQ_BCM_laserOff_N->Fill(beamNegQuad);
    }
  }

  //Beam Off Accumulator Histograms
  if(3==subquad && quartetValid &&quartetStable  && beamState==BEAM_OFF) {
    //Accumulator 0 histograms
    if(laserState==LASER_LEFT || laserState==LASER_LEFTOFF){
      hQ_aligned0_beamOff->Fill(accPosQuad[0]);
      hQ_anti0_beamOff->Fill(accNegQuad[0]);
    }else if(laserState==LASER_RIGHT  || laserState==LASER_RIGHTOFF){
      hQ_aligned0_beamOff->Fill(accNegQuad[0]);
      hQ_anti0_beamOff->Fill(accPosQuad[0]);
    }
    diff=accPosQuad[0]-accNegQuad[0];
    sum=accPosQuad[0]+accNegQuad[0];
    hQ_sum0_beamOff->Fill(sum);
    hQ_diff0_beamOff->Fill(diff);
    if(sum>1.0E-5)  hQ_Asym0_raw_beamOff->Fill(diff/sum);
    diff=accPosQuad[4]-accNegQuad[4];
    sum=accPosQuad[4]+accNegQuad[4];
    hQ_sum4_beamOff->Fill(sum);
    hQ_diff4_beamOff->Fill(diff);
  }
  if(3==subquad && quartetValid && quartetStable) {
    quartetWiseTree->Fill();
  }
  return 0;
}
  //
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
int fadcAccums::DoLaserTransition(int wrapup){
  // called when a laser state transition is detected
  // AND at end of file (called with wrapup==1)
  // it is called BEFORE the new MPS data is processed, so this routine
  // works with the previous laser cycles
  //  If the new state is a laser on state, the laser three laser states
  // should have been laser-off, laser-on, laser-off
  int laserState=theStatus->GetLaserState();
  if(wrapup==0){
    if(laserState==laserStateBeingSummed||laserState==LASER_UNKNOWN){
      prevLaserStateEncountered=laserState;
      return 0;   //just keep doing what we've been doing
    }
  }
  //encountered a new laserState that we want to sum.  Store what we've got
  //from previous laser state
  if(countQuadsLaserPeriod!=0){
    //      StoreAccumulatorInfo(laserStateThisQuad, beamStateThisQuad,
    //		       countQuadsLaserPeriod,
    //		       accDiffLaserPeriod,accSumLaserPeriod);
   }  
  //keep laser-off data for background subtraction in next laser on period

  if(beamStateBeingSummed==BEAM_ON) {
    if(laserStateBeingSummed==LASER_RIGHTOFF ||laserStateBeingSummed==LASER_LEFTOFF){
      double beamBCM=beamPosLaserPeriod+beamNegLaserPeriod;
      if(beamBCM>0.0){
	for(int i=0; i<NUM_ACCUM_TYPES; i++){
	  accSumLastLaserOff[i]=
	    (accPosLaserPeriod[i]+accNegLaserPeriod[i])/beamBCM;
	}
	prevLaserStateEncountered=true;
	bcmLastLaserOff=beamPosLaserPeriod+beamNegLaserPeriod;
	//	printf("Laser Off backgroundSaved %e\n",accSumLastLaserOff[0]);
      }else{
	if(countQuadsLaserPeriod>0){
	  printf("Warning: bcm=0 for laser period taged as BEAM ON \n");
	}
      }
    }
  }
  //do laser-wise bookkeeping
  if(countQuadsLaserPeriod!=0){
    double diff,sum;
    sum= (accPosLaserPeriod[0]+accNegLaserPeriod[0])/countQuadsLaserPeriod;
    diff=(accPosLaserPeriod[0]-accNegLaserPeriod[0])/countQuadsLaserPeriod;;
    //    printf(" laser period Acc0 diff %e  sum %e \n", diff, sum);
    if(laserStateBeingSummed==LASER_LEFT){
      hL_sum0_laserLeft->Fill(sum);
      hL_diff0_laserLeft->Fill(diff);
      hL_Asym0_raw_laserLeft->Fill(diff/sum);
    }else if(laserStateBeingSummed==LASER_LEFTOFF||laserState==LASER_RIGHTOFF){
      hL_sum0_laserOff->Fill(sum);
      hL_diff0_laserOff->Fill(diff);
      hL_Asym0_raw_laserOff->Fill(diff/sum);
    }    
  }
  // prepare summing variables for next cycle
  for(int i=0; i<NUM_ACCUM_TYPES; i++){
    accPosLaserPeriod[i]=0.;    
    accNegLaserPeriod[i]=0.;
  }
  beamPosLaserPeriod=0.;
  beamNegLaserPeriod=0.;
  countQuadsLaserPeriod=0;
  firstMPS=theStatus->GetCountMPS();  //record start of new laser-off
  laserStateBeingSummed=laserState;   //start summing these states in BuildQuad
  beamStateBeingSummed=theStatus->GetBeamState();  //don't sum if  beam status to change
  return 0;
} 
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
bool fadcAccums::GetSummedAccumulatorRecord(int record,int accumNumber,
					    int* laserState,  int* beamState,
					    int* countQuads,
					    double* diff, double* sum){
  //retrieves accumulator info for 1 type of accumulator summed over 1 quad
  // record =0 for current quad,  record=1 prevous quad, etc.
  // accumNumber=- for accumulator zero, etc.
  if(record>=0 &&record<NUM_HISTORY_LASER_PERIODS){
    int pointer=pointerHistory-record;
    if (pointer<0) pointer+=NUM_HISTORY_LASER_PERIODS;
    *laserState=laserStateHistory[pointer];
    *countQuads=countQuadsHistory[pointer];
    *beamState=beamStateHistory[pointer];
    *diff=accDiffHistory[accumNumber][pointer];
    *sum=accSumHistory[accumNumber][pointer];
    return true;
  }else{
    *laserState=LASER_UNKNOWN;
    *countQuads=0;
    *diff=0.;
    *sum=1.;
  return false;
  }
}
void fadcAccums::StoreAccumulatorInfo(
				      int laserState,  int beamState, int countQuads,
				       double diff[NUM_ACCUM_TYPES],
				       double sum[NUM_ACCUM_TYPES]){
   //stores accumulator info for all accumulator types summed over 1 quad
   pointerHistory++;   //bump circular pointer
   if(pointerHistory>=NUM_HISTORY_LASER_PERIODS) pointerHistory=0;
   int pointer=pointerHistory;
   laserStateHistory[pointer]=laserState;
   countQuadsHistory[pointer]=countQuads;
   beamStateHistory[pointer]=beamState;
   for(int i=0; i<NUM_ACCUM_TYPES; i++){
     accDiffHistory[i][pointer]=diff[i];
     accSumHistory[i][pointer]=sum[i];
   }
   return;
}
