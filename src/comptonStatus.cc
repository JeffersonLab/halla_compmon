//*****************************
//*  comptonStatus.cc
//*  created 02/11/2015 G.B. Franklin
//*  Designed for keeping track of Compton status (helicity, laser on etc.)
//*
//*****************************

#define MAX_ROOTFILE_SIZE 10000000000  // Roughly 10 GB

#include "comptonStatus.h"

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
comptonStatus::comptonStatus(){
  TTree::SetMaxTreeSize(MAX_ROOTFILE_SIZE);
  //constructor
  runWiseTree=0;
  runWiseTreeFilled=false;  //used to allow single fill of runWiseTree
  return;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
comptonStatus::comptonStatus(textParams* theParamsIn){
  TTree::SetMaxTreeSize(MAX_ROOTFILE_SIZE);
  //constructor
  theParams=theParamsIn;
  theHelicityTracker=new helicityTracker();
  return;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TString comptonStatus::DecodeLaserState(int laserState){
  if(laserState==LASER_RIGHT){
    return TString("Right");
  }else if (laserState==LASER_LEFT){
    return TString("Left");
  }else if (laserState==LASER_RIGHTOFF){
    return TString("RightOff");
  }else if (laserState==LASER_LEFTOFF){
    return TString("LeftOff");
  }else if (laserState==LASER_UNKNOWN){
    return TString("Unknown");
  }
  return TString("Undefined");
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TString comptonStatus::DecodeBeamState(int beamState){
  if(beamState==BEAM_ON){
    return TString("On");
  }else if (beamState==BEAM_OFF){
    return TString("Off");
  }else if (beamState==BEAM_UNKNOWN){
    return TString("Unknown");
  } else {
    return TString("Undefined");
  }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int comptonStatus::DefineTrees(){
  printf("DEBUG *****  runWiseTree created\n");

  runWiseTree=new TTree("runwise",
			"Info for entire run");
  runWiseTree->Branch("runNumber", &run_number, "runNumber/I");
  //Stuff from compmon parameters
  runWiseTree->Branch("FADC_ped_value",&ped_value);
  runWiseTree->Branch("cavityPowerOnMin",&cavityPowerOnMin);
  runWiseTree->Branch("cavityPowerOffMax",&cavityPowerOffMax);
  runWiseTree->Branch("BPMSumOnMin",&BPMSum_OnMin);
  runWiseTree->Branch("BPMSumOffMax",&BPMSum_OffMax);
  runWiseTree->Branch("BCMOnMin",&BCM_OnMin);
  runWiseTree->Branch("BCMOffMax",&BCM_OffMax);
  runWiseTree->Branch("useBPMSumCuts",&useBPMSumCuts,"useBPMSumCuts/I");
  
  //Stuff from CODA Event
  runWiseTree->Branch("FADC_ithrnear", &ithr_near, "FADC_ithrnear/I");
  runWiseTree->Branch("FADC_ithrfar", &ithr_far, "FADC_ithrfar/I");
  runWiseTree->Branch("FADC_n4before", &n4before, "FADC_n4before/I");
  runWiseTree->Branch("FADC_n4after", &n4after, "FADC_n4after/I");
  runWiseTree->Branch("FADC_n5before", &n5before, "FADC_n5before/I");
  runWiseTree->Branch("FADC_n5after", &n5after, "FADC_n5after/I");
  runWiseTree->Branch("FADC_dac", &dac, "FADC_dac/I");
  runWiseTree->Branch("HTB_rampdelay", &rampdelay, "HTB_rampdelay/I");
  runWiseTree->Branch("HTB_inttime", &inttime, "HTB_inttime/I");

  // Stuff from EPICS events
  epicsWiseTree = new TTree("epicswise","Info from just EPICS events");
  DefineEpicsBranches(epicsWiseTree);

  return 0;  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int comptonStatus::DefineEpicsBranches(TTree* mytree){
  mytree->Branch("countEpics", &countEpics, "countEpics/I");
  //SLOW epics stuff
  mytree->Branch("epics_evnum",&epics_evnum,"epics_evnum/I");
  mytree->Branch("epics_PMTRate",&epics_PMTRate,"epics_PMTRATE/F");
  mytree->Branch("epics_PMTRateHigh",&epics_PMTRateHigh,"epics_PMTRATEHigh/F");
  mytree->Branch("epics_bcm_average",&epics_hacbmf,"epics_bcm_average/F");
  mytree->Branch("epics_crystalHV", &epics_crystalHV, "epics_crystalHV/F");
  mytree->Branch("epics_crystal_current", &epics_crystalCurrent, "epics_crystal_current/F");
  mytree->Branch("epics_horizontalFingerHV", &epics_horizontalFingerHV, "epics_horizontalFingerHV/F");
  mytree->Branch("epics_horizontalFingerCurrent", &epics_horizontalFingerCurrent, "epics_horizontalFingerCurrent/F");
  mytree->Branch("epics_verticalFingerHV", &epics_verticalFingerHV, "epics_verticalFingerHV/F");
  mytree->Branch("epics_verticalFingerCurrent", &epics_verticalFingerCurrent, "epics_verticalFingerCurrent/F");
  mytree->Branch("epics_tablePosX", &epics_tablePosX, "epics_tablePosX/F");
  mytree->Branch("epics_tablePosY", &epics_tablePosY, "epics_tablePosY/F");
  mytree->Branch("epics_ihwp_in", &epics_ihwp_in, "epics_ihwp_in/I");
  mytree->Branch("epics_wein_right", &epics_wein_right, "epics_wein_right/I");
  mytree->Branch("epics_bpmAx",&epics_aPosX,"epics_bpmAx/F"); //BPM info
  mytree->Branch("epics_bpmAy",&epics_aPosY,"epics_bpmAy/F");
  mytree->Branch("epics_bpmBx",&epics_bPosX,"epics_bpmBx/F");
  mytree->Branch("epics_bpmBy",&epics_bPosY,"epics_bpmBy/F");
  mytree->Branch("epics_coda_deadtime", &epics_coda_deadtime, "epics_coda_deadtime/F");
  mytree->Branch("epics_Thermo1",&epics_Thermo1,"epics_Thermo1/F");
  mytree->Branch("epics_Thermo2",&epics_Thermo2,"epics_Thermo2/F");
  mytree->Branch("epics_TimeStamp",&epics_TimeStamp,"epics_TimeStamp/F");

  return 0;
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int comptonStatus::DefineStatusBranches(TTree* mytree){
// External-data branches containing status info shared by all tree
  mytree->Branch("runNumber", &run_number, "runNumber/I");
  mytree->Branch("helicityState", &helicityState, "helicityState/I");
  mytree->Branch("laserState", &currentLaserState, "laserstate/I");
  mytree->Branch("combinedSpinState", &currentSpinState, "CombinedSpinState/I");
  mytree->Branch("beamState", &beamState, "beamState/I");
  mytree->Branch("numTriggers", &numTriggers, "numTriggers/I");
  mytree->Branch("mpsCoda", &mpsCoda, "mpsCoda/I"); //mps # read from Coda header
  mytree->Branch("mpsSignal", &mpsSignal, "mpsSignal/I"); //mps signal state
  mytree->Branch("mpsScaler", &mpsScaler, "mpsScaler/I");  //scaler mps read
  mytree->Branch("mpsAnalyzed", &countMPS, "mpsAnalyzed/I"); //mps counter
  mytree->Branch("mpsSinceEpics",&countMPSsinceEPICS,"mpsSinceEpics/I");
  mytree->Branch("countLaserCycles", &countLaserCycles, "countLaserCycles/i");

  mytree->Branch("cavPowerCalibrated", &cavPowerCalibrated, "cavPowerCalibrated/F");
  mytree->Branch("rawCavPower", &rawCavPowFloat, "rawCavPower/F");
  mytree->Branch("s1power", &ip_s1power, "s1power/F");
  mytree->Branch("s2power", &ip_s2power, "s2power/F");

  mytree->Branch("bcm", &calbcm, "bcm/F");
  mytree->Branch("rawBCM", &rawBCMFloat, "rawBCM/F");
  mytree->Branch("bpmSum", &bpmsum, "bpmSum/F");

  mytree->Branch("bpmAx_raw", &ip_bpmAx, "bpmAx_raw/F");
  mytree->Branch("bpmAy_raw", &ip_bpmAy, "bpmAy_raw/F");
  mytree->Branch("bpmBx_raw", &ip_bpmBx, "bpmBx_raw/F");
  mytree->Branch("bpmBy_raw", &ip_bpmBy, "bpmBy_raw/F");

  // The lab coordinates
  mytree->Branch("bpmAx", &bpmAx, "bpmAx/F");
  mytree->Branch("bpmAy", &bpmAy, "bpmAy/F");
  mytree->Branch("bpmBx", &bpmBx, "bpmBx/F");
  mytree->Branch("bpmBy", &bpmBy, "bpmBy/F");

  mytree->Branch("clockRun", &clockscaler, "clockRun/i");
  mytree->Branch("clockIP", &clockIP, "clockIP/i");
  mytree->Branch("clockdiff", &clock_diff, "clockdiff/I");

  //SLOW epics stuff
  DefineEpicsBranches(mytree);

  //mytree->Branch("rtcavpow", &rtcavpow, "rtcavpow/F");
  mytree->Branch("epbcmu3", &epbcmu3, "epbcmu3/F");
  mytree->Branch("buflen", &buflen, "buflen/I");
  mytree->Branch("l1a", &l1a, "l1a/i");
  mytree->Branch("dithering", &dithering, "dithering/I");

  mytree->Branch("beameng", &epbeameng, "beameng/F");
  mytree->Branch("transmit", &eptransmit, "transmit/F");
  mytree->Branch("beam_trip", &beam_trip, "beam_trip/I");
  mytree->Branch("beam_burp", &beam_burp, "beam_burp/I");
  mytree->Branch("eppol_burp", &eppol_burp, "eppol_burp/I");
  mytree->Branch("HV_trip", &HV_trip, "HV_trip/I");
  mytree->Branch("rate_fluct", &rate_cut, "rate_fluct/I");
  mytree->Branch("epics_dead", &epics_dead, "epics_dead/I");
  mytree->Branch("test0", &test0, "test0/I");
  mytree->Branch("test1", &test1, "test1/I");
  mytree->Branch("test2", &test2, "test2/I");
  mytree->Branch("test3", &test3, "test3/I");
  mytree->Branch("test4", &test4, "test4/I");
  mytree->Branch("test5", &test5, "test5/I");
  mytree->Branch("test6", &test6, "test6/I");
  mytree->Branch("test7", &test7, "test7/I");
  mytree->Branch("test8", &test8, "test8/I");
  mytree->Branch("test9", &test9, "test9/I");
  mytree->Branch("test10", &test10, "test10/I");
  mytree->Branch("test11", &test11, "test11/I");
  mytree->Branch("test12", &test12, "test12/I");
  mytree->Branch("test13", &test13, "test13/I");
  mytree->Branch("test14", &test14, "test14/I");
  mytree->Branch("test15", &test15, "test15/I");
  return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
int comptonStatus::newRun(){
  //Init all counters and status
  countMPS=0;
  countEpics=0;
  countMPSsinceEPICS=-1;  //not EPICS event yet encountered
  epics_evnum=-1;  //receives EPICS event number (not valid till first EPICS event)
  SetEpicsDefaults();  //erase old EPICS values
  //helicity pattern decoding info
  //look up run-dependent helicity info
  run_number=theParams->getRunNumber();
  helicityDelay=theParams->getInt("helicity_delay");
  helicityShiftBits=theParams->getInt("helicity_bits");
  helicityStructure=theParams->getInt("helicity_structure");
  helicityPredictorBitFlip=theParams->getInt("helicity_predictor_bit_flip");
  ped_value=theParams->getFloat("ped_value");
  //
  cavPowerCalibration=theParams->getFloat("cavity_calibration");
  cavPowerPedestal=theParams->getFloat("cavity_power_pedestal");
  BCMCalibration=theParams->getFloat("BCM_calibration");
  BCMPedestal=theParams->getFloat("BCM_pedestal");
  BPMSumCalibration=theParams->getFloat("BPMSum_calibration");
  BPMSumPedestal=theParams->getFloat("BPMSum_pedestal");
  //
  BPMSum_OnMin=theParams->getFloat("beam_on_min_BPMSum");
  BPMSum_OffMax=theParams->getFloat("beam_off_max_BPMSum");
  BCM_OnMin=theParams->getFloat("beam_on_min");
  BCM_OffMax=theParams->getFloat("beam_off_max");
  useBPMSumCuts=theParams->getInt("use_BPMSumCuts");
  cavityPowerOnMin=theParams->getFloat("cavity_power_on_min");
  cavityPowerOffMax=theParams->getFloat("cavity_power_off_max");
  clockRateIP=theParams->getFloat("clockRateIP");  
  // BPM parameters
  BPM2Axm_pedestal=theParams->getFloat("BPM2Axm_pedestal");
  BPM2Axp_pedestal=theParams->getFloat("BPM2Axp_pedestal");
  BPM2Aym_pedestal=theParams->getFloat("BPM2Aym_pedestal");
  BPM2Ayp_pedestal=theParams->getFloat("BPM2Ayp_pedestal");
  BPM2Bxm_pedestal=theParams->getFloat("BPM2Bxm_pedestal");
  BPM2Bxp_pedestal=theParams->getFloat("BPM2Bxp_pedestal");
  BPM2Bym_pedestal=theParams->getFloat("BPM2Bym_pedestal");
  BPM2Byp_pedestal=theParams->getFloat("BPM2Byp_pedestal");
  BPM2A_alphax=theParams->getFloat("BPM2A_alphax");
  BPM2A_alphay=theParams->getFloat("BPM2A_alphay");
  BPM2B_alphax=theParams->getFloat("BPM2B_alphax");
  BPM2B_alphay=theParams->getFloat("BPM2B_alphay");
  BPM2A_xoff=theParams->getFloat("BPM2A_xoff");
  BPM2A_yoff=theParams->getFloat("BPM2A_yoff");
  BPM2B_xoff=theParams->getFloat("BPM2B_xoff");
  BPM2B_yoff=theParams->getFloat("BPM2B_yoff");
  BPM2A_sensitivity=theParams->getFloat("BPM2A_sensitivity");
  BPM2B_sensitivity=theParams->getFloat("BPM2B_sensitivity");
  BPM2A_angle=theParams->getFloat("BPM2A_angle")*0.0174533;
  BPM2B_angle=theParams->getFloat("BPM2B_angle")*0.0174533;
  //
  countLaserCycles=0;
  subcountMPSLaserCycles=0; //number MPS within current laser mode
  laserStateValid=false;
  currentLaserState=LASER_UNKNOWN;  //enum
  for(int i=0; i<6; i++){
    laserStateFound[i]=false;
  }
  beamStatePrev=BEAM_UNKNOWN;
  theHelicityTracker->SetHelicityMode(helicityDelay,helicityShiftBits,
				     helicityStructure,helicityPredictorBitFlip);
  theHelicityTracker->newRun();

  printf("Initializing Run                 %8d\n",run_number);
  printf("   Helicity Delay                %8d\n", helicityDelay);
  printf("   Helicity Shift Register Bits  %8d\n", helicityShiftBits);
  printf("   Helicity Structure            %8d\n",helicityStructure);
  return 0;
}
bool comptonStatus::newMPS(int codaEventNumber, fadcdata* theFADCdata, vmeauxdata* theAuxData){
  countMPS++;   //mps count via count of analyzed Event 1s
  if(countMPSsinceEPICS>=0) countMPSsinceEPICS++;
  mpsCoda=codaEventNumber;  //mps count via CODA Event 1 header
  
  mpsScaler=theAuxData->GetMPSScaler();  //mps count via VME scaler
  mpsSignal=theAuxData->GetMPSSignal();
  numTriggers=theFADCdata->GetSumsNumberTriggers(0); //number of calorimeter triggers
  dithering=theAuxData->GetDithering();

  rtcavpol= theAuxData->GetCavityPolBit();  
  currentHelicityBit=theAuxData->GetHelicityBit(); //from input register 
  inttime=theAuxData->GetIntTime();
  ip_s1power=theAuxData->GetVtoFPowLeftIPScaler();
  ip_s2power=theAuxData->GetVtoFPowRightIPScaler();

  //StuffFADC settings that should not actually be changing every event...
  ithr_near=theFADCdata->GetThresh(0,0); //Assume FADC channel 0 for now
  ithr_far=theFADCdata->GetThresh(0,1);
  n4before=theFADCdata->GetN4before(0);
  n4after=theFADCdata->GetN4after(0);
  n5before=theFADCdata->GetN5before(0);
  n5after=theFADCdata->GetN5after(0);
  dac=theFADCdata->GetDac(0); //FADC DAC  (not MiniMegan Pulse DAC)

  //update Helicity Tracker using autosynch mode since we don't have the
  // helicity set synch bit (use helicitySetSynchBit=-1
  //  last argument is the verbose control (Use 1 for lasertranistion logs and errors only)
  //1 after mpsScaler sets verbose level, 0 is off/quiet, 2 is a lot
  helicityState=
    theHelicityTracker->newHelicityWindow(currentHelicityBit,-1,mpsScaler,1);
  statusHW=theHelicityTracker->GetHelicityStatus(); //pointer to helicity info
  if(countMPS%100==0){
    if(!statusHW->helicitySynchValid){
       printf("Unable to synch to helicity pattern currentHelicityBit= %d\n",
	      currentHelicityBit);
    }
  }
  //  if( !statusHW->helicityBitInAgreement &&statusHW->helicityPredictionValid) {
  //   printf("DATA Error: Current helicity bit / prediction %d / %d \n",
  //	   statusHW->currentHelicityBit,statusHW->predictedCurrentBit);
  //}
  //laser cavity power
  old_clockscaler=clockscaler;
  clockscaler = theAuxData->GetClockScaler();
  clock_diff=clockscaler-old_clockscaler;
  //  int clockRate=4.0E7;
  clockIP=theAuxData->GetGatedClock();  //use if we have an IP scale clock
  bcmscaler =theAuxData->GetGatedBCM();  //bcm value from VTF gated scaler
  rawBCMFloat= bcmscaler;

  float freqConversion = float(clockRateIP)/float(clockIP);

  //backup info for beam on/off
  int bpmXP, bpmXM;
  int bpmYP, bpmYM;
  bpmYM=theAuxData->GetVtoFBPM2AymIPScaler();
  bpmYP=theAuxData->GetVtoFBPM2AypIPScaler();
  bpmXM=theAuxData->GetVtoFBPM2AxmIPScaler();
  bpmXP=theAuxData->GetVtoFBPM2AxpIPScaler();
  // quick BPM sum
  bpmsum=bpmXP+bpmXM+bpmYP+bpmYM;
  // Determine the rotated BPM values
  float rot_bpmAx=ComputeBPMPosition(bpmXP,bpmXM,BPM2Axp_pedestal,
      BPM2Axm_pedestal,BPM2A_alphax, BPM2A_sensitivity,freqConversion);
  float rot_bpmAy=ComputeBPMPosition(bpmYP,bpmYM,BPM2Ayp_pedestal,
      BPM2Aym_pedestal,BPM2A_alphay, BPM2A_sensitivity,freqConversion);
  // For backwards compatibility, also do the "raw" values, with no
  // calibrations
  if(bpmYM+bpmYP>0){
    ip_bpmAy= (bpmYP-bpmYM)/float(bpmYP+bpmYM);
  }
  if(bpmYM+bpmYP>0){
    ip_bpmAx= (bpmXP-bpmXM)/float(bpmXP+bpmXM);
  }
  // Now do the second BPM
  bpmYM=theAuxData-> GetVtoFBPM2BymIPScaler();
  bpmYP=theAuxData-> GetVtoFBPM2BypIPScaler();
  bpmXM=theAuxData-> GetVtoFBPM2BxmIPScaler();
  bpmXP=theAuxData-> GetVtoFBPM2BxpIPScaler();
  float rot_bpmBx=ComputeBPMPosition(bpmXP,bpmXM,BPM2Bxp_pedestal,
      BPM2Bxm_pedestal,BPM2B_alphax, BPM2B_sensitivity,freqConversion);
  float rot_bpmBy=ComputeBPMPosition(bpmYP,bpmYM,BPM2Byp_pedestal,
      BPM2Bym_pedestal,BPM2B_alphay, BPM2B_sensitivity,freqConversion);
  // For backwards compatibility, also do the "raw" values, with no
  // calibrations
  if(bpmYM+bpmYP>0){
    ip_bpmBy= (bpmYP-bpmYM)/float(bpmYP+bpmYM);
  }
  if(bpmXM+bpmXP>0){
    ip_bpmBx= (bpmXP-bpmXM)/float(bpmXP+bpmXM);
  }
  // Now determine the lab frame beam positions
  float BPM2A_sintheta = sin(BPM2A_angle);
  float BPM2A_costheta = cos(BPM2A_angle);
  float BPM2B_sintheta = sin(BPM2B_angle);
  float BPM2B_costheta = cos(BPM2B_angle);
  ComputeBPMPositionLab(rot_bpmAx,rot_bpmAy,BPM2A_sintheta,BPM2A_costheta,
      BPM2A_xoff,BPM2A_yoff,bpmAx,bpmAy);
  ComputeBPMPositionLab(rot_bpmBx,rot_bpmBy,BPM2B_sintheta,BPM2B_costheta,
      BPM2B_xoff,BPM2B_yoff,bpmBx,bpmBy);

  //assume 40 MHz clock, but this may not be correct for Spring 2016
  rawCavPowFloat=theAuxData->GetCavityPowerScaler();
  if(clockIP<0){
    cavPowerCalibrated=0;
    calbcm=0;
    bpmsum=0;
  }else{
    cavPowerCalibrated=clockRateIP*( cavPowerCalibration*
    theAuxData->GetCavityPowerScaler() )/clockIP;
    cavPowerCalibrated+=-cavPowerPedestal;
    calbcm = (float) bcmscaler/ (float)clockIP;
    calbcm *= clockRateIP*BCMCalibration;
    calbcm +=-BCMPedestal;
    bpmsum *=clockRateIP*BPMSumCalibration/(float)clockIP;
    bpmsum += -BPMSumPedestal;
  }
  beamStatePrev=beamState;      //keep previous MPS beam on/off status
  if(useBPMSumCuts==0){     //normal BCM beam state cuts
    if(calbcm > BCM_OnMin){
      beamState=BEAM_ON;
    }else if(calbcm<BCM_OffMax){
      beamState=BEAM_OFF;
    }else{
      beamState=BEAM_UNKNOWN;
    }
  }else{
    if(bpmsum > BPMSum_OnMin){   //BPM sums for runs with bad BCM scaler
      beamState=BEAM_ON;
    }else if(bpmsum<BPMSum_OffMax){
      beamState=BEAM_OFF;
    }else{
      beamState=BEAM_UNKNOWN;
    }
  }
  bool transition=SetLaserState(); //determine laser state
  if(transition) countLaserCycles++;  //count change in laser state
  if(currentLaserState==LASER_UNKNOWN || helicityState==2){
    currentSpinState=SPIN_UNKNOWN;
  } else  if(currentLaserState==LASER_LEFT){
    if(helicityState==0){
      currentSpinState=SPIN_LN;
    }else if(helicityState==1){
      currentSpinState=SPIN_LP;
    }
  }else if(currentLaserState==LASER_RIGHT){
    if(helicityState==0){
      currentSpinState=SPIN_RN;
    }else if (helicityState==1){
      currentSpinState=SPIN_RP;
    }
  }else if (currentLaserState==LASER_LEFTOFF){
    if(helicityState==0){
      currentSpinState=SPIN_LNOFF;
    }else if (helicityState==1){
      currentSpinState=SPIN_LPOFF;
    }
  }else if (currentLaserState==LASER_RIGHTOFF){	
    if(helicityState==0){
      currentSpinState=SPIN_RNOFF;
    }else if (helicityState==1){
      currentSpinState=SPIN_RPOFF;
    }
  }else{
    currentSpinState=SPIN_UNKNOWN;
  }
  //transfer in scalers.
  test0 = theAuxData->GetIPScaler(0);
  test1 = theAuxData->GetIPScaler(1);
  test2 = theAuxData->GetIPScaler(2);
  test3 = theAuxData->GetIPScaler(3);
  test4 = theAuxData->GetIPScaler(4);
  test5 = theAuxData->GetIPScaler(5);
  test6 = theAuxData->GetIPScaler(6);
  test7 = theAuxData->GetIPScaler(7);
  test8 = theAuxData->GetIPScaler(8);
  test9 = theAuxData->GetIPScaler(9);
  test10= theAuxData->GetIPScaler(10);
  test11= theAuxData->GetIPScaler(11);
  test12= theAuxData->GetIPScaler(12);
  test13= theAuxData->GetIPScaler(13);
  test14= theAuxData->GetIPScaler(14);
  test15= theAuxData->GetIPScaler(15);
  if( !runWiseTreeFilled){
    runWiseTree->Fill();
    runWiseTreeFilled=true;
  }
  return transition;
}
int comptonStatus::GetHelicityState(){
  return statusHW->helicityState;
 }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

void comptonStatus::DebugDump(int mode){
  // Dump Status info for debugging
  //  mode==0   lots of info
  // mode==1  single line for each MPS
  // mode==2 output only on laser or beam on/off ransition
  if(mode==0){
    printf("comptonStatus Dump \n");
    
    printf("mpsCoda  (from Coda Event Header) %d\n",mpsCoda);
    printf("mpsScaler)from VME Scalar       ) %d\n",mpsScaler);
    
    printf("Status bits \n");
    printf("   mpsSignal (mps bit) %d\n",mpsSignal);
    printf("   dithering           %d\n",dithering);
    printf("   rtcavpol             %f\n",rtcavpol);
    //    printf("   cavityPowerBit       %d\n",cavityPowerBit);
    printf("   beamState               %d\n",beamState);
    
    printf("integers: \n");
    printf("   bcmscaler                     %d\n",bcmscaler);
    printf("   calbcm (calibrated bcmscaler) %f\n",calbcm);
    printf("   rtcavPower                    %d\n",rtcavpow);
    //    printf("   cavityPowerScaler             %d\n",theAuxData->GetCavityPowerScaler());
    printf("   clock_diff                    %d\n",clock_diff);
    printf("   CavCalibration                %f\n",cavPowerCalibration);
  }else if(mode==1){
    if(statusHW->indexQuartet==0) printf("\n");
    printf("mps %5d LaserState %2d SpinState %2d index %2d helicity %d\n",
	   countMPS,currentLaserState,currentSpinState,
	   statusHW->indexQuartet,statusHW->helicityState);
  }else if(mode==2){
    if(beamState!=beamStatePrev || 
       (currentLaserState!=previousLaserState)){
      printf("%90d  BeamState %8s LaserState %8s \n",
	     countMPS, DecodeBeamState(beamState).Data(),
	     DecodeLaserState(currentLaserState).Data());
    }
  }else if(mode==3){
    if(beamState!=beamStatePrev || 
       (currentLaserState!=lastKnownLaserState&&currentLaserState!=LASER_UNKNOWN)){
    printf("%90d  BeamState %8s LaserState %8s \n",
	     countMPS, DecodeBeamState(beamState).Data(),
	     DecodeLaserState(currentLaserState).Data());
    }
  }
  return ;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Use TIR and EPICS data to figure out current/previous laser state
bool comptonStatus::SetLaserState(){
  bool transitioned = false;
  //get threshold parameters used for cavity power
  float rtcavpowon = cavityPowerOnMin;
  float rtcavpowoff = cavityPowerOffMax;
  // Set previous laser state
  previousLaserState = currentLaserState;
  lastKnownLaserState = knownLaserState;
  if(cavPowerCalibrated<rtcavpowoff){
    laser_on=0;
  }else if(cavPowerCalibrated>rtcavpowon){
    laser_on=1;
  }else{
    laser_on=-1;  //unkown state (in transistion):
  }
  if(laser_on==-1){
    currentLaserState=LASER_UNKNOWN;
  }else if(rtcavpol==0){
    if(laser_on==1){
      currentLaserState=LASER_LEFT;
    }else{
      currentLaserState=LASER_LEFTOFF;
    }
  }else{
    if(laser_on==1){
      currentLaserState=LASER_RIGHT;
    }else{
      currentLaserState=LASER_RIGHTOFF;
    }
  }    
  // Set new current laser state
  // Due to the unreliability of rtcavpow, we are eliminating it from the test
  // if(cavpow_burp){
  //   currentLaserState=LASER_UNKNOWN;
  // }
  // else if (cavPowerCalibrated < rtcavpowoff){   // Cavity is off!
  //     if(rtcavpol > 0.5){
  //     currentLaserState = LASER_LEFTOFF;
  //   } else {
  //     currentLaserState = LASER_RIGHTOFF;
  //   }
  // } else if (cavPowerCalibrated < rtcavpowon ){ 
  //   currentLaserState = LASER_UNKNOWN; // intermediate cavity power
  // } else {
  //   if ( rtcavpol > 0.5){ 
  //     currentLaserState = LASER_LEFT;
  //   } else {
  //     currentLaserState = LASER_RIGHT;
  //   }
  // }	

  // Now let's compare to the previous laser state

  if(currentLaserState != LASER_UNKNOWN){
    knownLaserState = currentLaserState;
  }
  if (knownLaserState != lastKnownLaserState ){
    transitioned = true;
    laser_good = 0;
    laserstate_clock_start=clockscaler;
    laser_count++;
  } else {
    laserstate_clock= (float)(clockscaler-laserstate_clock_start)/4.e7;
  }
  return transitioned;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Use TIR and EPICS data to figure out current/previous laser state
// private routine
bool comptonStatus::getEpicsValue(THaEpics *epics, const char* tag, float* pReturn, int verbose){
        if (epics->IsLoaded(tag)){
	  *pReturn =epics->GetData(tag);
	  if(verbose>1) printf("%30s = %f\n",tag,*pReturn);
	  //     *pReturn = epics->GetData(tag);
	}else{
	  if(verbose>0) printf("%s Not Found \n",tag);
	}
	return true;
}
void comptonStatus::SetEpicsDefaults(){
  //EPICS values aren't useful until first EPICS event is encounted
  //Set them to someting silly first.
  //Code should really check to see that epics_evnum>=0
  Double_t defValue=-10.;
  epics_PMTRate=defValue;
  epics_PMTRateHigh=defValue;
  epics_tablePosX=defValue;
  epics_tablePosY=defValue;
  epics_hacbmf=defValue;
  epics_s1=defValue;
  epics_s2=defValue;
  epics_cavpow=defValue;
  //
  defValue=-10.0;   //something silly but will showup  in histogram
  epics_aPosX=defValue;
  epics_aPosY=defValue;
  epics_bPosX=defValue;
  epics_bPosY=defValue;
  epics_Thermo1=defValue;
  epics_Thermo2=defValue;
  epics_TimeStamp=0;
 return;
}

int comptonStatus::UnpackEpics(THaEpics *epics, uint32_t* codadata){
      TString spol;
      int evtype = codadata[1]>>16;
      int evnum = codadata[4];
      /*      printf("DEBUG UnpackEpics evntype=%d, evnum=%d\n",evtype,evnum);
      for (int i=0; i<20; i++){
	printf("%2d %10u 0x%08x  ",i,codadata[i],codadata[i]);
	uint32_t tmp=codadata[i];
	for (int j=0; j<4; j++){
	  printf(" %c ",tmp&0xFF);
	  tmp=tmp>>8;
	}
	printf("\n");
      }
      */
      countEpics++;
      countMPSsinceEPICS=0;  //counter for MPS events since last EPICS event
      //  cout << "\nEvent #" << evnum << ": event type is " << evtype;
      if (evtype == 131)		// EPICS event
      {
        epics->LoadData(codadata,evnum);
	bool valid;
	int verbose=0;
	epics_evnum=evnum;
	valid=getEpicsValue(epics,"ComptonCentralRate",&epics_PMTRate,verbose);
	valid=getEpicsValue(epics,"ComptonCentralRateHigh",&epics_PMTRateHigh,verbose);
	valid=getEpicsValue(epics,"COM_DETPH_XCPOSai",&epics_tablePosX,verbose);
	valid=getEpicsValue(epics,"COM_DETPH_YCPOSai",&epics_tablePosY,verbose);
	valid=getEpicsValue(epics,"COMPTON_PW1PCAV_ca",&epics_cavpow,verbose);
	//	valid=getEpicsValue(epics,"COMPTON_SU_POLAR_mo",%epics_spol,verbose);
	valid=getEpicsValue(epics,"hac_bcm_average",&epics_hacbmf,verbose);
	valid=getEpicsValue(epics,"COMPTON_CAVPOLAR_ca",&epics_cavpolpercent,verbose);
	valid=getEpicsValue(epics,"COMPTON_PW1R_S1_ca",&epics_s1,verbose);
	valid=getEpicsValue(epics,"COMPTON_PW1R_S2_ca",&epics_s2,verbose);

	valid=getEpicsValue(epics,"IPM1P02A.YPOS",&epics_aPosY);
	valid=getEpicsValue(epics,"IPM1P02A.XPOS",&epics_aPosX);
	valid=getEpicsValue(epics,"IPM1P02B.YPOS",&epics_bPosY);
	valid=getEpicsValue(epics,"IPM1P02B.XPOS",&epics_bPosX);
	
	valid=getEpicsValue(epics,"HaComptonSIM900_P2T1",&epics_Thermo1);
	valid=getEpicsValue(epics,"HaComptonSIM900_P2T2",&epics_Thermo2);
	epics_TimeStamp = epics->GetTimeStamp("HaComptonSIM900_P2T1");

        if (epics->IsLoaded("IGL1I00OD16_16")){
          spol = epics->GetString("IGL1I00OD16_16");
          epics_ihwp_in =1;
          if (spol.Strip(TString::kTrailing, ' ') == "OUT") epics_cavpoldir = 0;
	}

        if (epics->IsLoaded("COMPTON_SU_POLAR_mo")){
          spol = epics->GetString("COMPTON_SU_POLAR_mo");
          epics_cavpoldir = -1;
          if (spol.Strip(TString::kTrailing, ' ') == "RIGHT") epics_cavpoldir = 1;
	}


	//	valid=getEpicsValue(epics,"MF2b1461_11ch1property.F",&epics_crystalHV,verbose);
	//	valid=getEpicsValue(epics,"MF2b1461_11ch1property.E",&epics_crystalCurrent,verbose);
	/*
        if (epics->IsLoaded("COMPTON_SU_POLAR_mo")){
          spol = epics->GetString("COMPTON_SU_POLAR_mo");
          cavpoldir = -1;
          if (spol.Strip(TString::kTrailing, ' ') == "RIGHT") cavpoldir = 1;
        }   // end cavpol loaded

        }
        if (epics->IsLoaded("MF2b1461_11ch1property.F")){
          crystalHV = epics->GetData("MF2b1461_11ch1property.F");
        }
        if (epics->IsLoaded("MF2b1461_11ch1property.E")){
          crystalCurrent = epics->GetData("MF2b1461_11ch1property.E");	
        }
        if (epics->IsLoaded("MF2b1461_11ch8property.F")){
          verticalFingerHV = epics->GetData("MF2b1461_11ch8property.F");
        }
        if (epics->IsLoaded("MF2b1461_11ch8property.E")){
          verticalFingerCurrent = epics->GetData("MF2b1461_11ch8property.E");
        }
        if (epics->IsLoaded("MF2b1461_11ch7property.F")){
          horizontalFingerHV = epics->GetData("MF2b1461_11ch7property.F");
        }
        if (epics->IsLoaded("MF2b1461_11ch7property.E")){
          horizontalFingerCurrent = epics->GetData("MF2b1461_11ch7property.E");
	*/
      }
      // Fill the EPICS tree
      epicsWiseTree->Fill();
      /*
      if(countMPSsinceEPICS%500==0&&countMPSsinceEPICS>0)
        epicsWiseTree->AutoSave("SaveSelf");
        */
      return 0;
}
int comptonStatus::EpicsDebugDump(){
      cout << "\n----------EPICS DATA DUMP-----------";
      cout <<"\nEpics Event Number: "<<epics_evnum;
      cout <<"\nEpics Event TimeStamp: "<<epics_TimeStamp;
      cout <<"\nEpics Counter "<<countEpics;
      cout << "\nCavity power: " << epics_cavpow;
      cout << "\nCavity polarization direction: " << epics_cavpoldir;
      cout << "\nCavity polarization percentage: " << epics_cavpolpercent;
      cout << "\nHigh voltage on crystal PMT: " << epics_crystalHV;
      cout << "\nCurrent on crystal PMT: " << epics_crystalCurrent;
      cout << "\nAverage beam current: " << epics_hacbmf;
      cout << "\nInsertable half-waveplate state: " << epics_ihwp_in;
      cout << "\nPhDet Temperature1: " << epics_Thermo1;
      cout << "\nPhDet Temperature2: " << epics_Thermo2;
      cout << endl;
       
      return 0;
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
float comptonStatus::ComputeBPMPosition( float x_pos, float x_neg,
    float x_pos_ped, float x_neg_ped, float alpha, float sensitivity,
    float freq_conversion){
  x_pos *= freq_conversion;
  x_neg *= freq_conversion;
  float denom = ((x_pos-x_pos_ped) + alpha*(x_neg-x_neg_ped));
  if(denom!=0)
    return sensitivity*((x_pos-x_pos_ped) - alpha*(x_neg-x_neg_ped))/denom;
  return 1e6;
}
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void comptonStatus::ComputeBPMPositionLab(float xrot, float yrot,
      float sintheta, float costheta, float xoff, float yoff,
      float &xlab, float &ylab){
  xlab = xrot*costheta-yrot*sintheta-xoff;
  ylab = xrot*sintheta+yrot*costheta-yoff;
}
