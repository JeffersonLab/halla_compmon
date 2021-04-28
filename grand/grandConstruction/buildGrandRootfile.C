#include "buildGrandRootfile.h"
#include "plot.h"
#include "burst.h"
#include "../../online/runs.h"

using namespace std;

void readKeysFile(TString expt){
  ifstream keysfile(Form("%s/%s_rmsCut.key", getenv("COMPMON_MAPS"), expt.Data()));
  if(keysfile.is_open()){
    string line;
    while(getline(keysfile, line)){
      vector<Float_t> keyRange; stringstream ss(line);
      for(Float_t i; ss >> i;){
        keyRange.push_back(i);
        if(ss.peek() == ',')
          ss.ignore();
      }
      keys.push_back(keyRange);
    }
  }
  keysfile.close();
}

void readCorrsFile(TString expt){
  ifstream keysfile(Form("%s/%s_corrSlopes.key", getenv("COMPMON_MAPS"), expt.Data()));
  if(keysfile.is_open()){
    string line;
    while(getline(keysfile, line)){
      vector<Float_t> keyRange; stringstream ss(line);
      for(Float_t i; ss >> i;){
        keyRange.push_back(i);
        if(ss.peek() == ',')
          ss.ignore();
      }
      corrs.push_back(keyRange);
    }
  }
  keysfile.close();
}

Bool_t acceptCycleRMS(Int_t runNum){
  if(keys.size() == 0){
    acc0OnMode = 0.0; acc0OnOffset = 0.0;
    acc0OffMode = 0.0; acc0OffOffset = 0.0;
    return true;
  }
  else{
    for(Int_t i = 0; i < keys.size(); i++){
      if(runNum >= (Int_t)keys[i][0] && runNum <= (Int_t)keys[i][1]){
        acc0OnMode = keys[i][3]; acc0OnOffset = keys[i][4];
        acc0OffMode = keys[i][2]; acc0OffOffset = keys[i][4];
        Float_t offLimit = keys[i][2] + keys[i][4];
        Float_t onLimit = keys[i][3] + keys[i][4];
        return cycMPSData[1].rms < offLimit && cycMPSData[2].rms < offLimit && cycMPSData[0].rms < onLimit;
      }
    }
    return true;
  }
}

Bool_t acceptCycleDoubleDiff(){
  Float_t diff = asym0LasOff1.mean - asym0LasOff2.mean;
  Float_t diffErr = TMath::Sqrt(TMath::Power(asym0LasOff1.meanErr, 2) + TMath::Power(asym0LasOff2.meanErr, 2));
  return TMath::Abs(diff*1.0/diffErr) < 3.0;
}

Bool_t acceptCycle(Int_t runNum){
  //return acceptCycleRMS(runNum);
  Int_t rmsCut = (Int_t)acceptCycleRMS(runNum);
  Int_t signalSizeCut = (Int_t)(cycMPSData[0].mean - (cycMPSData[1].mean + cycMPSData[2].mean)/2.0 > 0.7);
  Int_t doubleDiffCut = (Int_t)acceptCycleDoubleDiff();
  Int_t polErrCut = (Int_t)(pol0.meanErr < 0.3);
  Int_t backCut = (Int_t)(cycMPSData[16].mean < 0.8 && cycMPSData[17].mean < 2.8);
  if(runNum < 4330){
    polErrCut = (Int_t)(pol0.meanErr < 0.8);
  }
  Int_t offAsymCut = (Int_t)(TMath::Abs(asym0LasOff.mean) < 0.001656);
  cycleCut = 0;
  cycleCut += 0x01*(Int_t)(!rmsCut);
  cycleCut += 0x02*(Int_t)(!signalSizeCut);
  cycleCut += 0x04*(Int_t)(!doubleDiffCut);
  cycleCut += 0x08*(Int_t)(!polErrCut);
  //cycleCut += 0x10*(Int_t)(!backCut);

  return cycleCut == 0;
}

Bool_t isCloseTo(Float_t num, Float_t ref, Float_t pct=0.01){
  Bool_t val = TMath::Abs(num)>=(1.0-pct)*TMath::Abs(ref) && TMath::Abs(num)<=(1.0+pct)*TMath::Abs(ref);
  //printf("Close to checks: Num: %f, Ref: %f, Result: %s\n", num, ref, val ? "true" : "false");
  return val;
}

Float_t getLaserPol(Int_t snlNum, Float_t qw1Deg, Float_t hw1Deg){
  //PREX Conditions
  //printf("QW1 Deg: %.1f; HW1 Deg: %.1f\n", qw1Deg, hw1Deg);
  if(snlNum < 100 || snlNum == 500){
    if(isCloseTo(qw1Deg, 49.2)){
      if(isCloseTo(hw1Deg, 0.2)){return 0.9999;}
      else if(isCloseTo(hw1Deg, 15.2)){return 0.9880;}
      else if(isCloseTo(hw1Deg, 31.2)){return 0.9595;}
      else{return 0.0;}
    }
    else if(isCloseTo(qw1Deg, 47.7)){
      if(isCloseTo(hw1Deg, 0.2)){return 0.0;}
      else if(isCloseTo(hw1Deg, 19.1)){return 0.9887;}
      else{return 0.0;}
    }
    else{return 0.0;}
  }
  //CREX Conditions
  else{
    return 1.0;
  }
}

void calcCycLaserPol(Int_t snlNum){
  Float_t qw1Deg = 0.0;
  if(snlNum < 100 || snlNum == 500)
    qw1Deg = fmod((-0.004*cycQW1 - 1390.768), 360);
  else
    qw1Deg = fmod((-0.004*cycQW1), 360);
  Float_t hw1Deg = fmod((-0.004*cycHW1), 360);
  
  cycLaserPol.mean = getLaserPol(snlNum, qw1Deg, hw1Deg);
  if(snlNum < 100 || snlNum == 500){cycLaserPol.meanErr = 0.0038*cycLaserPol.mean;}
  else{cycLaserPol.meanErr = 0.0;}
}

void calcRunLaserPol(Int_t snlNum){
  Float_t qw1Deg = 0.0;
  if(snlNum < 100 || snlNum == 500)
    qw1Deg = fmod((-0.004*runEpcData[2].mean - 1390.768), 360);
  else
    qw1Deg = fmod((-0.004*runEpcData[2].mean), 360);
  Float_t hw1Deg = fmod((-0.004*runEpcData[3].mean), 360);
  
  runLaserPol.mean = getLaserPol(snlNum, qw1Deg, hw1Deg);
  if(snlNum < 100 || snlNum == 500){runLaserPol.meanErr = 0.0038*runLaserPol.mean;}
  else{runLaserPol.meanErr = 0.0;}
}

void calcSnailLaserPol(Int_t snlNum){
  Float_t qw1Deg = 0.0;
  if(snlNum < 100 || snlNum == 500)
    qw1Deg = fmod((-0.004*qw1 - 1390.768), 360);
  else
    qw1Deg = fmod((-0.004*qw1), 360);
  Float_t hw1Deg = fmod((-0.004*hw1), 360);
  
  laserPol = getLaserPol(snlNum, qw1Deg, hw1Deg);
  //if(snlNum < 100 || snlNum == 500){laserPol.meanErr = 0.0038*laserPol.mean;}
  //else{laserPol.meanErr = 0.0;}
}

void calcAsymCorr(Int_t snlNum, Int_t runNum){
  Float_t p0 = 0.0; Float_t p1 = 0.0;
  if(snlNum < 100 || snlNum == 500){
    p0 = 0.00521544;
    p1 = -0.000602352;
  }
  Float_t asymDiffCorr = 0.0;
  Float_t asymDiffCorrErr = 0.0;
  for(Int_t i = 0; i < corrs.size(); i++){
    if(runNum >= corrs[i][0] && runNum <= corrs[i][1]){
      Float_t diffAvg = (cycQrtData[33].mean + cycQrtData[35].mean)/2.0;
      asymDiffCorr = corrs[i][2]*diffAvg;
      asymDiffCorrErr = corrs[i][3]*diffAvg; 
    }
  }
  
  Float_t asymCorr = asym0NGC.mean - asymDiffCorr;
  Float_t asymCorrErr = TMath::Sqrt(TMath::Power(asym0NGC.meanErr, 2) + TMath::Power(asymDiffCorrErr, 2));
  Float_t meanAcc0LasOn = cycMPSData[0].mean;
  Float_t meanAcc0LasOff = (cycMPSData[1].mean + cycMPSData[2].mean)/2.;
  alphaOn = p0 + meanAcc0LasOn*p1;
  alphaOff = p0 + meanAcc0LasOff*p1;
  Float_t norm = meanSum0LasOn - meanSum0LasOff;
  Float_t diffCorr = (meanDiff0LasOn*alphaOn - meanDiff0LasOff*alphaOff)/norm;
  Float_t sumCorr = (meanSum0LasOn*alphaOn - meanSum0LasOff*alphaOff)/norm;
  asym0.mean = (asymCorr + diffCorr)/(1 + sumCorr);
  asym0.meanErr = asymCorrErr/(1 + sumCorr);
}

void calcBurstAsymCorr(Int_t snlNum, Int_t runNum){
  //This function assumes calcAsymCorr has
  //already been called once in the cycle

  Float_t burstDiff0LasOn = cycBurstData[2][0].mean;
  Float_t burstDiff0LasOff = cycBurstData[2][3].mean;
  Float_t burstSum0LasOn = cycBurstData[3][0].mean;
  Float_t burstSum0LasOff = cycBurstData[3][3].mean;
  Float_t asymDiffCorr = 0.0;
  Float_t asymDiffCorrErr = 0.0;
  for(Int_t i = 0; i < corrs.size(); i++){
    if(runNum >= corrs[i][0] && runNum <= corrs[i][1]){
      Float_t diffAvg = (cycQrtData[33].mean + cycQrtData[35].mean)/2.0;
      asymDiffCorr = corrs[i][2]*diffAvg;
      asymDiffCorrErr = corrs[i][3]*diffAvg;
    }
  }
  
  Float_t asymCorr = cycBurstComboData[0].mean - asymDiffCorr;
  Float_t asymCorrErr = TMath::Sqrt(TMath::Power(cycBurstComboData[0].meanErr, 2) + TMath::Power(asymDiffCorrErr, 2));
  Float_t norm = burstSum0LasOn - burstSum0LasOff;
  Float_t diffCorr = (burstDiff0LasOn*alphaOn - burstDiff0LasOff*alphaOff)/norm;
  Float_t sumCorr = (burstSum0LasOn*alphaOn - burstSum0LasOff*alphaOff)/norm;
  cycBurstComboData[1].mean = (asymCorr + diffCorr)/(1 + sumCorr);
  cycBurstComboData[1].meanErr = asymCorrErr/(1 + sumCorr);
}

void calcCyclePol(Int_t snlNum, Int_t runNum){
  calcCycLaserPol(snlNum);
  calcCollOffset(snlNum, runNum);
  setAnalyzingPower(snlNum, runNum);

  Int_t diff0On = 0; Int_t sum0On = 1; Int_t diff0Off1 = 2; Int_t sum0Off1 = 3; Int_t diff0Off2 = 4; Int_t sum0Off2 = 5;
  Int_t diff4On = 6; Int_t sum4On = 7; Int_t diff4Off1 = 8; Int_t sum4Off1 = 9; Int_t diff4Off2 = 10; Int_t sum4Off2 = 11;

  meanDiff0LasOn = cycQrtCalc[diff0On].mean; meanErrDiff0LasOn = cycQrtCalc[diff0On].meanErr; 
  meanSum0LasOn = cycQrtCalc[sum0On].mean; meanErrSum0LasOn = cycQrtCalc[sum0On].meanErr;
  meanDiff0LasOff1 = cycQrtCalc[diff0Off1].mean; meanErrDiff0LasOff1 = cycQrtCalc[diff0Off1].meanErr; 
  meanSum0LasOff1 = cycQrtCalc[sum0Off1].mean; meanErrSum0LasOff1 = cycQrtCalc[sum0Off1].meanErr; 
  meanDiff0LasOff2 = cycQrtCalc[diff0Off2].mean; meanErrDiff0LasOff2 = cycQrtCalc[diff0Off2].meanErr;
  meanSum0LasOff2 = cycQrtCalc[sum0Off2].mean; meanErrSum0LasOff2 = cycQrtCalc[sum0Off2].meanErr;
  meanDiff0LasOff = (meanDiff0LasOff1 + meanDiff0LasOff2)/2.0;
  meanErrDiff0LasOff = TMath::Sqrt(TMath::Power(meanErrDiff0LasOff1, 2) + TMath::Power(meanErrDiff0LasOff2, 2))/2.0;
  meanSum0LasOff = (meanSum0LasOff1 + meanSum0LasOff2)/2.;
  meanErrSum0LasOff = TMath::Sqrt(TMath::Power(meanErrSum0LasOff1, 2) + TMath::Power(meanErrSum0LasOff2, 2))/2.0;
  sigSubSum0.mean = meanSum0LasOn - meanSum0LasOff;
  sigSubSum0.meanErr = TMath::Sqrt(TMath::Power(meanErrSum0LasOn, 2) + TMath::Power(meanErrSum0LasOff, 2));
  asym0LasOnAlt.mean = meanDiff0LasOn/sigSubSum0.mean;
  asym0LasOnAlt.meanErr = TMath::Abs(asym0LasOnAlt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOn/meanDiff0LasOn, 2) + TMath::Power(sigSubSum0.meanErr/sigSubSum0.mean, 2));
  asym0LasOffAlt.mean = meanDiff0LasOff/sigSubSum0.mean;
  asym0LasOffAlt.meanErr = TMath::Abs(asym0LasOffAlt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOff/meanDiff0LasOff, 2) + TMath::Power(sigSubSum0.meanErr/sigSubSum0.mean, 2));
  asym0LasOff1Alt.mean = meanDiff0LasOff1/sigSubSum0.mean;
  asym0LasOff1Alt.meanErr = TMath::Abs(asym0LasOff1Alt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOff1/meanDiff0LasOff1, 2) + TMath::Power(sigSubSum0.meanErr/sigSubSum0.mean, 2));
  asym0LasOff2Alt.mean = meanDiff0LasOff2/sigSubSum0.mean;
  asym0LasOff2Alt.meanErr = TMath::Abs(asym0LasOff2Alt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff0LasOff2/meanDiff0LasOff2, 2) + TMath::Power(sigSubSum0.meanErr/sigSubSum0.mean, 2));
  //asym0LasOff1.mean = sigSubAsym0LasOff1.mean; asym0LasOff1.meanErr = sigSubAsym0LasOff1.meanErr;
  //asym0LasOff2.mean = sigSubAsym0LasOff2.mean; asym0LasOff2.meanErr = sigSubAsym0LasOff2.meanErr;
  asym0LasOff.mean = (asym0LasOff1.mean + asym0LasOff2.mean)/2.0;
  asym0LasOff.meanErr = TMath::Sqrt(TMath::Power(asym0LasOff1.meanErr, 2) + TMath::Power(asym0LasOff2.meanErr, 2))/2.0;
  asym0NGC.mean = asym0LasOn.mean - asym0LasOff.mean; 
  asym0NGC.meanErr = TMath::Sqrt(TMath::Power(asym0LasOn.meanErr, 2) + TMath::Power(asym0LasOff.meanErr, 2));
  calcAsymCorr(snlNum, runNum);
  Float_t norm = anPow.mean*cycLaserPol.mean;
  Float_t normErr = TMath::Abs(norm)*TMath::Sqrt(TMath::Power(anPow.meanErr/anPow.mean, 2) + TMath::Power(cycLaserPol.meanErr/cycLaserPol.mean, 2));
  pol0.mean = asym0.mean/norm;
  pol0.meanErr = TMath::Abs(pol0.mean)*TMath::Sqrt(TMath::Power(asym0.meanErr/asym0.mean, 2) + TMath::Power(normErr/norm, 2));

  meanDiff4LasOn = cycQrtCalc[diff4On].mean; meanErrDiff4LasOn = cycQrtCalc[diff4On].meanErr; 
  meanSum4LasOn = cycQrtCalc[sum4On].mean; meanErrSum4LasOn = cycQrtCalc[sum4On].meanErr;
  meanDiff4LasOff1 = cycQrtCalc[diff4Off1].mean; meanErrDiff4LasOff = cycQrtCalc[diff4Off1].meanErr; 
  meanSum4LasOff1 = cycQrtCalc[sum4Off1].mean; meanErrSum4LasOff1 = cycQrtCalc[sum4Off1].meanErr;
  meanDiff4LasOff2 = cycQrtCalc[diff4Off2].mean; meanErrDiff4LasOff = cycQrtCalc[diff4Off2].meanErr;
  meanSum4LasOff2 = cycQrtCalc[sum4Off2].mean; meanErrSum4LasOff2 = cycQrtCalc[sum4Off2].meanErr;
  meanDiff4LasOff = (meanDiff4LasOff1 + meanDiff4LasOff2)/2.0;
  meanErrDiff4LasOff = TMath::Sqrt(TMath::Power(meanErrDiff4LasOff1, 2) + TMath::Power(meanErrDiff4LasOff2, 2))/2.0;
  meanSum4LasOff = (meanSum4LasOff1 + meanSum4LasOff2)/2.;
  meanErrSum4LasOff = TMath::Sqrt(TMath::Power(meanErrSum4LasOff1, 2) + TMath::Power(meanErrSum4LasOff2, 2))/2.0;
  sigSubSum4.mean = meanSum4LasOn - meanSum4LasOff;
  sigSubSum4.meanErr = TMath::Sqrt(TMath::Power(meanErrSum4LasOn, 2) + TMath::Power(meanErrSum4LasOff, 2));
  asym4LasOnAlt.mean = meanDiff4LasOn/sigSubSum4.mean;
  asym4LasOnAlt.meanErr = TMath::Abs(asym4LasOnAlt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOn/meanDiff4LasOn, 2) + TMath::Power(sigSubSum4.meanErr/sigSubSum4.mean, 2));
  asym4LasOffAlt.mean = meanDiff4LasOff/sigSubSum4.mean;
  asym4LasOffAlt.meanErr = TMath::Abs(asym4LasOffAlt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOff/meanDiff4LasOff, 2) + TMath::Power(sigSubSum4.meanErr/sigSubSum4.mean, 2));
  asym4LasOff1Alt.mean = meanDiff4LasOff1/sigSubSum4.mean;
  asym4LasOff1Alt.meanErr = TMath::Abs(asym4LasOff1Alt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOff1/meanDiff4LasOff1, 2) + TMath::Power(sigSubSum4.meanErr/sigSubSum4.mean, 2));
  asym4LasOff2Alt.mean = meanDiff4LasOff2/sigSubSum4.mean;
  asym4LasOff2Alt.meanErr = TMath::Abs(asym4LasOff2Alt.mean)*TMath::Sqrt(TMath::Power(meanErrDiff4LasOff2/meanDiff4LasOff2, 2) + TMath::Power(sigSubSum4.meanErr/sigSubSum4.mean, 2));
  //asym4LasOff1.mean = sigSubAsym4LasOff1.mean; asym4LasOff1.meanErr = sigSubAsym4LasOff1.meanErr;
  //asym4LasOff2.mean = sigSubAsym4LasOff2.mean; asym4LasOff2.meanErr = sigSubAsym4LasOff2.meanErr;
  asym4LasOff.mean = (asym4LasOff1.mean + asym4LasOff2.mean)/2.0;
  asym4LasOff.meanErr = TMath::Sqrt(TMath::Power(asym4LasOff1.meanErr, 2) + TMath::Power(asym4LasOff2.meanErr, 2))/2.0;
  asym4NGC.mean = asym4LasOn.mean - asym4LasOff.mean; 
  asym4NGC.meanErr = TMath::Sqrt(TMath::Power(asym4LasOn.meanErr, 2) + TMath::Power(asym4LasOff.meanErr, 2));
  asym4.mean = asym4NGC.mean;
  asym4.meanErr = asym4NGC.meanErr;
  pol4.mean = asym4.mean/norm;
  pol4.meanErr = TMath::Abs(pol4.mean)*TMath::Sqrt(TMath::Power(asym4.meanErr/asym4.mean, 2) + TMath::Power(normErr/norm, 2));

  if(acceptCycle(runNum)){
    Float_t cycVars[nPolVars] = {asym0.mean, asym0NGC.mean, asym0LasOn.mean, asym0LasOff.mean, pol0.mean};
    Float_t cycErrs[nPolVars] = {asym0.meanErr, asym0NGC.meanErr, asym0LasOn.meanErr, asym0LasOff.meanErr, pol0.meanErr};
    for(Int_t i = 0; i < nPolVars; i++){
      snlPol0Avg[i].push_back(cycVars[i]); snlPol0Err[i].push_back(cycErrs[i]);
      runPol0Avg[i].push_back(cycVars[i]); runPol0Err[i].push_back(cycErrs[i]);
    }
    //runAsym0Avg.push_back(asym0.mean); runAsym0Err.push_back(asym0.meanErr);
    //runAsym0NGCAvg.push_back(asym0NGC.mean); runAsym0NGCErr.push_back(asym0NGC.meanErr);
    //runAsym0OnAvg.push_back(asym0LasOn.mean); runAsym0OnErr.push_back(asym0LasOn.meanErr);
    //runAsym0OffAvg.push_back(asym0LasOff.mean); runAsym0OffErr.push_back(asym0LasOff.meanErr);
    //runPol0Avg.push_back(pol0.mean); runPol0Err.push_back(pol0.meanErr);
    //snlAsym0Avg.push_back(asym0.mean); snlAsym0Err.push_back(asym0.meanErr);
    //snlAsym0NGCAvg.push_back(asym0NGC.mean); snlAsym0NGCErr.push_back(asym0NGC.meanErr);
    //snlAsym0OnAvg.push_back(asym0LasOn.mean); snlAsym0OnErr.push_back(asym0LasOn.meanErr);
    //snlAsym0OffAvg.push_back(asym0LasOff.mean); snlAsym0OffErr.push_back(asym0LasOff.meanErr);
    //snlPol0Avg.push_back(pol0.mean); snlPol0Err.push_back(pol0.meanErr);
    //snlAsym4Avg.push_back(asym4.mean); snlAsym4Err.push_back(asym4.meanErr);
    //snlPol4Avg.push_back(pol4.mean); snlPol4Err.push_back(pol4.meanErr);
    numCyclesAcc++;
  }
  else{
    printf("      Rejecting cycle for polarization with error %i!\n", cycleCut);
  }
}

void calcCycleBurstPol(Int_t snlNum, Int_t runNum, Int_t cNum, TFile *plotfile){
  Int_t nBursts = 0;
  for(Int_t i = 0; i < burstVars; i++){
    for(Int_t j = 0; j < burstStates; j++){
      TString hName = Form("h%i.%i_%s%s", runNum, cNum, burstTitles[i].Data(), burstLas[j].Data());
      TF1 *fit = new TF1(Form("f%i.%i_%s%s", runNum, cNum, burstTitles[i].Data(), burstLas[j].Data()), "pol0");
      TH1F *h = (TH1F *)plotfile->Get(hName.Data());
      h->Fit(fit, "Q", "goff", 0, h->GetNbinsX());
      cycBurstData[i][j].mean = fit->GetParameter(0);
      cycBurstData[i][j].meanErr = fit->GetParError(0);
      cycBurstData[i][j].Chi2 = fit->GetChisquare();
      cycBurstData[i][j].NDF = fit->GetNDF();
      if(i==0 && j==0) nBursts = h->GetNbinsX();
      delete fit;
    }
  }

  TF1 *fAsyms[burstAsyms];
  for(Int_t i = 0; i < burstAsyms; i++){
    TString hName = Form("h%i.%i_%s", runNum, cNum, burstAsymTitles[i].Data());
    TF1 *fit = new TF1(Form("f%i.%i_%s", runNum, cNum, burstAsymTitles[i].Data()), "pol0");
    TH1F *h = (TH1F *)plotfile->Get(hName.Data());
    h->Fit(fit, "Q", "goff", 0, h->GetNbinsX());
    cycBurstAsymData[i].mean = fit->GetParameter(0);
    cycBurstAsymData[i].meanErr = fit->GetParError(0);
    cycBurstAsymData[i].Chi2 = fit->GetChisquare();
    cycBurstAsymData[i].NDF = fit->GetNDF();
    fAsyms[i] = fit;
  }

  Float_t meanAsym0LasOn = fAsyms[0]->GetParameter(0);
  Float_t meanAsym0LasOff = fAsyms[3]->GetParameter(0);
  Float_t meanErrAsym0LasOn = fAsyms[0]->GetParError(0);
  Float_t meanErrAsym0LasOff = fAsyms[3]->GetParError(0);
  Float_t meanAsym0NGC = fAsyms[0]->GetParameter(0) - fAsyms[3]->GetParameter(0);
  Float_t meanErrAsym0NGC = TMath::Sqrt(TMath::Power(fAsyms[0]->GetParError(0), 2) + TMath::Power(fAsyms[3]->GetParError(0), 2));
  cycBurstComboData[0].mean = meanAsym0NGC;
  cycBurstComboData[0].meanErr = meanErrAsym0NGC;
  calcBurstAsymCorr(snlNum, runNum);
  Float_t norm = anPow.mean*cycLaserPol.mean;
  Float_t normErr = TMath::Abs(norm)*TMath::Sqrt(TMath::Power(anPow.meanErr/anPow.mean, 2) + TMath::Power(cycLaserPol.meanErr/cycLaserPol.mean, 2));
  cycBurstComboData[2].mean = cycBurstComboData[1].mean/norm;
  cycBurstComboData[2].meanErr = TMath::Abs(cycBurstComboData[2].mean) * 
                                 TMath::Sqrt(TMath::Power(cycBurstComboData[1].meanErr/cycBurstComboData[1].mean, 2) + TMath::Power(normErr/norm, 2));

  if(acceptCycle(runNum)){
    Float_t burstAvg[nPolVars] = {cycBurstComboData[1].mean, meanAsym0NGC, meanAsym0LasOn, meanAsym0LasOff, cycBurstComboData[2].mean};
    Float_t burstErr[nPolVars] = {cycBurstComboData[1].meanErr, meanErrAsym0NGC, meanErrAsym0LasOn, meanErrAsym0LasOff, cycBurstComboData[2].meanErr};
    for(Int_t i = 0; i < nPolVars; i++){
      snlBurstAvg[i].push_back(burstAvg[i]); snlBurstErr[i].push_back(burstErr[i]);
      runBurstAvg[i].push_back(burstAvg[i]); runBurstErr[i].push_back(burstErr[i]);
    }
  }
  
  //printf("Pol Info: %.4f +/- %.4f, with norm: %.4f\n", cycBurstComboData[2].mean, cycBurstComboData[2].meanErr, norm);
  for(Int_t i = 0; i < burstAsyms; i++){
    delete fAsyms[i];
  }
}

void calcCyclePedestals(TFile *plotfile){
  TString hNameF = Form("h%i.%i_pedF", runNum, cycleNum);
  TString hNameL = Form("h%i.%i_pedL", runNum, cycleNum);
  TString fNameF1 = Form("h%i.%i_pedF1", runNum, cycleNum);
  TString fNameF2 = Form("h%i.%i_pedF2", runNum, cycleNum);
  TString fNameL1 = Form("h%i.%i_pedL1", runNum, cycleNum);
  TString fNameL2 = Form("h%i.%i_pedL2", runNum, cycleNum);

  TH1F *hF = (TH1F *)plotfile->Get(hNameF.Data());
  TH1F *hL = (TH1F *)plotfile->Get(hNameL.Data());

  TF1 *fF1 = new TF1(fNameF1.Data(), "gaus");
  TF1 *fF2 = new TF1(fNameF2.Data(), "gaus");
  TF1 *fL1 = new TF1(fNameL1.Data(), "gaus");
  TF1 *fL2 = new TF1(fNameL2.Data(), "gaus");

  Float_t meanF1 = hF->GetMean(); Float_t rmsF1 = hF->GetRMS();
  Float_t meanL1 = hL->GetMean(); Float_t rmsL1 = hL->GetRMS();
  hF->Fit(fNameF1.Data(), "Q", "goff", meanF1 - rmsF1, meanF1 + rmsF1);
  hL->Fit(fNameL1.Data(), "Q", "goff", meanL1 - rmsL1, meanL1 + rmsL1);
  Float_t meanF2 = fF1->GetParameter(1); Float_t rmsF2 = fF1->GetParameter(2);
  Float_t meanL2 = fL1->GetParameter(1); Float_t rmsL2 = fF1->GetParameter(2);
  hF->Fit(fNameF2.Data(), "Q", "goff", meanF2 - rmsF2, meanF2 + rmsF2);
  hL->Fit(fNameL2.Data(), "Q", "goff", meanL2 - rmsL2, meanL2 + rmsL2);
  firstOffPedestal.mean = fF2->GetParameter(1); firstOffPedestal.meanErr = fF2->GetParError(1);
  firstOffPedestal.rms  = fF2->GetParameter(2); firstOffPedestal.rmsErr  = fF2->GetParError(2);
  lastOffPedestal.mean  = fL2->GetParameter(1); lastOffPedestal.meanErr  = fL2->GetParError(1);
  lastOffPedestal.rms   = fL2->GetParameter(2); lastOffPedestal.rmsErr   = fL2->GetParError(2);
}

void calcSnailSign(Int_t snailNum){
  Int_t ihwpSign = 1 - 2*ihwp;
  Int_t wienSign = 0;
  if(snailNum < 99){
    if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, -89.9008) && isCloseTo(PhiFG, 86.9014))
      wienSign = 1;
    else if(isCloseTo(HWienAngle, -13.0) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 91.1902))
      wienSign = -1;
  }
  else{
    if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, -90.5996) && isCloseTo(PhiFG, 88.0277))
      wienSign = 1;
    else if(isCloseTo(HWienAngle, -29.6402) && isCloseTo(VWienAngle, 88.0008) && isCloseTo(PhiFG, 89.9558))
      wienSign = -1;
  }
  snailSign = ihwpSign*wienSign;
  //printf("Snail sign: %i\n", snailSign);
}

void calcRunSign(Int_t runNum, Float_t runHWien, Float_t runVWien, Float_t runSolWien, Int_t runIHWP){
  Int_t ihwpSign = 1 - 2*runIHWP;
  Int_t wienSign = 0;
  if(runNum < 4800){
    if(isCloseTo(runHWien, -13.0) && isCloseTo(runVWien, -89.9008) && isCloseTo(runSolWien, 86.9014))
      wienSign = 1;
    else if(isCloseTo(runHWien, -13.0) && isCloseTo(runVWien, 88.0008) && isCloseTo(runSolWien, 91.1902))
      wienSign = -1;
  }
  else{
    if(isCloseTo(runHWien, -29.6402) && isCloseTo(runVWien, -90.5996) && isCloseTo(runSolWien, 88.0277))
      wienSign = 1;
    else if(isCloseTo(runHWien, -29.6402) && isCloseTo(runVWien, 88.0008) && isCloseTo(runSolWien, 89.9558))
      wienSign = -1;
  }
  runSign = ihwpSign*wienSign;
  //printf("Snail sign: %i\n", snailSign);
}

void initCycleTree(TTree *cyc){
  cyc->Branch("snailNum", &snailNum, "snailNum/I");
  cyc->Branch("runNum", &runNum, "runNum/I");
  cyc->Branch("cycleNum", &cycleNum, "cycleNum/I");
  cyc->Branch("firstOffStartMPS", &firstOffStartMPS, "firstOffStartMPS/I");
  cyc->Branch("firstOffEndMPS", &firstOffEndMPS, "firstOffEndMPS/I");
  cyc->Branch("onStartMPS", &onStartMPS, "onStartMPS/I");
  cyc->Branch("onEndMPS", &onEndMPS, "onEndMPS/I");
  cyc->Branch("lastOffStartMPS", &lastOffStartMPS, "lastOffStartMPS/I");
  cyc->Branch("lastOffEndMPS", &lastOffEndMPS, "lastOffEndMPS/I");
  cyc->Branch("cycleTime", &cycleTime, "cycleTime/F");
  cyc->Branch("sign", &cycSign, "sign/I");
  cyc->Branch("qw1", &cycQW1, "qw1/F");
  cyc->Branch("hw1", &cycHW1, "hw1/F");
  cyc->Branch("qw2", &cycQW2, "qw2/F");
  cyc->Branch("ihwp", &cycIHWP, "ihwp/F");
  cyc->Branch("VWienAngle", &cycVWien, "VWienAngle/F");
  cyc->Branch("HWienAngle", &cycHWien, "HWienAngle/F");
  cyc->Branch("PhiFG", &cycSolWien, "PhiFG/F");
  for(Int_t i = 0; i < cycMPSVars; i++){DataVar data; cycMPSData.push_back(data);}
  for(Int_t i = 0; i < cycMPSVars; i++){cyc->Branch(cycMPSTitles[i].Data(), &cycMPSData[i], "mean/F:meanErr/F:rms/F:rmsErr/F");}
  //cyc->Branch("PedestalMeanFirstOff", &firstOffPedestal, "mean/F:meanErr/F:rms/F:rmsErr/F");
  //cyc->Branch("PedestalMeanLastOff", &lastOffPedestal, "mean/F:meanErr/F:rms/F:rmsErr/F");
  for(Int_t i = 0; i < cycQrtVars; i++){DataVar data; cycQrtData.push_back(data);}
  for(Int_t i = 0; i < cycQrtVars; i++){cyc->Branch(cycQrtTitles[i].Data(), &cycQrtData[i], "mean/F:meanErr/F:rms/F:rmsErr/F");}
  cyc->Branch("CollimatorOffset", &collOffset, "mean/F:meanErr/F");
  cyc->Branch("AnalyzingPower", &anPow, "mean/F:meanErr/F");
  cyc->Branch("AlphaLasOn", &alphaOn, "AlphaLasOn/F");
  cyc->Branch("AlphaLasOff", &alphaOff, "AlphaLasOff/F");
  cyc->Branch("SigSubSum0", &sigSubSum0, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOn", &asym0LasOn, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff", &asym0LasOff, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff1", &asym0LasOff1, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff2", &asym0LasOff2, "mean/F:meanErr/F");
  cyc->Branch("Asym0NGC", &asym0NGC, "mean/F:meanErr/F");
  cyc->Branch("Asym0", &asym0, "mean/F:meanErr/F");
  cyc->Branch("Pol0", &pol0, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOnAlt", &asym0LasOnAlt, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOffAlt", &asym0LasOffAlt, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff1Alt", &asym0LasOff1Alt, "mean/F:meanErr/F");
  cyc->Branch("Asym0LasOff2Alt", &asym0LasOff2Alt, "mean/F:meanErr/F");
  cyc->Branch("SigSubSum4", &sigSubSum4, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOn", &asym4LasOn, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff", &asym4LasOff, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff1", &asym4LasOff1, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff2", &asym4LasOff2, "mean/F:meanErr/F");
  cyc->Branch("Asym4NGC", &asym4LasOff2, "mean/F:meanErr/F");
  cyc->Branch("Asym4", &asym4, "mean/F:meanErr/F");
  cyc->Branch("Pol4", &pol4, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOnAlt", &asym4LasOnAlt, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOffAlt", &asym4LasOffAlt, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff1Alt", &asym4LasOff1Alt, "mean/F:meanErr/F");
  cyc->Branch("Asym4LasOff2Alt", &asym4LasOff2Alt, "mean/F:meanErr/F");
  cyc->Branch("CycleCut", &cycleCut, "CycleCut/I");
  cyc->Branch("Acc0OnMode", &acc0OnMode, "Acc0OnMode/F");
  cyc->Branch("Acc0OnOffset", &acc0OnOffset, "Acc0OnOffset/F");
  cyc->Branch("Acc0OffMode", &acc0OffMode, "Acc0OffMode/F");
  cyc->Branch("Acc0OffOffset", &acc0OffOffset, "Acc0OffOffset/F");
  for(Int_t i = 0; i < burstVars; i++){
    vector<FitPolVar> laserState;
    for(Int_t j = 0; j < burstStates; j++){
      FitPolVar data; laserState.push_back(data);
    }
    cycBurstData.push_back(laserState);
  }
  for(Int_t i = 0; i < burstVars; i++){
    for(Int_t j = 0; j < burstStates; j++){
      cyc->Branch(Form("%s%s", burstTitles[i].Data(), burstLas[j].Data()), &cycBurstData[i][j], "mean/F:meanErr/F:Chi2/F:NDF/I");
    }
  }
  for(Int_t i = 0; i < burstAsyms; i++){FitPolVar data; cycBurstAsymData.push_back(data);}
  for(Int_t i = 0; i < burstAsyms; i++){cyc->Branch(burstAsymTitles[i].Data(), &cycBurstAsymData[i], "mean/F:meanErr/F:Chi2/F:NDF/I");}
  for(Int_t i = 0; i < comboAsyms; i++){PolVar data; cycBurstComboData.push_back(data);}
  for(Int_t i = 0; i < comboAsyms; i++){cyc->Branch(comboAsymTitles[i].Data(), &cycBurstComboData[i], "mean/F:meanErr/F");}
}

void initRunTree(TTree *run){
  for(Int_t i = 0; i < nPolVars; i++){
    vector<Float_t> avg, err, burstAvg, burstErr;
    runPol0Avg.push_back(avg);
    runPol0Err.push_back(err);
    runBurstAvg.push_back(burstAvg);
    runBurstErr.push_back(burstErr);
    FitPolVar varAvg, varBurstAvg;
    runPol0.push_back(varAvg);
    runBurst.push_back(varBurstAvg);
  }
  run->Branch("snailNum", &snailNum, "snailNum/I");
  run->Branch("runNum", &runNum, "runNum/I");
  run->Branch("numCycles", &numRunCycles, "numCycles/I");
  run->Branch("numCyclesAcc", &numRunCyclesAcc, "numCyclesAcc/I");
  run->Branch("runTime", &runTime);
  for(Int_t i = 0; i < runMPSVars; i++){DataVar data; runMPSData.push_back(data);}
  for(Int_t i = 0; i < runMPSVars; i++){run->Branch(runMPSTitles[i].Data(), &runMPSData[i], "mean/F:meanErr/F:rms/F:rmsErr/F");}
  for(Int_t i = 0; i < runEpcVars; i++){StdVar data; runEpcData.push_back(data);}
  for(Int_t i = 0; i < runEpcVars; i++){run->Branch(runEpcTitles[i].Data(), &runEpcData[i], "mean/F");}
  for(Int_t i = 0; i < runBPMVars; i++){PolVar data; runBPMData.push_back(data);}
  for(Int_t i = 0; i < runBPMVars; i++){run->Branch(runBPMTitles[i].Data(), &runBPMData[i], "mean/F:meanErr/F");}
  run->Branch("LaserPolarization", &runLaserPol, "mean/F:meanErr/F");
  run->Branch("sign", &runSign, "sign/I");
  //run->Branch("Asym0", &runAsym0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //run->Branch("Asym0NGC", &runAsym0NGC, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //run->Branch("Asym0LasOn", &runAsym0On, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //run->Branch("Asym0LasOff", &runAsym0Off, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //run->Branch("Pol0", &runPol0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  for(Int_t i = 0; i < nPolVars; i++){run->Branch(Form("%s", polNames[i].Data()), &runPol0[i], "mean/F:meanErr/F:Chi2/F:NDF/I");}
  for(Int_t i = 0; i < nPolVars; i++){run->Branch(Form("Burst%s", polNames[i].Data()), &runBurst[i], "mean/F:meanErr/F:Chi2/F:NDF/I");}
}

void initSnailTree(TTree *snl){
  for(Int_t i = 0; i < nPolVars; i++){
    vector<Float_t> avg, err, burstAvg, burstErr;
    snlPol0Avg.push_back(avg);
    snlPol0Err.push_back(err);
    snlBurstAvg.push_back(burstAvg);
    snlBurstErr.push_back(burstErr);
    FitPolVar varAvg, varBurstAvg;
    snlPol0.push_back(varAvg);
    snlBurst.push_back(varBurstAvg);
  }
  snl->Branch("snailNum", &snailNum, "snailNum/I");
  snl->Branch("numRuns", &numSnlRuns, "numRuns/I");
  snl->Branch("numCycles", &numSnlCycles, "numCycles/I");
  snl->Branch("numCyclesAcc", &numSnlCyclesAcc, "numCyclesAcc/I");
  snl->Branch("snailTime", &snailTime);
  //snl->Branch("Asym0", &snlAsym0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Asym0NGC", &snlAsym0NGC, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Asym0LasOn", &snlAsym0On, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Asym0LasOff", &snlAsym0Off, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Pol0", &snlPol0, "mean/F:meanErr/F:Chi2/F:NDF/I");
  for(Int_t i = 0; i < nPolVars; i++){snl->Branch(polNames[i].Data(), &snlPol0[i], "mean/F:meanErr/F:Chi2/F:NDF/I");}
  //snl->Branch("Asym4", &snlAsym4, "mean/F:meanErr/F:Chi2/F:NDF/I");
  //snl->Branch("Pol4", &snlPol4, "mean/F:meanErr/F:Chi2/F:NDF/I");
  snl->Branch("qw1", &qw1, "qw1/F");
  snl->Branch("hw1", &hw1, "hw1/F");
  snl->Branch("qw2", &qw2, "qw2/F");
  snl->Branch("LaserPolarization", &laserPol, "LaserPolarization/F");
  snl->Branch("ihwp", &ihwp, "ihwp/F");
  snl->Branch("VWienAngle", &VWienAngle, "VWienAngle/F");
  snl->Branch("HWienAngle", &HWienAngle, "HWienAngle/F");
  snl->Branch("PhiFG", &PhiFG, "PhiFG/F");
  snl->Branch("sign", &snailSign, "sign/I");
  for(Int_t i = 0; i < nPolVars; i++){snl->Branch(Form("Burst%s", polNames[i].Data()), &snlBurst[i], "mean/F:meanErr/F:Chi2/F:NDF/I");}
}

/**
  Called before run and cycle loops
**/
void snailIterSet(Int_t base, Int_t snlInd){
  printf("Analyzing snail %i...\n", base + snlInd);
  snailNum = base + snlInd;
  numSnlCycles = 0;
  numSnlCyclesAcc = 0;
  numSnlRuns = 0;
  snailTime = 0;
  snailSign = 0;
  //snlAsym0Avg.clear(); snlAsym0Err.clear();
  //snlAsym0NGCAvg.clear(); snlAsym0NGCErr.clear();
  //snlAsym0OnAvg.clear(); snlAsym0OnErr.clear();
  //snlAsym0OffAvg.clear(); snlAsym0OffErr.clear();  
  //snlPol0Avg.clear(); snlPol0Err.clear();
  for(Int_t i = 0; i < nPolVars; i++){
    snlPol0Avg[i].clear(); snlPol0Err[i].clear();
    snlBurstAvg[i].clear(); snlBurstErr[i].clear();
  }
  //snlAsym4Avg.clear(); snlAsym4Err.clear();
  //snlPol4Err.clear(); snlPol4Err.clear();
  qw1 = 0; hw1 = 0; qw2 = 0; ihwp = 0; HWienAngle = 0; VWienAngle = 0; PhiFG = 0;
}

void snailIterAfterSet(Int_t base, Int_t snlInd){
  printf("Finalizing snail %i...\n", base + snlInd);
  TH1F *hPol0[nPolVars]; TF1 *fPol0[nPolVars];
  TH1F *hBurst[nPolVars]; TF1 *fBurst[nPolVars];
  for(Int_t i = 0; i < nPolVars; i++){
    hPol0[i] = new TH1F(Form("h%s_%i", polNames[i].Data(), base + snlInd), polTitles[i].Data(), numSnlCyclesAcc, 0, numSnlCyclesAcc);
    hBurst[i] = new TH1F(Form("hBurst%s_%i", polNames[i].Data(), base + snlInd), Form("Burst %s", polTitles[i].Data()), numSnlCyclesAcc, 0, numSnlCyclesAcc);
    fPol0[i] = new TF1(Form("f%s_%i", polNames[i].Data(), base + snlInd), "pol0");
    fBurst[i] = new TF1(Form("fBurst%s_%i", polNames[i].Data(), base + snlInd), "pol0");
  }
  //TH1F *hAsym0 = new TH1F(Form("hAsym0_%i", base + snlInd), "Asym0", numSnlCyclesAcc, 0, numSnlCyclesAcc);
  //TH1F *hAsym0NGC = new TH1F(Form("hAsym0NGC_%i", base + snlInd), "Asym0 NGC", numSnlCyclesAcc, 0, numSnlCyclesAcc);
  //TH1F *hAsym0On = new TH1F(Form("hAsym0On_%i", base + snlInd), "Asym0 On", numSnlCyclesAcc, 0, numSnlCyclesAcc);
  //TH1F *hAsym0Off = new TH1F(Form("hAsym0Off_%i", base + snlInd), "Asym0 Off", numSnlCyclesAcc, 0, numSnlCyclesAcc);
  //TH1F *hPol0 = new TH1F(Form("hPol0_%i", base + snlInd), "Pol0", numSnlCyclesAcc, 0, numSnlCyclesAcc);
  //TH1F *hAsym4 = new TH1F(Form("hAsym4_%i", base + snlInd), "Asym4", numSnlCycles, 0, numSnlCycles);
  //TH1F *hPol4 = new TH1F(Form("hPol4_%i", base + snlInd), "Pol4", numSnlCycles, 0, numSnlCycles);
  //TF1 *fAsym0 = new TF1(Form("fAsym0_%i", base + snlInd), "pol0");
  //TF1 *fAsym0NGC = new TF1(Form("fAsym0NGC_%i", base + snlInd), "pol0");
  //TF1 *fAsym0On = new TF1(Form("fAsym0On_%i", base + snlInd), "pol0");
  //TF1 *fAsym0Off = new TF1(Form("fAsym0Off_%i", base + snlInd), "pol0");
  //TF1 *fPol0 = new TF1(Form("fPol0_%i", base + snlInd), "pol0");
  //TF1 *fAsym4 = new TF1(Form("fAsym4_%i", base + snlInd), "pol0");
  //TF1 *fPol4 = new TF1(Form("fPol4_%i", base + snlInd), "pol0");
  for(Int_t i = 0; i < snlPol0Avg.size(); i++){
    for(Int_t j = 0; j < snlPol0Avg[i].size(); j++){
      hPol0[i]->SetBinContent(j + 1, snlPol0Avg[i][j]); hPol0[i]->SetBinError(j + 1, snlPol0Err[i][j]);
      hBurst[i]->SetBinContent(j + 1, snlBurstAvg[i][j]); hBurst[i]->SetBinError(j + 1, snlBurstErr[i][j]);
      //hAsym0->SetBinContent(i + 1, snlAsym0Avg[i]); hAsym0->SetBinError(i + 1, snlAsym0Err[i]);
      //hAsym0NGC->SetBinContent(i + 1, snlAsym0NGCAvg[i]); hAsym0NGC->SetBinError(i + 1, snlAsym0NGCErr[i]);
      //hAsym0On->SetBinContent(i + 1, snlAsym0OnAvg[i]); hAsym0On->SetBinError(i + 1, snlAsym0OnErr[i]);
      //hAsym0Off->SetBinContent(i + 1, snlAsym0OffAvg[i]); hAsym0Off->SetBinError(i + 1, snlAsym0OffErr[i]);
      //hPol0->SetBinContent(i + 1, snlPol0Avg[i]); hPol0->SetBinError(i + 1, snlPol0Err[i]);
      //hAsym4->SetBinContent(i + 1, snlAsym4Avg[i]); hAsym4->SetBinError(i + 1, snlAsym4Err[i]);
      //hPol4->SetBinContent(i + 1, snlPol4Avg[i]); hPol4->SetBinError(i + 1, snlPol4Err[i]);
      //printf("  Asyms, Cycle %i: %f +/- %f, %f +/- %f\n", i + 1, snlAsym0Avg[i], snlAsym0Err[i], snlAsym4Avg[i], snlAsym4Err[i]);
      //printf("  Pol, Cycle %i: %f +/- %f, %f +/- %f\n", i + 1, snlPol0Avg[i], snlPol0Err[i], snlPol4Avg[i], snlPol4Err[i]);
      //printf("  Asyms, Cycle %i: %f +/- %f\n", i + 1, snlAsym0Avg[i], snlAsym0Err[i]);
      //printf("  Pol, Cycle %i: %f +/- %f\n", i + 1, snlPol0Avg[i], snlPol0Err[i]);
    }
  }
  for(Int_t i = 0; i < nPolVars; i++){
    hPol0[i]->Fit(fPol0[i], "Q", "", 0, numSnlCyclesAcc);
    hBurst[i]->Fit(fBurst[i], "Q", "", 0, numSnlCyclesAcc);
    snlPol0[i].mean = fPol0[i]->GetParameter(0);
    snlPol0[i].meanErr = fPol0[i]->GetParError(0);
    snlPol0[i].Chi2 = fPol0[i]->GetChisquare();
    snlPol0[i].NDF = fPol0[i]->GetNDF();
    snlBurst[i].mean = fBurst[i]->GetParameter(0);
    snlBurst[i].meanErr = fBurst[i]->GetParError(0);
    snlBurst[i].Chi2 = fBurst[i]->GetChisquare();
    snlBurst[i].NDF = fBurst[i]->GetNDF();
  }
  //hAsym0->Fit(fAsym0, "Q", "", 0, numSnlCyclesAcc);
  //hAsym0NGC->Fit(fAsym0NGC, "Q", "", 0, numSnlCyclesAcc);
  //hAsym0On->Fit(fAsym0On, "Q", "", 0, numSnlCyclesAcc);
  //hAsym0Off->Fit(fAsym0Off, "Q", "", 0, numSnlCyclesAcc);
  //hPol0->Fit(fPol0, "Q", "", 0, numSnlCyclesAcc);
  //snlAsym0.mean = fAsym0->GetParameter(0); snlAsym0.meanErr = fAsym0->GetParError(0);
  //snlAsym0.Chi2 = fAsym0->GetChisquare(); snlAsym0.NDF = fAsym0->GetNDF();
  //snlAsym0NGC.mean = fAsym0NGC->GetParameter(0); snlAsym0NGC.meanErr = fAsym0NGC->GetParError(0);
  //snlAsym0NGC.Chi2 = fAsym0NGC->GetChisquare(); snlAsym0NGC.NDF = fAsym0NGC->GetNDF();
  //snlAsym0On.mean = fAsym0On->GetParameter(0); snlAsym0On.meanErr = fAsym0On->GetParError(0);
  //snlAsym0On.Chi2 = fAsym0On->GetChisquare(); snlAsym0On.NDF = fAsym0On->GetNDF();
  //snlAsym0Off.mean = fAsym0Off->GetParameter(0); snlAsym0Off.meanErr = fAsym0Off->GetParError(0);
  //snlAsym0Off.Chi2 = fAsym0Off->GetChisquare(); snlAsym0Off.NDF = fAsym0Off->GetNDF();
  //snlPol0.mean = fPol0->GetParameter(0); snlPol0.meanErr = fPol0->GetParError(0);
  //snlPol0.Chi2 = fPol0->GetChisquare(); snlPol0.NDF = fPol0->GetNDF();
  //hAsym4->Fit(fAsym4, "Q", "", 0, numSnlCycles);
  //hPol4->Fit(fPol4, "Q", "", 0, numSnlCycles);
  //snlAsym4.mean = fAsym4->GetParameter(0); snlAsym4.meanErr = fAsym4->GetParError(0);
  //snlAsym4.Chi2 = fAsym4->GetChisquare(); snlAsym4.NDF = fAsym4->GetNDF();
  //snlPol4.mean = fPol4->GetParameter(0); snlPol4.meanErr = fPol4->GetParError(0);
  //snlPol4.Chi2 = fPol4->GetChisquare(); snlPol4.NDF = fPol4->GetNDF();

  //printf("Sign Vars snail %i: IHWP %f, VWien: %f, HWien: %f, PhiFG: %f; Num Runs: %i\n", base + snlInd, ihwp, VWienAngle, HWienAngle, PhiFG, numSnlRuns);

  qw1 /= numSnlRuns; hw1 /= numSnlRuns; qw2 /= numSnlRuns; ihwp /= numSnlRuns;
  HWienAngle /= numSnlRuns; VWienAngle /= numSnlRuns; PhiFG /= numSnlRuns;
  calcSnailSign(base + snlInd);
  calcSnailLaserPol(snailNum);
}

/**
  Called before cycle loop
**/
void runIterSet(Int_t runNumber, Int_t numCycles, TFile *plotFile){
  printf("  Analyzing run %i...\n", runNumber);
  //runAsym0Avg.clear(); runAsym0Err.clear();
  //runAsym0NGCAvg.clear(); runAsym0NGCErr.clear();
  //runAsym0OnAvg.clear(); runAsym0OnErr.clear();
  //runAsym0OffAvg.clear(); runAsym0OffErr.clear();
  //runPol0Avg.clear(); runPol0Err.clear();
  for(Int_t i = 0; i < nPolVars; i++){
    runPol0Avg[i].clear(); runPol0Err[i].clear();
    runBurstAvg[i].clear(); runBurstErr[i].clear();
  }
  runNum = runNumber;
  numRunCycles = numCycles;
  numCyclesAcc = 0;
  numRunCyclesAcc = 0;
  numSnlCycles += numCycles;
  numSnlRuns += 1;
  runSign = 0;
  //runTime = 1.0*mpswise->GetEntries()/helicityFreq(runNumber);
  for(int i = 0; i < runMPSVars; i++){
    TString hName = Form("h%i_%s", runNumber, runMPSTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    runMPSData[i].mean = h->GetMean(); runMPSData[i].meanErr = h->GetMeanError();
    runMPSData[i].rms = h->GetRMS(); runMPSData[i].rmsErr = h->GetRMSError();
  }
  for(int i = 0; i < runEpcVars; i++){
    TString hName = Form("h%i_%s", runNumber, runEpcTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    runEpcData[i].mean = h->GetMean();
  }
  runTime = 1.0*((TH1F *)plotFile->Get(Form("h%i_mps", runNumber)))->GetEntries()/helicityFreq(runNumber);
  snailTime += runTime;
  Float_t ihwpAvg = ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "ihwp")))->GetMean();
  Int_t ihwpVar;
  if(ihwpAvg > 0.5) ihwpVar = 1;
  else ihwpVar = 0;
  qw1 += ((TH1F *)plotFile->Get(Form("h%i_qw1", runNumber)))->GetMean();
  hw1 += ((TH1F *)plotFile->Get(Form("h%i_hw1", runNumber)))->GetMean();
  qw2 += ((TH1F *)plotFile->Get(Form("h%i_qw2", runNumber)))->GetMean();
  ihwp += ihwpVar;
  Float_t hWienMean = ((TH1F *)plotFile->Get(Form("h%i_HWienAngle", runNumber)))->GetMean();
  Float_t vWienMean = ((TH1F *)plotFile->Get(Form("h%i_VWienAngle", runNumber)))->GetMean();
  Float_t solWienMean = ((TH1F *)plotFile->Get(Form("h%i_PhiFG", runNumber)))->GetMean();
  HWienAngle += hWienMean;
  VWienAngle += vWienMean;
  PhiFG += solWienMean;
  //printf("  Sign settings: IHWP: %i, VWien: %f, HWien: %f, Sol: %f\n", 
  //       ihwpVar, ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "VWienAngle")))->GetMean(), 
  //       ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "HWienAngle")))->GetMean(), ((TH1F *)plotFile->Get(Form("h%i_%s", runNumber, "PhiFG")))->GetMean());
  calcRunSign(runNumber, hWienMean, vWienMean, solWienMean, ihwpVar);
  calcRunLaserPol(snailNum);
}

void runIterAfterSet(Int_t runNumber){
  printf("  Finalizing run %i...\n", runNumber);
  TH1F *hPol0[nPolVars]; TF1 *fPol0[nPolVars];
  TH1F *hBurst[nPolVars]; TF1 *fBurst[nPolVars];
  numRunCyclesAcc = numCyclesAcc;
  //TH1F *hAsym0 = new TH1F(Form("hAsym0_%i", runNumber), "Asym0", numRunCyclesAcc, 0, numRunCyclesAcc);
  //TH1F *hAsym0NGC = new TH1F(Form("hAsym0NGC_%i", runNumber), "Asym0 NGC", numRunCyclesAcc, 0, numRunCyclesAcc);
  //TH1F *hAsym0On = new TH1F(Form("hAsym0On_%i", runNumber), "Asym0 On", numRunCyclesAcc, 0, numRunCyclesAcc);
  //TH1F *hAsym0Off = new TH1F(Form("hAsym0Off_%i", runNumber), "Asym0 Off", numRunCyclesAcc, 0, numRunCyclesAcc);
  //TH1F *hPol0 = new TH1F(Form("hPol0_%i", runNumber), "Pol0", numRunCyclesAcc, 0, numRunCyclesAcc);
  //TF1 *fAsym0 = new TF1(Form("fAsym0_%i", runNumber), "pol0");
  //TF1 *fAsym0NGC = new TF1(Form("fAsym0NGC_%i", runNumber), "pol0");
  //TF1 *fAsym0On = new TF1(Form("fAsym0On_%i", runNumber), "pol0");
  //TF1 *fAsym0Off = new TF1(Form("fAsym0Off_%i", runNumber), "pol0");
  //TF1 *fPol0 = new TF1(Form("fPol0_%i", runNumber), "pol0");
  for(Int_t i = 0; i < nPolVars; i++){
    hPol0[i] = new TH1F(Form("h%s_%i", polNames[i].Data(), runNumber), polTitles[i].Data(), numRunCyclesAcc, 0, numRunCyclesAcc);
    fPol0[i] = new TF1(Form("f%s_%i", polNames[i].Data(), runNumber), "pol0");
    hBurst[i] = new TH1F(Form("hBurst%s_%i", polNames[i].Data(), runNumber), Form("Burst %s", polTitles[i].Data()), numRunCyclesAcc, 0, numRunCyclesAcc);
    fBurst[i] = new TF1(Form("fBurst%s_%i", polNames[i].Data(), runNumber), "pol0");
  }
  //for(Int_t i = 0; i < runAsym0Avg.size(); i++){
  //  hAsym0->SetBinContent(i + 1, runAsym0Avg[i]); hAsym0->SetBinError(i + 1, runAsym0Err[i]);
  //  hAsym0NGC->SetBinContent(i + 1, runAsym0NGCAvg[i]); hAsym0NGC->SetBinError(i + 1, runAsym0NGCErr[i]);
  //  hAsym0On->SetBinContent(i + 1, runAsym0OnAvg[i]); hAsym0On->SetBinError(i + 1, runAsym0OnErr[i]);
  //  hAsym0Off->SetBinContent(i + 1, runAsym0OffAvg[i]); hAsym0Off->SetBinError(i + 1, runAsym0OffErr[i]);
  //  hPol0->SetBinContent(i + 1, runPol0Avg[i]); hPol0->SetBinError(i + 1, runPol0Err[i]);
  //}
  for(Int_t i = 0; i < nPolVars; i++){
    for(Int_t j = 0; j < runPol0Avg[i].size(); j++){
      hPol0[i]->SetBinContent(j + 1, runPol0Avg[i][j]); hPol0[i]->SetBinError(j + 1, runPol0Err[i][j]);
      hBurst[i]->SetBinContent(j + 1, runBurstAvg[i][j]); hBurst[i]->SetBinError(j + 1, runBurstErr[i][j]);
    }
  }
  //hAsym0->Fit(fAsym0, "Q", "", 0, numRunCyclesAcc);
  //hAsym0NGC->Fit(fAsym0NGC, "Q", "", 0, numRunCyclesAcc);
  //hAsym0On->Fit(fAsym0On, "Q", "", 0, numRunCyclesAcc);
  //hAsym0Off->Fit(fAsym0Off, "Q", "", 0, numRunCyclesAcc);
  //hPol0->Fit(fPol0, "Q", "", 0, numRunCyclesAcc);
  //runAsym0.mean = fAsym0->GetParameter(0); runAsym0.meanErr = fAsym0->GetParError(0);
  //runAsym0.Chi2 = fAsym0->GetChisquare(); runAsym0.NDF = fAsym0->GetNDF();
  //runAsym0NGC.mean = fAsym0NGC->GetParameter(0); runAsym0NGC.meanErr = fAsym0NGC->GetParError(0);
  //runAsym0NGC.Chi2 = fAsym0NGC->GetChisquare(); runAsym0NGC.NDF = fAsym0NGC->GetNDF();
  //runAsym0On.mean = fAsym0On->GetParameter(0); runAsym0On.meanErr = fAsym0On->GetParError(0);
  //runAsym0On.Chi2 = fAsym0On->GetChisquare(); runAsym0On.NDF = fAsym0On->GetNDF();
  //runAsym0Off.mean = fAsym0Off->GetParameter(0); runAsym0Off.meanErr = fAsym0Off->GetParError(0);
  //runAsym0Off.Chi2 = fAsym0Off->GetChisquare(); runAsym0Off.NDF = fAsym0Off->GetNDF();
  //runPol0.mean = fPol0->GetParameter(0); runPol0.meanErr = fPol0->GetParError(0);
  //runPol0.Chi2 = fPol0->GetChisquare(); runPol0.NDF = fPol0->GetNDF();
  for(Int_t i = 0; i < nPolVars; i++){
    hPol0[i]->Fit(fPol0[i], "Q", "", 0, numRunCyclesAcc);
    hBurst[i]->Fit(fBurst[i], "Q", "", 0, numRunCyclesAcc);
    runPol0[i].mean = fPol0[i]->GetParameter(0);
    runPol0[i].meanErr = fPol0[i]->GetParError(0);
    runPol0[i].Chi2 = fPol0[i]->GetChisquare();
    runPol0[i].NDF = fPol0[i]->GetNDF();
    runBurst[i].mean = fBurst[i]->GetParameter(0);
    runBurst[i].meanErr = fBurst[i]->GetParError(0);
    runBurst[i].Chi2 = fBurst[i]->GetChisquare();
    runBurst[i].NDF = fBurst[i]->GetNDF();
  }
  numSnlCyclesAcc += numRunCyclesAcc;
}

/**
  Main body of cycle loop
**/
void cycleIterSet(vector<vector<int>> cycles, Int_t cycInd, Int_t runNumber, TFile *plotFile){
  printf("    Analyzing cycle %i...\n", cycInd + 1);
  cycleNum = cycInd + 1;
  firstOffStartMPS = cycles[cycInd][0];
  firstOffEndMPS = cycles[cycInd][1];
  onStartMPS = cycles[cycInd][2];
  onEndMPS = cycles[cycInd][3];
  lastOffStartMPS = cycles[cycInd][4];
  lastOffEndMPS = cycles[cycInd][5];
  Int_t totMPS = (firstOffEndMPS - firstOffStartMPS) + (onEndMPS - onStartMPS) + (lastOffEndMPS - lastOffStartMPS);
  cycleTime = 1.0*totMPS/helicityFreq(runNumber);
  cycQW1 = runEpcData[2].mean;
  cycHW1 = runEpcData[3].mean;
  cycQW2 = runEpcData[4].mean;
  cycSign = runSign;
  cycIHWP = runEpcData[5].mean;
  cycVWien = runEpcData[7].mean;
  cycHWien = runEpcData[8].mean;
  cycSolWien = runEpcData[9].mean;
  cycQrtCalc.clear();
  for(Int_t i = 0; i < cycMPSVars; i++){
    TString hName = Form("h%i.%i_%s", runNumber, cycleNum, cycMPSTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    cycMPSData[i].mean = h->GetMean(); cycMPSData[i].meanErr = h->GetMeanError();
    cycMPSData[i].rms = h->GetRMS(); cycMPSData[i].rmsErr = h->GetRMSError();
  }
  for(Int_t i = 0; i < cycQrtVars; i++){
    TString hName = Form("h%i.%i_%s", runNumber, cycleNum, cycQrtTitles[i].Data());
    TH1F *h = (TH1F *)plotFile->Get(hName.Data());
    cycQrtData[i].mean = h->GetMean(); cycQrtData[i].meanErr = h->GetMeanError();
    cycQrtData[i].rms = h->GetRMS(); cycQrtData[i].rmsErr = h->GetRMSError();
    for(Int_t j = 0; j < cycPolVars; j++){
      if(cycQrtTitles[i].CompareTo(cycPolTitles[j]) == 0){
        PolVar pol; pol.mean = h->GetMean(); pol.meanErr = h->GetMeanError();
        cycQrtCalc.push_back(pol);
      }
    }
  }
  TString hAsym0Name = Form("h%i.%i_AsymAcc0LasOn", runNumber, cycleNum);
  TH1F *h0 = (TH1F *)plotFile->Get(hAsym0Name.Data());
  asym0LasOn.mean = h0->GetMean(); asym0LasOn.meanErr = h0->GetMeanError();
  TString hAsym0Off1Name = Form("h%i.%i_AsymAcc0LasOff1", runNumber, cycleNum);
  TH1F *h0Off1 = (TH1F *)plotFile->Get(hAsym0Off1Name.Data());
  asym0LasOff1.mean = h0Off1->GetMean(); asym0LasOff1.meanErr = h0Off1->GetMeanError();
  TString hAsym0Off2Name = Form("h%i.%i_AsymAcc0LasOff2", runNumber, cycleNum);
  TH1F *h0Off2 = (TH1F *)plotFile->Get(hAsym0Off2Name.Data());
  asym0LasOff2.mean = h0Off2->GetMean(); asym0LasOff2.meanErr = h0Off2->GetMeanError();
  TString hAsym4Name = Form("h%i.%i_AsymAcc4LasOn", runNumber, cycleNum);
  TH1F *h4 = (TH1F *)plotFile->Get(hAsym4Name.Data());
  asym4LasOn.mean = h4->GetMean(); asym4LasOn.meanErr = h4->GetMeanError();
  TString hAsym4Off1Name = Form("h%i.%i_AsymAcc4LasOff1", runNumber, cycleNum);
  TH1F *h4Off1 = (TH1F *)plotFile->Get(hAsym4Off1Name.Data());
  asym4LasOff1.mean = h4Off1->GetMean(); asym4LasOff1.meanErr = h4Off1->GetMeanError();
  TString hAsym4Off2Name = Form("h%i.%i_AsymAcc4LasOff2", runNumber, cycleNum);
  TH1F *h4Off2 = (TH1F *)plotFile->Get(hAsym4Off2Name.Data());
  asym4LasOff2.mean = h4Off2->GetMean(); asym4LasOff2.meanErr = h4Off2->GetMeanError();
  calcCyclePol(snailNum, runNumber);
  calcCycleBurstPol(snailNum, runNumber, cycleNum, plotFile);
  //calcCyclePedestals(plotFile);
}

/**
  Main function
**/
void buildGrandRootfile(Int_t prexOrCrex){
  TString expt("prex");
  if(prexOrCrex == 2){
    expt = "crex";
  }
  else if(prexOrCrex != 1 && prexOrCrex != 2){
    printf("Please enter correct parameter for experiment:\n");
    printf("    PREX == 1\n");
    printf("    CREX == 2\n\n");
    printf("Get it right this time, %s.\n", randInsult());
  }
  readKeysFile(expt);
  readCorrsFile(expt);

  vector<vector<int>> runList = productionRunList(prexOrCrex);
  TFile *out = new TFile(Form("%s/aggregates/%sGrandCompton.root", getenv("COMPMON_WEB"), expt.Data()), "RECREATE");
  TTree *cyc = new TTree("cyc", "cycle testing tree");
  TTree *run = new TTree("run", "run number testing tree");
  TTree *snl = new TTree("snl", "snail number testing tree");

  //initVars();
  initCycleTree(cyc);
  initRunTree(run);
  initSnailTree(snl);

  Int_t base = 100*(prexOrCrex - 1) + 1;
  for(Int_t s = 0; s < runList.size(); s++){
    //if(s < 29 || s > 39) continue;
    snailIterSet(base, s);
    for(Int_t r = 0; r < runList[s].size(); r++){
      TFile *plotFile = new TFile(Form("%s/Run%i_Plots.root", getenv("COMPMON_RUNPLOTS"), runList[s][r]), "READ");
      vector<vector<int>> cycles = findCycles(runList[s][r]);
      runIterSet(runList[s][r], (Int_t)cycles.size(), plotFile);
      for(Int_t c = 0; c < cycles.size(); c++){
        cycleIterSet(cycles, c, runList[s][r], plotFile);
        cyc->Fill();
      }
      runIterAfterSet(runList[s][r]);
      run->Fill();
      plotFile->Close();
    }
    snailIterAfterSet(base, s);
    snl->Fill();
  }

  out->cd();
  //printf("Alright to here!\n");
  cyc->Write();
  //printf("Still okay by me!\n");
  run->Write();
  //printf("No problems detected!\n");
  snl->Write();
  //printf("We're still going!\n");
  out->Close();
  //printf("Made it all the way!\n");
}
