#include <TChain.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <TH2F.h>

//const double kGlobalPed=3788.03;
const double kGlobalPed=3790.17;
const Float_t kADCToV = (5000.0/4095.);


const int kMinSample = 50;
const int kMaxSample = 200*1000;
const int kMinCheckSample = 25;
const int kMaxCheckSample = 275;
const double kMaxSampleDeviation = 3;

TChain *gChain = 0;
int gLoadedRun = 0;
TH1F *gHist = 0;
TH1F *gHist2 = 0;
TH1 *gFreq = 0;
TProfile *gProf = 0;
double gMaxStdDeviation = 3;
double gTimeConversion = 1; // leave it in terms of samples

int numSamples;
int test15;
const int maxSamples = 1000*0+220*0+1000;
const int maxPedSamples = 30;
float snap[maxSamples]; // Large enough to accomodate any reasonable window we toss at it
int snapClock;
float bcm;
int laserState;
int randomTime;
int mpsCoda;

void snapshot_ch1()
{
}

void config(int run)
{
  if(gLoadedRun != run) {
    if(gChain)
      delete gChain;

    gChain = new TChain("snapshots");
    //gChain->Add(TString::Format("../rootfiles/FadcCalo2016_%d.root",run));
    //gChain->Add(TString::Format("/data/cmuwork/rootfiles/compmon_%d.root",run));
    //gChain->Add(TString::Format("../rootfiles/compmon_%d.root",run));
    //gChain->Add(TString::Format("/home/compton/cornejo/rootfiles/compmon_%d.root",run));
    //gChain->Add(TString::Format("/raid5/cornejo/compton/dvcs/rootfiles/compmon_%d.root",run));
    //gChain->Add(TString::Format("/raid5/cornejo/compton/dvcs/rootfiles/compmon_%d.root",run));
    gChain->Add(TString::Format("/data/cmuwork/rootfiles/prex/compmon_ch1_%d.root",run));

    gChain->SetBranchStatus("*",0);
    gChain->SetBranchStatus("numSamples",1);
    gChain->SetBranchStatus("test15",1);
    gChain->SetBranchStatus("randomTime",1);
    gChain->SetBranchStatus("mpsCoda",1);
    gChain->SetBranchStatus("snap",1);
    gChain->SetBranchStatus("snapClock",1);
    //gChain->SetBranchStatus("bcm",1);
    gChain->SetBranchStatus("laserState",1);
    gChain->SetBranchAddress("numSamples",&numSamples);
    gChain->SetBranchAddress("test15",&test15);
    gChain->SetBranchAddress("randomTime",&randomTime);
    gChain->SetBranchAddress("snap",&snap);
    gChain->SetBranchAddress("snapClock",&snapClock);
    gChain->SetBranchAddress("bcm",&bcm);
    gChain->SetBranchAddress("laserState",&laserState);
    gChain->SetBranchAddress("mpsCoda",&mpsCoda);
  }
  if(gHist)
    delete gHist;
  if(gProf)
    delete gProf;

  gLoadedRun = run;
}

void averageSnapshot(int run, bool randomOnly = false)
{
  config(run);
  int entries = gChain->GetEntries();
  gProf = new TProfile("gProf","",//
      maxSamples,0,maxSamples*gTimeConversion);
  gProf->SetStats(false);
  float data[maxSamples] = {0.0};
  float count[maxSamples] = {0.0};

  int goodEntries = 0;
  bool good = false;
  int count0 = 0;;
  int count4095 = 0;
  for(int entry = 0; entry < entries; entry++) {
    gChain->GetEntry(entry);

    if(randomOnly&&!randomTime)
      continue;
    if(!randomOnly&&randomTime)
      continue;

    double mean = 0;
    double counts = 0;
    double stdDev = 0;
    count0 = count4095 = 0;
    for(int s = 0; s < maxSamples; s++) {
      mean += snap[s];
      counts++;
      if(snap[s]==0)
        count0++;
      if(snap[s]==4095)
        count4095++;
    }
    mean /= counts;
    if(mean < 3791) {
      //continue;
    }
    if(count0>5&&count4095>5) {
      std::cerr << "Bad snapshot " << entry << "!" << std::endl;
      continue;
    }

    good = true;

    // Compute the standard deviation
    if(randomOnly){
      for(int s = 0; s < maxSamples; s++) {
        stdDev += TMath::Power(snap[s]-mean,2);
      }
      stdDev = TMath::Sqrt(stdDev/(counts-1));
      for(int s = 0; s < maxSamples; s++) {
        if(TMath::Abs(snap[s]-mean)>gMaxStdDeviation*stdDev) {
          good = false;
        }
      }
    }

    if(!good&&false)
      continue;

    //std::cerr << "Entry=" << entry << std::endl;
    for(int s = 0; s < numSamples && s < maxSamples; s++) {
      data[s] += snap[s];
      count[s]++;
      //gProf->Fill(s,data[s],1);
      gProf->Fill(s*gTimeConversion,snap[s],1);
      if(randomOnly&&snap[s]<3800) {
        std::cerr << "i=" << entry << " snap[" << s << "]= " << snap[s] << std::endl;
      }
    }
    goodEntries++;
    //std::cout << "goodEntries=" << goodEntries << std::endl;
  }
  /*
  for(int s = 0; s < maxSamples; s++) {
    data[s] /= count[s];
    gHist->SetBinContent(s+1,data[s]);
  }*/
  gProf->SetTitle(TString::Format("Run %d Avg %s (%d entries)",run,
        (randomOnly?" (random) ":""),
        goodEntries));
  //gProf->SetMinimum(3840);

  TCanvas *c = new TCanvas("c","c",800,400);
  gProf->Draw();
  c->SaveAs(TString::Format("results/average_snapshot_%d.png",gLoadedRun));
  //gFreq = gHist->FFT(0,"RE");
  //gFreq->Draw();
  // Save out the average snapshot
  double ped = 0.;
  for(int bin = 10; bin < 40; bin++) {
    ped += gProf->GetBinContent(bin);
  }
  ped /= (40.-10.);
  int num_bins = gProf->GetNbinsX();
  fstream outAvg;
  outAvg.open(TString::Format("results/average_snapshot_%d.dat",gLoadedRun),
      std::ios::out);
  outAvg << ped << std::endl;
  for(int bin = 0; bin < num_bins; bin++) {
    outAvg << bin*5 << " " << gProf->GetBinContent(bin+1) << " "
      << gProf->GetBinError(bin+1) << std::endl;
  }
  outAvg.close();
}

int gStepRandomSnapshotNextEntry = 0;
TCanvas *canvRandom = 0;
/*
void stepRandomSnapshot(int run, int startEntry)
{
  config(run);
  if(startEntry<0)
    startEntry = gStepRandomSnapshotNextEntry;

  if(gHist2)
    delete gHist2;

  if(!canvRandom) {
    canvRandom = new TCanvas("canvRandom","canvRandom",600,500*2);
    canvRandom->Divide(1,2);
  }

  for(int entry = startEntry; entry < gChain->GetEntries(); entry++) {
    gChain->GetEntry(entry);
    if(randomTime>0) {
      gHist = new TH1F("gHist",TString::Format("Run %d Snapshot %d (mps=%d)"
            ,run,entry,mpsCoda),
          numSamples,0,numSamples);
      gHist2 = new TH1F("gHist2","",3860-3830,3830,3860);
      gHist2->SetStats(false);
      gHist->SetStats(false);

      for(int s = 0; s < numSamples && s < maxSamples; s++) {
        gHist->SetBinContent(s+1,snap[s]);
        gHist2->Fill(snap[s]);
      }
      canvRandom->cd(1);
      gHist->Draw();
      canvRandom->cd(2);
      gHist2->Draw();
      gStepRandomSnapshotNextEntry=entry+1;

      for(Int_t bin = 1; bin <= gHist2->GetNbinsX(); bin++) {
        std::cout << "Bin: " << bin << ", lowEdge: " << 
          gHist2->GetBinLowEdge(bin) << ", content: "
         << gHist2->GetBinContent(bin) << std::endl;
      }
      return;
    }
  }
  std::cout << "Apparently nothing got plotted" << std::endl;
}
*/
void stepRandomSnapshot(int run, int startEntry = -1)
{
  config(run);
  if(startEntry<0)
    startEntry = gStepRandomSnapshotNextEntry;

  for(int entry = startEntry; entry < gChain->GetEntries(); entry++) {

    int maxY = 0;
    int minY = 1e5; 
    int minX = -100;
    int sumsum = 0;
    double avg = 0;
    double lped = 0;
    gChain->GetEntry(entry);
    if(numSamples>maxSamples)
      numSamples = maxSamples;
    if(randomTime==1) {
      gHist = new TH1F("gHist",TString::Format("Run %d Snapshot %d (mps=%d,entry=%d,clock=%d,#S:%d)"
            ,run,entry,mpsCoda,entry,snapClock,test15),
          numSamples,0,numSamples);
      gHist->SetStats(false);
      lped = 0.;
      for(int s = 0; s < 30; s++) {
        lped += snap[s];
      }
      lped /= 30.;
      for(int s = 0; s < numSamples && s < maxSamples; s++) {
        gHist->SetBinContent(s+1,snap[s]);
        //sumsum += snap[s];
        sumsum += lped-snap[s];
        avg += snap[s];
        //lped+= kGlobalPed;
        if ( minY > snap[s] ) {
           minY = snap[s];
           minX = s;
        }
        if ( maxY < snap[s] ) {
           maxY = snap[s];
        }
      }
      avg /= Double_t(numSamples);
      //sumsum=lped-sumsum;
      gHist->Draw();
      gStepRandomSnapshotNextEntry=entry+1;
      //std::cout << "Plotted!! MinY: " << minY << ", minX: " << minX << ", maxY: " << maxY << ", sum: " <<(835381- sumsum) << std::endl;
      std::cout << "Plotted!! MinY: " << minY << ", minX: " << minX << ", maxY: " << maxY << ", sum: " <<(sumsum)  << ", avg: " << avg << std::endl;
      return;
    }
  }
 std::cout << "Apparently nothing got plotted" << std::endl;
}


int gStepTriggSnapshotNextEntry = 1;
void stepTrigSnapshot(int run, int startEntry = -1)
{
  config(run);
  if(startEntry<0)
    startEntry = gStepTriggSnapshotNextEntry;

  for(int entry = startEntry; entry < gChain->GetEntries(); entry++) {

    int maxY = 0;
    int minY = 1e5; 
    int minX = -100;
    int sumsum = 0;
    double avg = 0;
    double lped = 0;
    gChain->GetEntry(entry);
    if(numSamples>maxSamples)
      numSamples = maxSamples;
    if(randomTime==0) {
      gHist = new TH1F("gHist",TString::Format("Run %d Snapshot %d (mps=%d,entry=%d,clock=%d,#S:%d)"
            ,run,entry,mpsCoda,entry,snapClock,test15),
          numSamples,0,numSamples);
      gHist->SetStats(false);
      lped = 0.;
      for(int s = 0; s < 30; s++) {
        lped += snap[s];
      }
      lped /= 30.;
      for(int s = 0; s < numSamples && s < maxSamples; s++) {
        gHist->SetBinContent(s+1,snap[s]);
        //sumsum += snap[s];
        sumsum += lped-snap[s];
        avg += snap[s];
        //lped+= kGlobalPed;
        if ( minY > snap[s] ) {
           minY = snap[s];
           minX = s;
        }
        if ( maxY < snap[s] ) {
           maxY = snap[s];
        }
      }
      avg /= Double_t(numSamples);
      //sumsum=lped-sumsum;
      gHist->Draw();
      gStepTriggSnapshotNextEntry=entry+1;
      //std::cout << "Plotted!! MinY: " << minY << ", minX: " << minX << ", maxY: " << maxY << ", sum: " <<(835381- sumsum) << std::endl;
      std::cout << "Plotted!! MinY: " << minY << ", minX: " << minX << ", maxY: " << maxY << ", sum: " <<(sumsum)  << ", avg: " << avg << std::endl;
      return;
    }
  }
 std::cout << "Apparently nothing got plotted" << std::endl;
}

void stepTrigSnapshotSum(int run, double sum_gt, double sum_lt, int startEntry = -1, int randomOnly = 0)
{
  config(run);
  if(startEntry<0)
    startEntry = gStepTriggSnapshotNextEntry;

  for(int entry = startEntry; entry < gChain->GetEntries(); entry++) {

    int maxY = 0;
    int minY = 1e5; 
    int minX = -100;
    int sumsum = 0;
    double lped=0.0;
    gChain->GetEntry(entry);
    //if(gHist)
    //  delete gHist;

    if(numSamples>maxSamples)
      numSamples = maxSamples;
    if(randomTime==randomOnly) {
      lped = 0;
      for(int s = 0; s < 30; s++) {
        lped += snap[s];
      }
      lped /= 30.;
      for(int s = 0; s < numSamples && s < maxSamples; s++) {
        sumsum += snap[s];
        //lped += kGlobalPed;
        if ( minY > snap[s] ) {
           minY = snap[s];
           minX = s;
        }
        if ( maxY < snap[s] ) {
           maxY = snap[s];
        }
      }
      sumsum = lped*numSamples - sumsum;
      if(sumsum >= sum_gt && sumsum <= sum_lt) {
        gHist = new TH1F("gHist",TString::Format("Run %d Snapshot %d (mps=%d,entry=%d,clock=%d,#S:%d)"
              ,run,entry,mpsCoda,entry,snapClock,test15),
            numSamples,0,numSamples);
        gHist->SetStats(false);
        for(int s = 0; s < numSamples && s < maxSamples; s++) {
          gHist->SetBinContent(s+1,snap[s]);
        }
      	gHist->Draw();
      	gStepTriggSnapshotNextEntry=entry+1;
      	std::cout << "Plotted!! MinY: " << minY << ", minX: " << minX << ", maxY: " << maxY << ", sum: " <<(sumsum) << ", V=" << kADCToV*(kGlobalPed-minY) << " mV" << ", randomTime=" << randomTime << std::endl;
      	return;
       }
    }
  }
 std::cout << "Apparently nothing got plotted" << std::endl;
}



void stepTrigCal(int run, int startEntry = -1)
{
  config(run);
  if(startEntry<0)
    startEntry = gStepTriggSnapshotNextEntry;

  for(int entry = startEntry; entry < gChain->GetEntries(); entry++) {

    int maxY = 0;
    int minY = 1e5; 
    int minX = -100;
    int sumsum = 0;
    gChain->GetEntry(entry);
    if(numSamples>maxSamples)
      numSamples = maxSamples;
    if(randomTime==0) {
      gHist = new TH1F("gHist",TString::Format("Run %d Snapshot %d (mps=%d,entry=%d,clock=%d,#S:%d);Sample Number;Signal (mV)"
            ,run,entry,mpsCoda,entry,snapClock,test15),
          numSamples,0,numSamples);
      gHist->SetStats(false);
      float pedSum = 0;
      int nPedSamps = 0;
      for(int s = 0; s < numSamples && s < maxPedSamples; s++) {
        pedSum += snap[s];
        nPedSamps++;
      }
      if(nPedSamps>0) {
        pedSum /= float(nPedSamps);
      }
      float nSamps = 0; 
      for(int s = 0; s < numSamples && s < maxSamples; s++) {
        gHist->SetBinContent(s+1,(snap[s]-pedSum)*kADCToV);
        sumsum += snap[s];
        nSamps++;
        if ( minY > snap[s] ) {
           minY = snap[s];
           minX = s;
        }
        if ( maxY < snap[s] ) {
           maxY = snap[s];
        }
      }
      gHist->Draw();
      gStepTriggSnapshotNextEntry=entry+1;
      Float_t voltage = kADCToV*(minY-pedSum);
      sumsum = 835381-sumsum;
      std::cout << "Plotted!! MinY: " << minY << ", minX: " << minX << ", maxY: " << maxY << ", v=" << voltage << "(" << voltage/0.75 << "), sum: " << sumsum << "(" << sumsum/0.75 << ")"  << std::endl;
      return;
    }
  }
 std::cout << "Apparently nothing got plotted" << std::endl;
}


void displaySnapshot(int run, int entry)
{
  config(run);
  gChain->GetEntry(entry);
  if(numSamples>maxSamples)
    numSamples = maxSamples;
  gHist = new TH1F("gHist",TString::Format("Run %d Snapshot %d",run,entry),
      numSamples,0,numSamples);
  gHist->SetStats(false);
  double tmpMean = 0;
  double tmpStdDev = 0;
  for(int s = 0; s < numSamples && s < maxSamples; s++) {
    gHist->SetBinContent(s+1,snap[s]);
    tmpMean += snap[s];
    //std::cout << "snap[" << s << "] = " << snap[s] << std::endl;
  }
  tmpMean /= double(numSamples);

  for(int s = 0; s < numSamples && s < maxSamples; s++) {
    tmpStdDev += TMath::Power(tmpMean-snap[s],2);
  }
  tmpStdDev = TMath::Sqrt(tmpStdDev/double(numSamples-1));
  std::cout << "Full window mean = " << tmpMean << " stdDev = " << tmpStdDev 
    << " errorMean = " << tmpStdDev / TMath::Sqrt(numSamples) << std::endl;
  TCanvas *c = new TCanvas("c","c",600,400);
  gHist->Draw();
  c->SaveAs("results/snapshot.png");
}

void findBad(int run, int startEntry)
{
  config(run);
  int count0;
  int count4095;
  for(int entry = startEntry; entry < gChain->GetEntries(); entry++) {
    gChain->GetEntry(entry);
    count0 = count4095 = 0;
    for(int s = 200; s < 300; s++) {
      if(snap[s]==0)
        count0++;
      if(snap[s]==4095)
        count4095++;
    }
    if(count0>5&&count4095>5&&!randomTime) {
      std::cout << "Encountered bad snapshot! Entry=" << entry
       << " mpsCoda=" << mpsCoda << " randomTime=" << randomTime  << std::endl;
      displaySnapshot(run,entry);
      return;
    }
  }
}

void pedestalVsCurrent(int run, int laser = 3)
{
  config(run);
  double mean;
  double stdDev;
  double n;
  bool isGood;
  TH2F *hPedVsBCM = new TH2F("hPedVsBCM","Pedestal vs Current;"
      "Current (#muA);Pedestal (rau)",
      1000,-2.,20.,1000,3710,3970);
  hPedVsBCM->SetStats(false);

  TCanvas *c = new TCanvas("c","c",800,400);
  for(int entry = 25; entry < gChain->GetEntries(); entry++) {
    // Reset temp variables
    mean = stdDev = n = 0.0;
    isGood = true;

    // Get entry from tree
    gChain->GetEntry(entry);

    // Skip triggered snapshots
    if(!randomTime)
      continue;

    // Skip bad snapshots!
    if(snapClock>6600200)
      continue;

    // If laser state selected, skip others
    if(laser>=0&&laserState!=laser)
      continue;

    for(int s = kMinSample; s < kMaxSample; s++) {
      if(snap[s]<100) {
        isGood = false;
      }
      mean += snap[s];
      n++;
    }
    mean /= n; // compute mean
    if(mean < 3700) {
      //std::cerr << "Bad snapshot: " << entry << std::endl;
      continue;
    }
    // Now compute the standard deviation
    for( int s = kMinSample; s < kMaxSample; s++) {
      stdDev += TMath::Power(mean-snap[s],2);
    }
    stdDev = TMath::Sqrt(stdDev/(n-1));

    // Now see if there is anything in there that deviates more
    // than 3 times the standard deviation, if so, it indicates
    // a possible photon in that region
    for(int s = kMinCheckSample; s < kMaxCheckSample && isGood; s++) {
      if(TMath::Abs(snap[s]-mean)>kMaxSampleDeviation*stdDev) {
        isGood = false;
      }
    }
    if(isGood)
      hPedVsBCM->Fill(bcm,mean);
  }
  gPad->SetGrid(true,true);
  hPedVsBCM->Draw("COLZ");
  hPedVsBCM->Fit("pol1");
  c->SaveAs(TString::Format("results/snapshot_ped_vs_bcm_%d.png",run));
}

void pedestalVsTime(int run, int laser = 3)
{
  config(run);
  double mean;
  double stdDev;
  double n;
  bool isGood;
  TH2F *hPedVsTime = new TH2F("hPedVsTime","Pedestal vs MPS;"
      "MPS Entry;Pedestal (rau)",
      1000,0,gChain->GetEntries()+1,1000,3820,3850);
  hPedVsTime->SetStats(false);

  TCanvas *c = new TCanvas("c","c",800,400);
  double sumCoda;
  double sumMean;
  double count;
  for(int entry = 25; entry < gChain->GetEntries(); entry++) {
    // Reset temp variables
    mean = stdDev = n = 0.0;
    isGood = true;

    if(count>10) {
      hPedVsTime->Fill(sumCoda/count,sumMean/count);
      sumCoda = 0.;
      sumMean = 0.;
      count = 0.;
    }

    // Get entry from tree
    gChain->GetEntry(entry);

    // Skip triggered snapshots
    if(!randomTime)
      continue;

    // Skip bad snapshots!
    if(snapClock>6600200)
      continue;

    // If laser state selected, skip others
    if(laser>=0&&laserState!=laser)
      continue;

    for(int s = kMinSample; s < kMaxSample; s++) {
      if(snap[s]<100) {
        isGood = false;
      }
      mean += snap[s];
      n++;
    }
    mean /= n; // compute mean
    if(mean < 3700) {
      //std::cerr << "Bad snapshot: " << entry << std::endl;
      continue;
    }
    // Now compute the standard deviation
    for( int s = kMinSample; s < kMaxSample; s++) {
      stdDev += TMath::Power(mean-snap[s],2);
    }
    stdDev = TMath::Sqrt(stdDev/(n-1));

    // Now see if there is anything in there that deviates more
    // than 3 times the standard deviation, if so, it indicates
    // a possible photon in that region
    for(int s = kMinCheckSample; s < kMaxCheckSample && isGood; s++) {
      if(TMath::Abs(snap[s]-mean)>kMaxSampleDeviation*stdDev) {
        isGood = false;
      }
    }
    if(isGood) {
      sumCoda += mpsCoda;
      sumMean += mean;
      count++;
    }
  }
  gPad->SetGrid(true,true);
  hPedVsTime->Draw();
  c->SaveAs(TString::Format("results/snapshot_ped_vs_time_%d.png",run));
}

void plotPeakValue(int run, int laser = 1)
{
  config(run);
  TCanvas *c = new TCanvas("c","c",800,800);

  double peakVal;
  int    peakPos;
  int minPeakPos = 0+60;
  int maxPeakPos = 300*0+90;
  int minPeak = 3100*0;
  int maxPeak = 3900;
  TH1F *hPeakValue = new TH1F("hPeakValue","",
      maxPeak-minPeak,minPeak,maxPeak);

  //TProfile *hPeakValueVsPos = new TProfile("hPeakValueVsPos","",maxPeakPos-minPeakPos,minPeakPos,maxPeakPos);
  //hPeakValueVsPos->SetMinimum(minPeak);
  //hPeakValueVsPos->SetMarkerStyle(20);
  //hPeakValueVsPos->SetMarkerColor(kRed+1);

  for(int entry = 0; entry < gChain->GetEntries(); entry++) {
    // Get entry from tree
    gChain->GetEntry(entry);

    // Skip random snapshots
    if(randomTime)
      continue;

    // Skip bad snapshots!
    if(snapClock>6600200*0+830e3)
      continue;

    // If laser state selected, skip others
    if(laser>=0&&laserState!=laser)
      continue;

    peakVal = 1e6;
    peakPos = -1;
    for(int s = 0; s < 300; s++) {
      if(snap[s]<3800 && snap[s]<peakVal) {
        peakVal = snap[s];
        peakPos = s;
      }
    }

    if(peakPos>=0) {
	   hPeakValue->Fill(peakVal);
    }
    if(peakPos>78)
      std::cout << "PeakPos: " << peakPos << " for entry " << entry << std::endl;
  }
  gPad->SetGrid(true,true);
  c->cd(1);
  hPeakValue->Draw();
  //hPeakValueVsPos->Fit("pol1");
  //hPeakPos[0]->Fit("gaus","","",50.,100.);
  //c->SaveAs(TString::Format("results/snapshot_peak_position_vs_time_%d.png",run));

}

void plotPeakValueVsPos(int run, int laser = 1)
{
  config(run);
  TCanvas *c = new TCanvas("c","c",800,800);

  double peakVal;
  int    peakPos;
  int minPeakPos = 0+60;
  int maxPeakPos = 300*0+90;
  int minPeak = 3100;
  int maxPeak = 3900;
  TH2F *hPeakValueVsPos = new TH2F("hPeakValueVsPos","",maxPeakPos-minPeakPos,minPeakPos,maxPeakPos,
      maxPeak-minPeak,minPeak,maxPeak);

  //TProfile *hPeakValueVsPos = new TProfile("hPeakValueVsPos","",maxPeakPos-minPeakPos,minPeakPos,maxPeakPos);
  //hPeakValueVsPos->SetMinimum(minPeak);
  //hPeakValueVsPos->SetMarkerStyle(20);
  //hPeakValueVsPos->SetMarkerColor(kRed+1);

  for(int entry = 0; entry < gChain->GetEntries(); entry++) {
    // Get entry from tree
    gChain->GetEntry(entry);

    // Skip random snapshots
    if(randomTime)
      continue;

    // Skip bad snapshots!
    if(snapClock>6600200)
      continue;

    // If laser state selected, skip others
    if(laser>=0&&laserState!=laser)
      continue;

    peakVal = 1e6;
    peakPos = -1;
    for(int s = 0; s < 300*0+220; s++) {
      if(snap[s]<3800 && snap[s]<peakVal) {
        peakVal = snap[s];
        peakPos = s;
      }
    }

    if(peakPos>=0) {
	   hPeakValueVsPos->Fill(peakPos,peakVal);
    }
    if(peakPos>78)
      std::cout << "PeakPos: " << peakPos << " for entry " << entry << std::endl;
  }
  gPad->SetGrid(true,true);
  c->cd(1);
  hPeakValueVsPos->Draw();
  //hPeakValueVsPos->Fit("pol1");
  //hPeakPos[0]->Fit("gaus","","",50.,100.);
  //c->SaveAs(TString::Format("results/snapshot_peak_position_vs_time_%d.png",run));

}

void plotPeakPosVsTime(int run, int laser = 1)
{
  config(run);
  TCanvas *c = new TCanvas("c","c",800,800);
  c->Divide(1,2);
  c->cd(1);

  double peakVal;
  int    peakPos;
  int minPeakPos = -1*0+50;
  int maxPeakPos = 301*0+100;
  TH2F *hPeakPosVsTime = new TH2F("hPeakPosVsTime","Peak Position vs MPS;"
      "MPS Entry;Peak Position at Sample No.)",
      300,0,300,maxPeakPos-minPeakPos,minPeakPos,maxPeakPos);
  hPeakPosVsTime->SetStats(false);
  TH1F *hPeakPos[2];
  hPeakPos[0] = new TH1F("hPeakPos0","Peak Position",maxPeakPos-minPeakPos,minPeakPos,maxPeakPos);
  hPeakPos[0]->SetStats(false);
  hPeakPos[0]->SetLineColor(kBlue+1);

  for(int entry = 25; entry < gChain->GetEntries(); entry++) {
    // Get entry from tree
    gChain->GetEntry(entry);

    // Skip random snapshots
    if(randomTime)
      continue;

    // Skip bad snapshots!
    if(snapClock>6600200)
      continue;

    // If laser state selected, skip others
    if(laser>=0&&laserState!=laser)
      continue;

    peakVal = 1e6;
    peakPos = -1;
    for(int s = 10; s < 290; s++) {
      if(snap[s]<3400 && snap[s]<peakVal) {
        peakVal = snap[s];
        peakPos = s;
      }
    }

    if(peakPos>=0) {
      hPeakPosVsTime->Fill(mpsCoda,peakPos);
      if(snapClock > 0 && snapClock < (6601042-700)/ 2.) {
        hPeakPos[0]->Fill(peakPos);
      } else if ( snapClock > (6601042-700)/2. && snapClock < 6601042-700 ) {
        hPeakPos[1]->Fill(peakPos);
      }
      if(peakPos > 100) {
        std::cout << "Peak pos = " << peakPos << " for entry = " << entry << std::endl;
      }
    }
  }
  gPad->SetGrid(true,true);
  c->cd(1);
  hPeakPosVsTime->Draw();
  c->cd(2);
  hPeakPos[0]->Draw();
  hPeakPos[1]->Draw("SAME");
  hPeakPos[0]->Fit("gaus","","",50.,100.);
  hPeakPos[1]->Fit("gaus","","",50.,100.);
  c->SaveAs(TString::Format("results/snapshot_peak_position_vs_time_%d.png",run));

}

void plotPeakValueVsSum(int run)
{
  config(run);
  //TCanvas *c = new TCanvas("c","c",800,800);
  //TCanvas *c2 = new TCanvas("c2","c2",800,800);
  TCanvas *c = new TCanvas("c","c",600,600);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  TCanvas *c4 = new TCanvas("c4","c4",600,600);
  TCanvas *c5 = new TCanvas("c5","c5",600,600);

  double peakVal;
  double sum;
  double width = 0;
  double wbefore = 0;
  double wafter = 0;
  int peakPos = 0;
  int minSum  = 0;
  int maxSum  = 100000;
  float minWidth  = 0;
  float maxWidth  = 40.0*5;
  int minPeak = 0;
  int maxPeak = 4095;
  TH1F *hWidthFit = (TH1F*)gDirectory->Get("hWidthFit");
  if(!hWidthFit) {
    gChain->GetEntry(0); // Find first entry to set limits
    if(numSamples>maxSamples) {
      numSamples = maxSamples;
    }
    hWidthFit = new TH1F("hWidthFit","",numSamples,0,numSamples*5.);
  }
  TH2F *hPeakValueVsSum = new TH2F("hPeakValueVsSum","",
    1000,minSum,maxSum,
    maxPeak-minPeak,minPeak,maxPeak);

  TH2F *hPeakValueVsWidth = new TH2F("hPeakValueVsWidth","",
    100,minWidth,maxWidth,
    maxPeak-minPeak,minPeak,maxPeak);

  TH2F *hWidthValueVsSum = new TH2F("hWidthValueVsSum","",
    1000,minSum,maxSum,
    maxWidth-minWidth,minWidth,maxWidth);

  TH1F *hWidth = new TH1F("hWidth","hWidth",100.,minWidth,maxWidth);
  TH1F *hPeak = new TH1F("hPeak","hPeak",100.,minPeak,maxPeak);

  // A quick model for the pulse to see if it gets a better fit
  //TF1 *f1 = new TF1("f1Pulse1","[0]*TMath::Max(0.,x/([1]**2))*TMath::Exp(-x/[1])");


  //TProfile *hPeakValueVsSum = new TProfile("hPeakValueVsSum","",maxPeakPos-minPeakPos,minPeakPos,maxPeakPos);
  //
  //hPeakValueVsSum->SetMinimum(minPeak);
  //hPeakValueVsSum->SetMarkerStyle(20);
  //hPeakValueVsSum->SetMarkerColor(kRed+1);

  int foundPeak = false;
  bool foundWidth = false;
  bool foundwafter = false;
  bool foundwbefore = false;
  float pedsum = 0;
  TFitResultPtr res;
  for(int entry = 0; entry < gChain->GetEntries(); entry++) {
    // Get entry from tree
    gChain->GetEntry(entry);
    hWidthFit->Reset();

    // Skip random snapshots
    if(randomTime)
      continue;

    // Skip bad snapshots!
    if(snapClock>6600200)
      continue;

    //if(bcm<0.1)
    //  continue;

    peakVal = 1e6;
    sum = 0;
    foundPeak = false;
    foundWidth = false;
    foundwafter = false;
    foundwbefore = false;
    width = wafter = wbefore = -4095;
    pedsum = 0.0;
    for(int s = 0; s < 30; s++) {
        pedsum += snap[s];
    }
    pedsum /= 30.0;
    bool notSat = true;
    for(int s = 0; s < 300*0+220&&notSat; s++) {
      //sum += kGlobalPed-snap[s];
      sum += pedsum-snap[s];
      if(snap[s]<peakVal) {
        peakVal = snap[s];
        foundPeak = true;
        peakPos =  s;
      }
      if(snap[s]<1)
        notSat = false;
      hWidthFit->SetBinContent(1+s,pedsum-snap[s]);
      // Now fit the snapshot
    }
    if(foundPeak && peakVal < 3700 && notSat) {
      res = hWidthFit->Fit("landau","QS","",(peakPos-40.)*5,(peakPos+1)*5);
      Float_t dx = 0.1;
      //hWidthFit->GetListOfFunctions()->Print();
      //Float_t pvh = 0.5*(kGlobalPed-peakVal);
      Float_t pvh = 0.5*(pedsum-peakVal);
      for(Float_t tt = peakPos*5; tt > 0. && !foundwbefore; tt-=dx ) {
        if (res->Parameter(0)*TMath::Landau(tt,res->Parameter(1),res->Parameter(2),false) <= pvh) {
        foundwbefore = true;
        wbefore = tt;
        //std::cout << "wbefore: " << wbefore << ", peakPos: " << peakPos << ", phv: " << pvh << ", fit: " << res->Parameter(0)*TMath::Landau(tt,res->Parameter(1),res->Parameter(2),false) << std::endl;
        }
      }

      res = hWidthFit->Fit("landau","QS+","",peakPos*5,(peakPos+60.)*5);
      for(Float_t tt = peakPos*5; tt<numSamples*5 && !foundwafter; tt+=dx ) {
        if (res->Parameter(0)*TMath::Landau(tt,res->Parameter(1),res->Parameter(2),false) <= pvh) {
        foundwafter = true;
        //width = s-peakPos;
        wafter = tt;
        //std::cout << "wafter: " << wafter << ", peakPos: " << peakPos << ", phv: " << pvh << ", fit: " << res->Parameter(0)*TMath::Landau(tt,res->Parameter(1),res->Parameter(2),false) << std::endl;
        }
      }

    }
    /*
    for(int s = peakPos; s < 220 && !foundwafter && snap[s] < kGlobalPed; s++) {
      if ((kGlobalPed-snap[s])<(kGlobalPed-peakVal)*.5) {
        foundwafter = true;
        //width = s-peakPos;
        wafter = s;
        //std::cout << "Width: " << width << ", peakPos: " << peakPos << ", s: " << s << ", peakVal: " << peakVal << ", : " << (kGlobalPed-snap[s]) << ", : " << (kGlobalPed-peakVal)/2. << ", snap[s]: " << snap[s] << std::endl;
      }
    }
    for(int s = peakPos; s>=0 && !foundwbefore && snap[s] < kGlobalPed; s--) {
      if ((kGlobalPed-snap[s])<(kGlobalPed-peakVal)*.5) {
        foundwbefore = true;
        //width = s-peakPos;
        wbefore = s;
        //std::cout << "Width: " << width << ", peakPos: " << peakPos << ", s: " << s << ", peakVal: " << peakVal << ", : " << (kGlobalPed-snap[s]) << ", : " << (kGlobalPed-peakVal)/2. << ", snap[s]: " << snap[s] << std::endl;
      }
    }
    */
    if(foundwafter&&foundwbefore) {
      foundWidth = true;
      width = wafter-wbefore;
    }

    if(foundPeak&&notSat) {
	   hPeakValueVsSum->Fill(sum,peakVal);
     hPeak->Fill(pedsum-peakVal);
        //std::cout << "sum: " << sum << " peak: " << peakVal << std::endl;
    }
    if(foundWidth&&notSat) {
      hPeakValueVsWidth->Fill(width,peakVal);
      hWidthValueVsSum->Fill(sum,width);
      //if(width<75.) {
        //std::cout << "width: " << width  << ", sum: " << sum << " peak: " << peakVal << std::endl;
        //return;
      //}
      hWidth->Fill(width);
    }
  }
  c->cd(0);
  gPad->SetGrid(true,true);
  hPeakValueVsSum->Draw();
  hPeakValueVsSum->Fit("pol1");

  c2->cd(0);
  gPad->SetGrid(true,true);
  hPeakValueVsWidth->Draw();

  c3->cd(0);
  gPad->SetGrid(true,true);
  hWidthValueVsSum->Draw();
  gPad->Update();

  c4->cd(0);
  hWidth->Draw();
  gPad->Update();

  c5->cd(0);
  hPeak->Draw();
  gPad->Update();


  //hPeakPos[0]->Fit("gaus","","",50.,100.);
  //c->SaveAs(TString::Format("results/snapshot_peak_position_vs_time_%d.png",run));

}


