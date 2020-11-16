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

void snapshot_test()
{
}

void config(int run)
{
  if(gLoadedRun != run) {
    if(gChain)
      delete gChain;

    gChain = new TChain("snapshots");
    gChain->Add(TString::Format("/data/cmuwork/rootfiles/prex/compmon_%d.root",run));

    gChain->SetBranchStatus("*",0);
    gChain->SetBranchStatus("numSamples",1);
    gChain->SetBranchStatus("test15",1);
    gChain->SetBranchStatus("randomTime",1);
    gChain->SetBranchStatus("mpsCoda",1);
    gChain->SetBranchStatus("snap",1);
    gChain->SetBranchStatus("snapClock",1);
    //gChain->SetBranchStatus("bcm",1);
    gChain->SetBranchStatus("laserState",1);
    gChain->SetBranchStatus("beamState",1);
    gChain->SetBranchAddress("numSamples",&numSamples);
    gChain->SetBranchAddress("test15",&test15);
    gChain->SetBranchAddress("randomTime",&randomTime);
    gChain->SetBranchAddress("snap",&snap);
    gChain->SetBranchAddress("snapClock",&snapClock);
    gChain->SetBranchAddress("bcm",&bcm);
    gChain->SetBranchAddress("laserState",&laserState);
    //  gChain->SetBranchAddress("beamState",&beamState);
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
    if(randomTime==0 && snapClock < 1666000) {
      std::cerr << "SnapClock = " << snapClock << std::endl;
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

