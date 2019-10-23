const Bool_t kSubHelDiff = true;

TH1F *hON[2] = {0};
TH1F *hOFF[2] = {0};
TCanvas *canvas = 0;
TCanvas *canvasTime = 0;
TCanvas *canvasNAcc = 0;
TCanvas *canvasAsym = 0;
TH1F *hONSum = 0;
TH1F *hONDiff = 0;
TH1F *hONAsym = 0;
TH1F *hOFFSum = 0;
TH1F *hOFFDiff = 0;
TH1F *hOFFAsym = 0;
TLegend *leg[100];
TH2F *hTimeON[5] = {0};
TH2F *hTimeOFF[5] = {0};
TH1F *hNAccON[5] = {0};
TH1F *hNAccOFF[5] = {0};

const char *HelNames[2] = { "POS", "NEG" };
const char *kRootfilePrefix = "/data/cmuwork/rootfiles/prex/compmon_";

TCut kLaserState[2] = { "laserState==1", "laserState==3" };
//TCut kBeamState = "bcm>95&&abs(bpmBy_raw-0.191)<0.002&&dithering==0";
TCut kBeamState = "bcm>95&&dithering==0&&mpsCoda>30e3";
//TCut kBeamState = "bcm>65&&mpsCoda<350e3&&epics_bpmAy<-.400&&epics_bpmAy>-.455";
//TCut kBeamState = "bcm>65&&epics_bpmAy<-.400&&epics_bpmAy>-.455";

void GetAnalyzingPower(int run, double &apower, double &apowerE)
{
  if( run < 4000) {
    apower = 0.1105;
    apowerE = 0.0;
  } else {// if ( run < 4290 ) {
    apower = 0.0166;
    apowerE = 0.0001;
  }
}

void MyRatio(double *v1, double *v2, double *results, double mult = 1.0)
{
  results[0] = mult*v1[0]/v2[0];
  results[1] = results[0]*TMath::Sqrt(
    TMath::Power(v1[1]/v1[0],2)+
    TMath::Power(v2[1]/v2[0],2));
}

void add(double *v1, double *v2, double *results, double mult = 1.0)
{
  results[0] = v1[0] + mult*v2[0];
  results[1] = TMath::Sqrt(TMath::Power(v1[1],2)+TMath::Power(v2[1],2));
}

void printResult(const char *title, double *r)
{
  TString out = TString::Format("%-20s: %15.5f +/- %15.5f", title,r[0],r[1]);
  std::cout << out.Data() << std::endl;
}


void printResult2(const char *title, double val, double err = -1)
{
  TString out = TString::Format("%-40s: %10.5f +/- %10.5f", title,val,err);
  std::cout << out.Data() << std::endl;
}

void MakeHisto(const char *name, const char *title, int nbins, double min, double max, TH1F **h1, TH1F **h2, int l)
{
  TH1F *h = new TH1F(name,title,nbins,min,max);
  if(l==0) {
    *h1 = h;
  } else {
    *h2 = h;
  }
}

void MakeHisto(const char *name, const char *title, int nbinsx, double minx, double maxx,int nbinsy, double miny, double maxy, TH2F **h1, TH2F **h2, int l)
{
  TH2F *h = new TH2F(name,title,nbinsx,minx,maxx,nbinsy,miny,maxy);
  if(l==0) {
    *h1 = h;
  } else {
    *h2 = h;
  }
}


void DrawPair(TH1* h1, TH1* h2, const char *l1, const char *l2, size_t idx)
{
  h1->SetMarkerColor(kGreen+1);
  h2->SetMarkerColor(kRed+1);
  h1->SetLineColor(kGreen+1);
  h2->SetLineColor(kRed+1);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h1->Draw();
  h2->Draw("SAME");
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)gPad->GetPrimitive("stats");
  leg[idx] = new TLegend(ps->GetX1NDC(),ps->GetY1NDC()-0.1,ps->GetX2NDC(),ps->GetY1NDC());
  leg[idx]->AddEntry(h1,l1);
  leg[idx]->AddEntry(h2,l2);
  leg[idx]->Draw();
}

void quickQrtAsym(int run = 2810, int doAccum = 0)
{
  //double  NAcc0 = 6601042;
  TChain *chain = new TChain("quartetwise");
  //chain->Add(TString::Format("/data/cmuwork/rootfiles/compmon_%d.root",run));
  chain->Add(TString::Format("%s%d.root",kRootfilePrefix,run));
  //chain->Add(TString::Format("%s%d.root",kRootfilePrefix,run-1));
  //chain->Add(TString::Format("%s%d.root",kRootfilePrefix,run-2));
  if(!canvas) {
    canvas = new TCanvas("canvas","canvas",1200,1000);
    canvas->Divide(2,2);
  }
  if(!canvasAsym) {
    canvasAsym = new TCanvas("canvasAsym","canvasAsym",600,500*2);
    canvasAsym->Divide(1,2);
  }
  if(!canvasTime) {
    canvasTime = new TCanvas("canvasTime","canvasTime",1200,1000);
    canvasTime->Divide(2,2);
  }
  if(!canvasNAcc && doAccum > 0) {
    canvasNAcc = new TCanvas("canvasNAcc","canvasNAcc",1200,1000);
    canvasNAcc->Divide(2,2);
  }


  TString draw[5] = {"","","","",""};
  TString drawN[5] = {"","","","",""};
  TString units = "";
  int binsTime =  100;
  int minTime = 0;
  int maxTime = chain->GetEntries();
  int binsAcc =  100;
  double minAcc =  -156e6*0 + -2;
  double maxAcc =  -146e6*0 +  2;
  double minSum = 0;
  double maxSum = 0;
  double minDiff = 0;
  double maxDiff = 0;
  double minAsym = 0;
  double maxAsym = 0;
  //chain->GetEntries()/100;
  double minNAcc=0;
  double maxNAcc=0;
  double minNAccSum=0;
  double maxNAccSum=0;
  double minNAccDiff=0;
  double maxNAccDiff=0;
  if(doAccum ==0) {
    minAcc = -0.5*0-1*0+5*0;
    maxAcc =  2.0*0+6*0+20;
    minSum = -.5*0+10*0;
    maxSum = 5.0*0+10*0+35;
    minDiff = -0.9*0-2;
    maxDiff =  0.9*0+2*2;
    units = "rau";
    draw[0] = "PosHelAcc0/PosHelNSamples0";
    draw[1] = "NegHelAcc0/NegHelNSamples0";
    //draw[0] = "100*PosHelAcc0/PosHelNSamples0/bcm";
    //draw[1] = "100*NegHelAcc0/NegHelNSamples0/bcm";
  } else if(doAccum <= 4) {
    minAcc = 15e6*0+100e3*0+1.3e7*0+5e6*0;
    maxAcc = 26e6*0+7000e3*0+2.5e7*2;
    draw[0] = TString::Format("PosHelAcc%d",doAccum);
    draw[1] = TString::Format("NegHelAcc%d",doAccum);
    units = "srau";
    minSum = 100e3*0+2.8e7*0+5e6*0;
    maxSum=  7000e3*4*0+5e7*2;
    minNAcc=0;
    maxNAcc=833e3;
    minNAccSum=900;
    maxNAccSum=833e3*2;
  } else {
    return;
  }

  drawN[0] = Form("PosHelNSamples%d",doAccum);
  drawN[1] = Form("NegHelNSamples%d",doAccum);
  drawN[2] = TString::Format("(%s+%s)",drawN[0].Data(),drawN[1].Data());
  drawN[3] = TString::Format("(%s-%s)",drawN[0].Data(),drawN[1].Data());
  drawN[4] = TString::Format("%s/%s",drawN[3].Data(),drawN[2].Data());

  minAsym =  -0.13*0-0.4;
  maxAsym =   0.13*0-0.4;
  draw[2] = TString::Format("(%s+%s)",draw[0].Data(),draw[1].Data());
  draw[3] = TString::Format("(%s-%s)",draw[0].Data(),draw[1].Data());
  draw[4] = TString::Format("%s/%s",draw[3].Data(),draw[2].Data());
  canvas->cd(1);
  // Check if these histograms are already created
  TString laserNames[2] = { "ON", "OFF" };

  TH1F *tmpH1;
  TH2F *tmpH2;
  for(Int_t l = 0; l < 2; l++) {
    for(Int_t h = 0; h < 2; h++) {
      MakeHisto(TString::Format("h%s%d",laserNames[l].Data(),h),
          //TString::Format("Run %d Acc%d Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d %s-Hel; Acc%d [%s]",run,doAccum,HelNames[h],doAccum,units.Data()),
          binsAcc,minAcc,maxAcc,&hON[h],&hOFF[h],l);
      MakeHisto(TString::Format("hTime%s%d",laserNames[l].Data(),h),
          //TString::Format("Run %d Acc%d Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d %s-Hel; Quartet Number; Acc%d [%s]",run,doAccum,HelNames[h],doAccum,units.Data()),
          binsTime,minTime,maxTime,binsAcc,minAcc,maxAcc,&hTimeON[h],&hTimeOFF[h],l);

      MakeHisto(TString::Format("hNAcc%s%d",laserNames[l].Data(),h),
          TString::Format("Run %d NAcc%d %s-Hel; NAcc%d",run,doAccum,HelNames[h],doAccum),
          binsAcc,minNAcc,maxNAcc,&hNAccON[h],&hNAccOFF[h],l);

      // Draw yields
      chain->Draw(TString::Format("%s>>h%s%d",draw[h].Data(),laserNames[l].Data(),h),
          kLaserState[l]&&kBeamState);
      chain->Draw(TString::Format("%s:Entry$>>hTime%s%d",draw[h].Data(),laserNames[l].Data(),h),
          kLaserState[l]&&kBeamState);
      chain->Draw(TString::Format("%s>>hNAcc%s%d",drawN[h].Data(),laserNames[l].Data(),h),
          kLaserState[l]&&kBeamState);

    }
    MakeHisto(TString::Format("h%sSum",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Sum Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Sum; Sum [%s]",run,doAccum,units.Data()),
          binsAcc,minSum,maxSum, &hONSum, &hOFFSum, l);
    MakeHisto(TString::Format("h%sDiff",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Diff Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Diff;Diff [%s]",run,doAccum,units.Data()),
          binsAcc,minDiff,maxDiff, &hONDiff, &hOFFDiff, l);
    MakeHisto(TString::Format("h%sAsym",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Asym Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Asym; Asymmetry",run,doAccum),
          binsAcc,minAsym,maxAsym, &hONAsym, &hOFFAsym, l);
    MakeHisto(TString::Format("hTime%sSum",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Sum Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Sum;Quartet Number;Sum [%s]",run,doAccum,units.Data()),
          binsTime,minTime,maxTime,binsAcc,minSum,maxSum, &hTimeON[2], &hTimeOFF[2], l);
    MakeHisto(TString::Format("hTime%sDiff",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Diff Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Diff;Quartet Number;Diff [%s]",run,doAccum,units.Data()),
          binsTime,minTime,maxTime,binsAcc,minDiff,maxDiff, &hTimeON[3], &hTimeOFF[3], l);
    MakeHisto(TString::Format("hTime%sAsym",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Asym Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Asym;Quartet Number;Asym",run,doAccum),
          binsTime,minTime,maxTime,binsAcc,minAsym,maxAsym, &hTimeON[4], &hTimeOFF[4], l);

    MakeHisto(TString::Format("hNAcc%sSum",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Sum Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d NAcc%d Sum; Sum",run,doAccum),
          binsAcc,minNAccSum,maxNAccSum, &hNAccON[2], &hNAccOFF[2], l);
    MakeHisto(TString::Format("hNAcc%sDiff",laserNames[l].Data()),
          //TString::Format("Run %d Acc%d Diff Laser %s",run,doAccum,laserNames[l].Data()),
          TString::Format("Run %d Acc%d Diff;Diff",run,doAccum),
          binsAcc,minNAccDiff,maxNAccDiff, &hNAccON[3], &hNAccOFF[3], l);


    //chain->Draw(TString::Format("%s>>hOFF%d",draw[h].Data(),h),kLaserState[1]&&kBeamState);
    chain->Draw(TString::Format("%s>>h%sSum",draw[2].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);
    chain->Draw(TString::Format("%s:Entry$>>hTime%sSum",draw[2].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);
    chain->Draw(TString::Format("%s>>hNAcc%sSum",drawN[2].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);

    chain->Draw(TString::Format("%s>>h%sDiff",draw[3].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);
    chain->Draw(TString::Format("%s:Entry$>>hTime%sDiff",draw[3].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);
    chain->Draw(TString::Format("%s>>hNAcc%sDiff",drawN[3].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);

    chain->Draw(TString::Format("%s>>h%sAsym",draw[4].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);
    chain->Draw(TString::Format("%s:Entry$>>hTime%sAsym",draw[4].Data(),laserNames[l].Data()),
        kLaserState[l]&&kBeamState);

  }

  hON[0]->SetLineColor(kRed+1);
  hON[1]->SetLineColor(kGreen+1);
  hOFF[0]->SetLineColor(kRed+1);
  hOFF[1]->SetLineColor(kGreen+1);
  hONSum->SetLineColor(kGreen+1);
  hONDiff->SetLineColor(kGreen+1);
  hONAsym->SetLineColor(kGreen+1);
  hOFFSum->SetLineColor(kRed+1);
  hOFFDiff->SetLineColor(kRed+1);
  hOFFAsym->SetLineColor(kRed+1);

  canvas->cd(1);
  int idx = 0;
  DrawPair(hON[0],hOFF[0], "Laser ON","Laser OFF",idx++);
  canvas->cd(2);
  DrawPair(hON[1],hOFF[1], "Laser ON","Laser OFF",idx++);
  canvas->cd(3);
  DrawPair(hONSum,hOFFSum,"Laser ON", "Laser OFF",idx++);
  canvas->cd(4);
  DrawPair(hONDiff,hOFFDiff,"Laser ON", "Laser OFF",idx++);
  //
  canvasAsym->cd(1);
  DrawPair(hONAsym,hOFFAsym,"Laser ON", "Laser OFF",idx++);
  canvasAsym->cd(2);
  DrawPair(hTimeON[4],hTimeOFF[4],"Laser ON", "Laser OFF",idx++);

  canvasTime->cd(1);
  DrawPair(hTimeON[0],hTimeOFF[0],"Laser ON","Laser OFF",idx++);
  canvasTime->cd(2);
  DrawPair(hTimeON[1],hTimeOFF[1],"Laser ON", "Laser OFF",idx++);
  canvasTime->cd(3);
  DrawPair(hTimeON[2],hTimeOFF[2],"Laser ON", "Laser OFF",idx++);
  canvasTime->cd(4);
  DrawPair(hTimeON[3],hTimeOFF[3],"Laser ON", "Laser OFF",idx++);

  if( doAccum > 0 ) {
    canvasNAcc->cd(1);
    DrawPair(hNAccON[0],hNAccOFF[0],"Laser ON","Laser OFF",idx++);
    canvasNAcc->cd(2);
    DrawPair(hNAccON[1],hNAccOFF[1],"Laser ON", "Laser OFF",idx++);
    canvasNAcc->cd(3);
    DrawPair(hNAccON[2],hNAccOFF[2],"Laser ON", "Laser OFF",idx++);
    canvasNAcc->cd(4);
    DrawPair(hNAccON[3],hNAccOFF[3],"Laser ON", "Laser OFF",idx++);
  }



  double vOn[2][2];
  double vOnSum[2];
  double vOnDiff[2];
  double vOnAsym[2];
  
  double vOff[2][2];
  double vOffSum[2];
  double vOffDiff[2];
  double vOffAsym[2];
  double vOnNAccSum[2];
  double vOnNAccDiff[2];
  double vOffNAccSum[2];
  double vOffNAccDiff[2];
  

  for(int i = 0; i < 2; i++) {
    vOn[i][0] = hON[i]->GetMean();
    vOn[i][1] = hON[i]->GetMeanError();
    vOff[i][0] = hOFF[i]->GetMean();
    vOff[i][1] = hOFF[i]->GetMeanError();
  }

  //add(vOn[1],vOn[0],vOnSum,1.0);
  //add(vOn[1],vOn[0],vOnDiff,-1.0);
  vOnSum[0]=hONSum->GetMean();
  vOnSum[1]=hONSum->GetMeanError();
  vOnDiff[0]=hONDiff->GetMean();
  vOnDiff[1]=hONDiff->GetMeanError();
  vOnNAccSum[0]=hNAccON[2]->GetMean();
  vOnNAccDiff[0]=hNAccON[3]->GetMean();
  TSpectrum *sOnSum = new TSpectrum(2);
  //sOnSum->Search(hONSum);
  //TFitResultPTr *rOnSum=hONSum->Fit("gaus","S","",sOnSum->GetPeaks()[0]<
  MyRatio(vOnDiff,vOnSum,vOnAsym);

  //add(vOff[1],vOff[0],vOffSum,1.0);
  //add(vOff[1],vOff[0],vOffDiff,-1.0);
  vOffSum[0]=hOFFSum->GetMean();
  vOffSum[1]=hOFFSum->GetMeanError();
  vOffDiff[0]=hOFFDiff->GetMean();
  vOffDiff[1]=hOFFDiff->GetMeanError();
  vOffNAccSum[0]=hNAccOFF[2]->GetMean();
  vOffNAccDiff[0]=hNAccOFF[3]->GetMean();
  MyRatio(vOffDiff,vOffSum,vOffAsym);

  double vBkSubSum[2];
  double vBkSubDiff[2];
  double vBkSubAsym[2];
  add(vOnSum,vOffSum,vBkSubSum,-1.0);
  if(kSubHelDiff) {
    add(vOnDiff,vOffDiff,vBkSubDiff,-1.0);
  } else {
    vBkSubDiff[0] = vOnDiff[0];
    vBkSubDiff[1] = vOnDiff[1];
  }
  MyRatio(vBkSubDiff,vBkSubSum,vBkSubAsym);

  std::cout << "----------------------" << std::endl;
  printResult("Laser-ON  H0",vOn[0]);
  printResult("Laser-ON  H1",vOn[1]);
  printResult("Laser-ON  H1+H0",vOnSum);
  printResult("Laser-ON  H1-H0",vOnDiff);
  printResult("Laser-ON  Asym",vOnAsym);
  std::cout << "----------------------" << std::endl;
  printResult("Laser-OFF H0",vOff[0]);
  printResult("Laser-OFF H1",vOff[1]);
  printResult("Laser-OFF H1+H0",vOffSum);
  printResult("Laser-OFF H1-H0",vOffDiff);
  printResult("Laser-OFF Asym",vOffAsym);
  std::cout << "----------------------" << std::endl;
  printResult("BkSub H1+H0",vBkSubSum); 
  printResult("BkSub H1-H0",vBkSubDiff); 
  printResult("BkSub Asym",vBkSubAsym);

  //double apower = .1105;
  double apower[2];
  GetAnalyzingPower(run,apower[0],apower[1]);
  double pol[2];
  MyRatio(vBkSubAsym,apower,pol,100.);
  printResult("Polarization",pol);

  if( doAccum > 0 ) {
    // Now see what happens if we get the pedestal off by:
    std::vector<double> pedOff;
    pedOff.push_back(-1.0);
    pedOff.push_back(+1.0);
    pedOff.push_back(-0.5);
    pedOff.push_back(+0.5);
    pedOff.push_back(+0.1);
    pedOff.push_back(-0.1);
    pedOff.push_back(+0.05);
    pedOff.push_back(-0.05);
    double ped = 0;
    double corr = 0;
    double vBkSubSumPed[2];
    double vBkSubDiffPed[2];
    double vBkSubAsymPed[2];
    double vBkNAccSubSum[2];
    double vBkNAccSubDiff[2];
    double polPed[2];
    vBkNAccSubSum[0]=vOnNAccSum[0]-vOffNAccSum[0];
    vBkNAccSubSum[1] =0;
    vBkSubSumPed[1] = vBkSubSum[1];
    vBkSubDiffPed[1] = vBkSubDiff[1];
    vBkNAccSubDiff[0]=vOnNAccDiff[0]-vOffNAccDiff[0];
    vBkNAccSubDiff[1] =0;
    std::cout << "NAcc ON:     " << int(vOnNAccSum[0]) << std::endl;
    std::cout << "NAcc OFF:    " << int(vOffNAccSum[0]) << std::endl;
    std::cout << "Nacc ON-OFF: " << int(vBkNAccSubSum[0]) << std::endl;
    for(size_t k = 0 ; k < pedOff.size(); k++) {
      ped=pedOff[k];
      corr = vBkNAccSubSum[0]*ped;
      vBkSubDiffPed[0] = vBkSubDiff[0]-vBkNAccSubDiff[0]*ped;
      vBkSubSumPed[0] = vBkSubSum[0]-vBkNAccSubSum[0]*ped;
      MyRatio(vBkSubDiffPed,vBkSubSumPed,vBkSubAsymPed);
      MyRatio(vBkSubAsymPed,apower,polPed,100.);
      polPed[1] = 0;
      std::cout << "If pedestal is off by: " << ped << std::endl;
      //printResult("Polarization",polPed);
      printf("Polarization %15.5f \t--> %15.5f%% (absolute pol change)\n",polPed[0],(polPed[0]-pol[0]));
      //vBkSubSumPed[0] = vBkSubSum[
    }
  }
  
}
