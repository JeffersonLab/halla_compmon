
#include "TF1.h"


//---- a class to parameterize a histogram -------------

class HistoFunctionParam{
public:
  void HistoFunctionParam(TH1F* datahisto){buffhisto=datahisto;};
  double Evaluate(double *x, double* p){
    double xvalue = p[1]*(x[0]+p[3]);
    double yvalue = p[2]+p[0]*yInterpolate(buffhisto,xvalue);
    return yvalue;
  };
  void SetHisto(TH1F* datahisto){
    buffhisto = datahisto;
  };
private:
  TH1F* buffhisto;
};
//----- a function to turn a histogram's data into a function y(x)-----

Double_t yInterpolate(TH1F* datahisto, Double_t x){
  Double_t delx = x - datahisto->GetBinLowEdge(1);  //finds the difference between x value and first bin x value
  Double_t width = datahisto->GetBinWidth(1); //finds the bin width
  Int_t nbins = datahisto->GetNbinsX(); //finds total number of bins
  Int_t xbin = delx/width+1;  //finds the bin that x values is located in
  Double_t xbincenter = datahisto->GetBinCenter(xbin); //finds the center of bin x is located in
  Double_t xinterp = x - xbincenter;  //finds difference between x and xbincenter
  Double_t y;

  if(xinterp<0 && xbin>1 && xbin<nbins){
    y = ((width + xinterp)*datahisto->GetBinContent(xbin) + (-xinterp)*datahisto->GetBinContent(xbin-1))/width;
  }elseif(xinterp>0.0 && xbin>0 && xbin<nbins-1){
    y = ((width - xinterp)*datahisto->GetBinContent(xbin) + (xinterp)*datahisto->GetBinContent(xbin+1))/width;
  }elseif(xbin>0 && xbin<nbins){
    y = datahisto->GetBinContent(xbin);
  }else{
    y = 0.;
  }
  return y;
}

//--------- a function to establish errors for the histogram hasymdata-------------

void SetAsymErrors(TH1F* htotal, TH1F* hasymdata, Double_t factor){
  Int_t nbins = htotal->GetNbinsX();
  Double_t ncounts;
  for(Int_t i=0; i<nbins; i++){
    ncounts = factor*htotal->GetBinContent(i);
    if(ncounts>0){
      hasymdata->SetBinError(i,1.0/sqrt(ncounts));
    }else{
      hasymdata->SetBinError(i,0.0);
      hasymdata->SetBinContent(i,0.0);
    }
  }
  return;
}

void FitAsymmetry(TF1* Fit, TH1F* hData, Double_t xscale, Double_t xoffset, Double_t CavityPolarization, Double_t fadcHi, Double_t AsymFitRange_lo, Double_t AsymFitRange_hi){
  //fits the function Fit, to the data in the histogram hData
  //begin by setting up the hData histogram

  hData->Rebin(10);
  hData->Scale(1./10.);
  hData->GetXaxis()->SetRangeUser(1000.,fadcHi);
  hData->SetMaximum(0.1);
  hData->SetMinimum(-0.05);
  hData->SetLineColor(kGreen);

  //set up the fit function

  Fit->SetParameters(1.0,1.0,0.0,0.0);
  Fit->SetParNames("vScale","xScale","yOffset","xOffset");
  Fit->SetRange(AsymFitRange_lo,AsymFitRange_hi);
  Fit->FixParameter(1,xscale);
  Fit->FixParameter(3,xoffset);
  Fit->SetParLimits(0,0.1,2.);
  Fit->SetLineWidth(1);

  cout<<"************************************"<<endl;

  hData->Fit(Fit,"rn");
  Double_t yscale = Fit->GetParameter(0);
  Double_t yerror = Fit->GetParError(0);
  Double_t xscale = Fit->GetParameter(1);
  printf("\n");
  printf("     Fit to Combined Polarization %10.4f +- %10.4f\n",yscale,yerror);
  printf("     Cavity Polarization (input) %10.4f \n",CavityPolarization);
  printf("     Beam Polarization (extracted)  %10.4f +- %10.4f\n",yscale/CavityPolarization,yerror/CavityPolarization);
  printf("     xscale =                  %10.4f\n",Fit->xscale);
  printf("     Chisqr =                  %10.4f\n",Fit->GetChisquare());
  printf("     Number of data points in fit is  %10.4f\n",Fit->GetNumberFitPoints());
  hdata->Draw();
  Fit->DrawClone("same");
  return;
}

void ClipZeros(char* vartoclip, double numtoclip, int bufferlength){
  char*pos;
  int length;

  snprintf(vartoclip,bufferlength,"%lf\n",numtoclip);
  length=strlen(vartoclip);
  pos = vartoclip+length-1;

  while(*pos == '0')
    *pos-- = '\0';
  if(*pos == '.')
    *pos = '\0';
}

//----------------- now begins the main function, ComptonFit_Alexa, that calls the above functions.------
//---2509 is the run number from Compton commissioning
 void  ComptonFit_Asymmetry(Int_t RunNumber){
  
  void ClipZeros(char*vartoclip, double numtoclip, int bufferlength);
  Double_t xhi = 250;
  Double_t fadcHi = 25000;
  Double_t xRescale = 4000.0; 
  Int_t nbins = 500;
  Double_t asym_min = -0.3;
  Double_t asym_max = 0.5;
  //----------- pick different fitting parameters for different kinematics -----------
  if (RunNumber <= 2570){
  Double_t Fitrange_lo = 1800;
  Double_t Fitrange_hi = 4200;
  Double_t RightCavityPolarization;  //  laser only left polarized in spring 2016
  Double_t LeftCavityPolarization = 1.0;
  Double_t shieldthickness = 2; // mm 
  Double_t beamenergy = 4.48; //GeV 
  Double_t tail_lo = 6000;
  Double_t tail_hi = 10000;
  Double_t half_edge = 3500;
  Double_t mcxscale = 1/4000.;
  Double_t mcxoffset = -300.;
  Double_t ymax = 8;
  Double_t ymin = -0.2;
  Double_t legxmin = 0.60;
  Double_t legxmax = 0.90;
  Double_t mc_halfedge = 0.95;
  Double_t res = 0.11;


  }
  if(RunNumber>=2577 && RunNumber<2614){
  Double_t Fitrange_lo = 8300;
  Double_t Fitrange_hi = 20000;
  Double_t RightCavityPolarization;  //  laser only left polarized in spring 2016
  Double_t LeftCavityPolarization = 1.0;
  Double_t shieldthickness = 2; // mm 
  Double_t beamenergy = 8.8; //GeV 
  Double_t tail_lo = 23000;
  Double_t tail_hi = 25000;
  Double_t half_edge = 18750;
  Double_t mcxscale = 1/13000.;
  Double_t mcxoffset = -7000.;
  Double_t ymax = 6;
  Double_t ymin = -0.5;
  Double_t legxmin = 0.15;
  Double_t legxmax = 0.45;
  Double_t mc_halfedge = 0.95;
  Double_t res = 0.06;

  }
  if(RunNumber>=2614 && RunNumber < 2647){
  Double_t Fitrange_lo = 2000;
  Double_t Fitrange_hi = 13000;
  Double_t RightCavityPolarization;  //  laser only left polarized in spring 2016
  Double_t LeftCavityPolarization = 1.0;
  Double_t shieldthickness = 6; // mm 
  //Double_t shieldthickness = 2; // see if this fits data better
  Double_t beamenergy = 8.8; //GeV
  Double_t tail_lo = 15000;
  Double_t tail_hi = 20000;
  Double_t half_edge = 13000;
  Double_t mcxscale = 1/15000.;
  Double_t mcxoffset = 0.;
  Double_t ymax = 6;
  Double_t ymin = -0.5;
  Double_t legxmin = 0.15;
  Double_t legxmax = 0.45;
  Double_t mc_halfedge = 0.925;
  Double_t res = 0.055;

	 if(RunNumber == 2640)
	   res = 0.065;
	 if(RunNumber == 2641)
	   res = 0.065;
	 if(RunNumber == 2642)
	   res = 0.08;
	 if(RunNumber == 2643)
	   res = 0.08;
	 if(RunNumber == 2644)
	   res = 0.08;
	 if(RunNumber == 2645)
	   res = 0.105;
	 if(RunNumber == 2646)
	   res = 0.103;

  }
 if(RunNumber>=2647){
  Double_t Fitrange_lo = 2000;
  Double_t Fitrange_hi = 13000;
  Double_t RightCavityPolarization;  //  laser only left polarized in spring 2016
  Double_t LeftCavityPolarization = 1.0;
  Double_t shieldthickness = 6; // mm 
  //Double_t shieldthickness = 2; // see if this fits data better
  Double_t beamenergy = 8.8; //GeV
  Double_t tail_lo = 15000;
  Double_t tail_hi = 20000;
  Double_t half_edge = 12000;
  Double_t mcxscale = 1/15000.;
  Double_t mcxoffset = 0.;
  Double_t ymax = 6;
  Double_t ymin = -0.5;
  Double_t legxmin = 0.15;
  Double_t legxmax = 0.45;
  Double_t mc_halfedge = 0.925;
  Double_t res = 0.055;
  }

if(RunNumber>=2751 && RunNumber <= 2753){
  Double_t Fitrange_lo = 2000;
  Double_t Fitrange_hi = 18000;
  Double_t RightCavityPolarization;  //  laser only left polarized in spring 2016
  Double_t LeftCavityPolarization = 1.0;
  Double_t shieldthickness = 6; //mm
  //Double_t shieldthickness = 2; // see if this fits data better
  Double_t beamenergy = 11.0; //GeV
  Double_t tail_lo = 15000;
  Double_t tail_hi = 30000;
  Double_t half_edge = 17500;
  Double_t mcxscale = 1/15000.;
  Double_t mcxoffset = 0.;
  Double_t ymax = 6;
  Double_t ymin = -0.5;
  Double_t legxmin = 0.15;
  Double_t legxmax = 0.45;
  Double_t mc_halfedge = 0.925;
  Double_t res = 0.08;
  }



//char beam[8];
  char shield[8];
  char mc_res[8];
  ClipZeros(shield,shieldthickness,8);
  //ClipZeros(beam,beamenergy,8);
  ClipZeros(mc_res,res,8);


  //-------------------------------------------------------------------
  //-------- get MC data from the output of MC_hist_generator.C ------------
  //-------------------------------------------------------------------

  char mcfilename[60];
  sprintf(mcfilename,"/home/compton/lthorne/CompMon/AnalyzingPower/root/MC_%s_11_%s.root", shield, mc_res);
  TFile* MCFile = new TFile(mcfilename);
  if(! MCFile->IsOpen()){
    printf("MC file is already open somewhere else, unable to open MC file.\n");
    return;
  }
  TH1F* MCSum = (TH1F*) MCFile->Get("sum");
  TH1F* MCDiff = (TH1F*) MCFile->Get("diff");
  TH1F* MCAsym = (TH1F*) MCFile->Get("div");
  TH1F* MCAlgnd = (TH1F*) MCFile->Get("algnd");
  TH1F* MCAnti = (TH1F*) MCFile->Get("anti");

  Double_t y = MCSum->GetBinContent(5); //cut out spike near zero
  MCSum->SetMaximum(1.5*y);
  printf("Successfully opened MC file and retrieved histograms beam 11GeV, shield %smm \n", shield);

  /* TCanvas*c_mc =new TCanvas("c_mc");
  c_mc->Divide(1,3);
  c_mc->cd(1);
  MCSum->Draw();
  MCSum->GetXaxis()->SetRangeUser(-0.1,1.5);
  c_mc->cd(2);
  MCDiff->Draw();
  c_mc->cd(3);
  MCAsym->Draw(); */

  /* TCanvas*mc = new TCanvas("mc");
  mc->Divide(1,2);
  mc->cd(1);
  MCAlgnd->Draw();
  mc->cd(2);
  MCAnti->Draw(); */


  
  //-----------------------------------------------------------------------
  //---- get the real data from root file created after using compmon-----------------------------------------
  //--- looks for the rootfile that compmon in /home/compton/franklin/CompMon/ creates in /rootfiles dirctory
  //----------------------------------------------------------------------
  char datafilename[50];
  sprintf(datafilename,"/home/compton/franklin/CompMon/rootfiles/compmon_%d.root",RunNumber);
  TFile* DataFile = new TFile(datafilename);
  if(!DataFile->IsOpen()){
    printf("Unable to open data file, data file is open somewhere else.\n");
    return;
  }
  printf("Successfully opened the data file for run number %d \n", RunNumber);

  DataFile->cd("triggeredHistos");
  DataSumOn_L_P = (TH1F*) hTrig_sums_L_P->Clone("");
  DataSumOn_L_N = (TH1F*) hTrig_sums_L_N->Clone("");
  DataSumOff = (TH1F*) hTrig_sums_laserOff->Clone("");
  DataSumOn_R_P = (TH1F*) hTrig_sums_R_P->Clone("");
  DataSumOn_R_N = (TH1F*) hTrig_sums_R_N->Clone("");
  DataSumOn_L_P->Sumw2();
  DataSumOn_L_N->Sumw2();
  DataSumOff->Sumw2();

  DataFile->cd("quartetHistos");
  //BeamAsymHisto = (TH1F*) hQ_BCM_Asym->Clone("");
  //Double_t beam_asym = BeamAsymHisto->GetMean();
  Double_t beam_asym = 0;
  //cout<< "beam charge asymmetry is = "<<beam_asym<<endl;
  
  Int_t bin_lo = DataSumOn_L_P->FindBin(tail_lo);
  Int_t bin_hi = DataSumOn_L_P->FindBin(tail_hi);
  Double_t Tail_Laser_On_L_P = DataSumOn_L_P->Integral(bin_lo,bin_hi);
  Double_t Tail_Laser_On_L_N = DataSumOn_L_N->Integral(bin_lo,bin_hi);
  Double_t Tail_Laser_Off = DataSumOff->Integral(bin_lo,bin_hi);
  //Double_t bg_scale_L_P = Tail_Laser_On_L_P/Tail_Laser_Off;
  //Double_t bg_scale_L_N = Tail_Laser_On_L_N/Tail_Laser_Off;

  //-------- Finding the correct background scaling  (scaling by accepted triggers)--------

  DataFile->cd("mpsHistos");
  int OnState[] = {0,1,1,0,2,2,2,2};
  int MPS[3];
  int Triggers[3];
  int AcceptedTriggers[3];

  for(int i=0; i<3; i++){
    MPS[i]=0;
    Triggers[i]=0;
    AcceptedTriggers[i]=0;
  }

  TH1F*h_Trigger = hM_Trigs_Prescaled;
  int state;
  float scale;
  for(int i=0;i<=7;i++){
    state = OnState[i];
    if(state>=0){
      MPS[state]+=hM_SpinState_BeamOn->GetBinContent(i+1);
      AcceptedTriggers[state]+=hM_Trig_Accepted->GetBinContent(i+1);
      Triggers[state]+=h_Trigger->GetBinContent(i+1);
    }
  }
  DataSumOn_L_P->Sumw2();
  DataSumOn_L_N->Sumw2();
  for(int i=0; i<3; i++){
    if(MPS[i]>0 && AcceptedTriggers[i]>0){
      float scale = Triggers[i];
      scale = float(Triggers[i])/(float(MPS[i])*float(AcceptedTriggers[i]));
      printf(" triggers %d  mps %d   acceptedtriggers  %d and i is %d\n", Triggers[i], MPS[i], AcceptedTriggers[i],i);
      if(i==0){
	scale = (1-beam_asym)*scale;
	DataSumOn_L_N->Scale(scale);
	DataSumOn_L_N->GetXaxis()->SetRangeUser(-1000.,25000);
	Double_t bg_scale_L_N_error = scale*sqrt(1./(Tail_Laser_On_L_N + Tail_Laser_Off)); 
      }     
      if(i==1){
	scale = (1+beam_asym)*scale;
	DataSumOn_L_P->Scale(scale);
	DataSumOn_L_P->GetXaxis()->SetRangeUser(-1000.,25000);
	Double_t bg_scale_L_P_error = scale*sqrt(1./(Tail_Laser_On_L_P + Tail_Laser_Off));  
      }    
      if(i==2){
	DataSumOff->Scale(scale);
      }
    }       
  }
  
  /* TCanvas*c_data = new TCanvas("c_data");
  c_data->Divide(1,2);
  c_data->cd(1);
  DataSumOn_L_P->Draw();
  c_data->cd(2);
  DataSumOn_L_N->Draw(); */
 
  if(RunNumber>2567){
  TH1F* DataSum = (TH1F*) DataSumOn_L_N->Clone("");
  TH1F* DataDiff = (TH1F*) DataSumOn_L_N->Clone("");
  DataDiff->Add(DataSumOn_L_P,-1.0);
  DataSum->Add(DataSumOn_L_P,+1.0);
  }
  if(RunNumber<=2567){
  TH1F* DataSum = (TH1F*) DataSumOn_L_P->Clone("");
  TH1F* DataDiff = (TH1F*) DataSumOn_L_P->Clone("");
  DataDiff->Add(DataSumOn_L_N,-1.0);
  DataSum->Add(DataSumOn_L_N,+1.0);
  }
  if(RunNumber>2614){
  TH1F* DataSum = (TH1F*) DataSumOn_L_N->Clone("");
  TH1F* DataDiff = (TH1F*) DataSumOn_L_N->Clone("");
  DataDiff->Add(DataSumOn_L_P,-1.0);
  DataSum->Add(DataSumOn_L_P,+1.0);
  }



  DataSum->Add(DataSumOff,-2.0);
  
  TH1F* DataSum1 = (TH1F*) DataSum->Clone("");
  //DataDiff->Add(DataSumOff,-2.0);
  TH1F* DataAsym = (TH1F*) DataDiff->Clone("");
  DataAsym->Divide(DataSum1);


  TCanvas*c_test = new TCanvas("c_test");
  c_test->Divide(1,3);
  c_test->cd(1);
  DataSumOn_L_P->Draw();

  c_test->cd(2);
  DataSumOn_L_N->Draw();

  c_test->cd(3);
  DataSumOff->Draw();



  /*TCanvas*c_asym = new TCanvas("c_asym");
  c_asym->Divide(1,3);
  c_asym->cd(1);
  DataSum1->Draw();
  DataSum1->SetMinimum(0.);
  DataSum1->SetMaximum(2.);
  c_asym->cd(2);
  DataDiff->Draw();
  DataDiff->SetMinimum(-0.1);
  DataDiff->SetMaximum(1.5);
  c_asym->cd(3);
  DataAsym->Draw();
  DataAsym->SetMinimum(-0.2);
  DataAsym->SetMaximum(0.2); */

  //------------------first fit unpolrized data (N+P)-----------------------------------------
  TF1* UnpolarizedFit->Clear();
  HistoFunctionParam* UnpolarizedFunction = new HistoFunctionParam(MCSum);
  TF1* UnpolarizedFit = new TF1("UnpolarizedFit",UnpolarizedFunction, &HistoFunctionParam::Evaluate,Fitrange_lo,Fitrange_hi,4,"HistoFunctionParam","Evaluate");

 
  Double_t yDataUnpolarized = yInterpolate(DataSum,half_edge);
  Double_t yMCUnpolarized = yInterpolate(MCSum,mc_halfedge);
  Double_t UnpolarizedScale = yDataUnpolarized/yMCUnpolarized;
  printf(" The predicted y scale is = %f\n",UnpolarizedScale);


  //Double_t UnpolarizedScale = 2/3000;
  UnpolarizedFit->SetParNames("yScale","xScale","yOffset","xOffset");
  UnpolarizedFit->SetParameter(0,UnpolarizedScale);
  UnpolarizedFit->SetParLimits(0,0.5*UnpolarizedScale,1.5*UnpolarizedScale);
  UnpolarizedFit->FixParameter(2,0.0);
  UnpolarizedFit->SetParameter(3,mcxoffset);
  UnpolarizedFit->SetParameter(1,mcxscale);
  UnpolarizedFit->SetParLimits(1,0.5*mcxscale,1.5*mcxscale);
  UnpolarizedFit->SetParLimits(3,0.5*mcxoffset,1.5*mcxoffset);
  UnpolarizedFit->SetRange(Fitrange_lo,Fitrange_hi);
  UnpolarizedFit->SetLineColor(kMagenta);
  UnpolarizedFit->SetLineWidth(1.0);

 TCanvas* c_fit = new TCanvas("c_fit");
  c_fit->Divide(1,3);
  c_fit->cd(1);


  DataSum->GetXaxis()->SetRangeUser(0.,tail_hi);
  DataSum->SetMinimum(ymin);
  DataSum->SetMaximum(ymax);
  cout<<"************************************"<<endl;
  cout<<"Fit to Summed Spectrum"<<endl;
  DataSum->Fit("UnpolarizedFit","rn");

  Double_t yscale = UnpolarizedFit->GetParameter(0);
  Double_t xscale = UnpolarizedFit->GetParameter(1);
  Double_t xoffset = UnpolarizedFit->GetParameter(3);

  printf(" Unpolarized fit y scaling = %f\n",yscale);
  printf(" Unpolarized fit x scaling = %f\n",xscale);
  printf(" Unpolarized fit x offset = %f\n",xoffset);
  printf(" Unpolarized fit Chi sqr = %f\n",UnpolarizedFit->GetChisquare());
  printf(" Unpolarized fit number of fit points = %f\n",UnpolarizedFit->GetNumberFitPoints());

  DataSum->Draw();
  UnpolarizedFit->DrawCopy("same");
  // UnpolarizedFit->SetRange(Fitrange_hi,25000);
  //UnpolarizedFit->SetLineColor(kRed);
  //UnpolarizedFit->DrawCopy("same");
  //UnpolarizedFit->SetRange(500,Fitrange_lo);
  //UnpolarizedFit->DrawCopy("same"); 




  /* c_unpolarized -> cd(2);
  Double_t x,y,yerr,yfit,p[4];
  UnpolarizedFit->GetParameters(p);
  TH1F*UnpolarizedResiduals = (TH1F*)DataSum->Clone("UnpolarizedResiduals");
  Int_t nbins = UnpolarizedResiduals->GetNbinsX();
  Double_t yfitMidValue = UnpolarizedFit->Eval(UnpolarizedFitMiddle);
  printf(" Unpolarized fit middle %f, y fit middle value %f \n",yfitMidValue,UnpolarizedFitMiddle);
  for(Int_t i=1;i<nbins+1;i++){
    x = UnpolarizedResiduals->GetBinCenter(i);
    yfit = UnpolarizedFit->Eval(x);
    y = UnpolarizedResiduals->GetBinContent(i) - yfit;
    yerr = UnpolarizedResiduals->GetBinError(i);
    if(yfit>0){
      UnpolarizedResiduals->SetBinContent(i,y/yfitMidValue);
      UnpolarizedResiduals->SetBinError(i,yerr/yfitMidValue);
    }else{
      UnpolarizedResiduals->SetBinContent(i,0.);
      UnpolarizedResiduals->SetBinError(i,0.);
    }
  }
  UnpolarizedResiduals->SetMaximum(0.1);
  UnpolarizedResiduals->SetMinimum(-0.1);
  UnpolarizedResiduals->Draw(); */
  

  // c_unpolarized->cd(1);
  //UnpolarizedFit->FixParameter(3,0.0);
  //cout<<"Fit to Unpolarized Spectrum with x0 = 0"<<endl;
  //UnpolarizedFit->SetRange(0.,fadcHi);
  //UnpolarizedFit->SetLineColor(kGreen);

  //------- Fit for hist differences  (N - P) ------------

  TF1* DiffFit->Clear();
  HistoFunctionParam* DiffFunction = new HistoFunctionParam(MCDiff);
  TF1* DiffFit = new TF1("DiffFit",DiffFunction, &HistoFunctionParam::Evaluate,Fitrange_lo,Fitrange_hi,4,"HistoFunctionParam","Evaluate");

 
  Double_t yDataDiff = yInterpolate(DataDiff,half_edge);
  Double_t yMCDiff = yInterpolate(MCDiff,mc_halfedge);
  Double_t DiffScale = yDataDiff/yMCDiff;
  printf(" The predicted y scale is = %f\n",DiffScale);
  Double_t xscale_diff = xscale;
  Double_t xoffset_diff = xoffset;
  DiffFit->SetParNames("yScale_Diff","xScale_Diff","yOffset_Diff","xOffset_Diff");
  DiffFit->SetParameter(0,DiffScale);
  DiffFit->SetParLimits(0,0.5*DiffScale,3*DiffScale);
  DiffFit->FixParameter(2,0.0);
  DiffFit->FixParameter(3,xoffset_diff);
  DiffFit->FixParameter(1,xscale_diff);
  DiffFit->SetRange(Fitrange_lo,Fitrange_hi);
  DiffFit->SetLineColor(kMagenta);
  DiffFit->SetLineWidth(1.0);


  DataDiff->GetXaxis()->SetRangeUser(0.,tail_hi);
  DataDiff->SetMinimum(ymin);
  DataDiff->SetMaximum(1.5);
  cout<<"************************************"<<endl;
  cout<<"Fit to Summed Spectrum"<<endl;
  DataDiff->Fit("DiffFit","rn");

  Double_t yscale_diff = DiffFit->GetParameter(0);

  printf(" Unpolarized fit y scaling = %f\n",yscale_diff);
  printf(" Unpolarized fit x scaling = %f\n",xscale_diff);
  printf(" Unpolarized fit x offset = %f\n",xoffset_diff);
  printf(" Unpolarized fit Chi sqr = %f\n",DiffFit->GetChisquare());
  printf(" Unpolarized fit number of fit points = %f\n",DiffFit->GetNumberFitPoints());

  c_fit->cd(2);

  DataDiff->Draw();
  DiffFit->DrawCopy("same");
  //DiffFit->SetRange(Fitrange_hi,6000);
  //DiffFit->SetLineColor(kRed);
  //DiffFit->DrawCopy("same");


  //---- Fit for Asymmetry (N-P)/(N+P)--------------

  TF1* AsymFit->Clear();
  HistoFunctionParam* AsymFunction = new HistoFunctionParam(MCAsym);
  TF1* AsymFit = new TF1("AsymFit",AsymFunction, &HistoFunctionParam::Evaluate,Fitrange_lo,Fitrange_hi,4,"HistoFunctionParam","Evaluate");

 
  Double_t yDataAsym = yInterpolate(DataAsym,half_edge);
  Double_t yMCAsym = yInterpolate(MCAsym,mc_halfedge);
  Double_t AsymScale = yDataAsym/yMCAsym;
  printf(" The predicted y scale is = %f\n",AsymScale);
  Double_t xscale_asym = xscale;
  Double_t xoffset_asym = xoffset;
  AsymFit->SetParNames("yScale_Asym","xScale_Asym","yOffset_Asym","xOffset_Asym");
  AsymFit->SetParameter(0,AsymScale);
  AsymFit->SetParLimits(0,0.5*AsymScale,3*AsymScale);
  AsymFit->FixParameter(2,0.0);
  AsymFit->FixParameter(3,xoffset_asym);
  AsymFit->FixParameter(1,xscale_asym);
  AsymFit->SetRange(Fitrange_lo,Fitrange_hi);
  AsymFit->SetLineColor(kMagenta);
  AsymFit->SetLineWidth(1.0);


  DataAsym->GetXaxis()->SetRangeUser(0.,tail_hi);
  DataAsym->SetMinimum(asym_min);
  DataAsym->SetMaximum(asym_max);
  cout<<"************************************"<<endl;
  cout<<"Fit to Summed Spectrum"<<endl;
  DataAsym->Fit("AsymFit","rn");

  Double_t yscale_asym = AsymFit->GetParameter(0);
  Double_t fit_error = AsymFit->GetParError(0);

  printf(" Asymmetry fit y scaling = %f\n",yscale_asym);
  printf(" Asymmetry fit x scaling = %f\n",xscale_asym);
  printf(" Asymmetry fit x offset = %f\n",xoffset_asym);
  printf(" Asymmetry fit Chi sqr = %f\n",AsymFit->GetChisquare());
  printf(" Asymmetry fit number of fit points = %f\n",AsymFit->GetNumberFitPoints());

  c_fit->cd(3);

  DataAsym->Draw();
  AsymFit->DrawCopy("same");
  //AsymFit->SetRange(Fitrange_hi,6000);
  //AsymFit->SetLineColor(kRed);
  //AsymFit->DrawCopy("same");

  //-------- create a nice histogram to be saved --------------
  char canvasname[100];
  char savetext[60];
  char unpol_yfit[60];
  char unpol_xoffset[60];
  char unpol_xscale[60];
  char beampol[60];
  char error[60];

  sprintf(canvasname,"Compton Spectrum Fits for run number %d, beam energy 11Gev, and shield thickness %s",RunNumber, shield);
  sprintf(unpol_yfit,"Vertical fit scaling is = %f",yscale);
  sprintf(unpol_xoffset,"Horizontal fit offset is = %.02f",xoffset);
  sprintf(unpol_xscale,"Horizontal fit scaling is = %f",xscale);
  sprintf(beampol,"Beam Polarization = %f",yscale_asym);
  sprintf(savetext,"Fit_%d_11GeV_%smm.pdf",RunNumber, shield);
  sprintf(error,"Polarization error =  %.06f",fit_error);

  TCanvas*c_save = new TCanvas("c_save",canvasname,20,20,900,700);
  c_save->Divide(1,2);
  gStyle->SetOptStat(0);
  c_save->cd(1);

  UnpolarizedFit->SetRange(Fitrange_lo,Fitrange_hi);
  UnpolarizedFit->SetLineColor(kMagenta);
  UnpolarizedFit->SetLineWidth(1.0);
  DataSum->SetTitle("Unpolarized Compton Spectrum");
  DataSum->Draw();
  UnpolarizedFit->DrawCopy("same");
  //UnpolarizedFit->SetRange(Fitrange_hi,25000);
  //UnpolarizedFit->SetLineColor(kRed);
  //UnpolarizedFit->DrawCopy("same");
  UnpolarizedFit->SetLineColor(kMagenta);
  leg_unpol = new TLegend(legxmin,0.55,legxmax,0.85);
  leg_unpol->AddEntry(DataSum,"Unpolarized Compton Data","lep");
  leg_unpol->AddEntry(UnpolarizedFit,"Unpolarized Fit","l");
  leg_unpol->AddEntry((TObject*)0,unpol_yfit,"");
  leg_unpol->AddEntry((TObject*)1,unpol_xscale,"");
  leg_unpol->AddEntry((TObject*)2,unpol_xoffset,"");
  leg_unpol->Draw("SAME");

  cout<<fit_error<<endl;
  c_save->cd(2);
  AsymFit->SetRange(Fitrange_lo,Fitrange_hi);
  AsymFit->SetLineColor(kMagenta);
  AsymFit->SetLineWidth(1.0);
  DataAsym->SetTitle("Compton Asymmetry Spectrum");
  DataAsym->Draw();
  AsymFit->DrawCopy("same");
  //AsymFit->SetRange(Fitrange_hi,25000);
  //AsymFit->SetLineColor(kRed);
  //AsymFit->DrawCopy("same");
  AsymFit->SetLineColor(kMagenta);
  leg_asym = new TLegend(0.15,0.55,0.45,0.85);
  leg_asym->AddEntry(DataAsym,"Compton Asymmetry Data","lep");
  leg_asym->AddEntry(AsymFit,"Asymmetry Fit","l");
  leg_asym->AddEntry((TObject*)0,beampol,"");
  leg_asym->AddEntry((TObject*)1,error,"");
  leg_asym->AddEntry((TObject*)2,unpol_xscale,"");
  leg_asym->AddEntry((TObject*)3,unpol_xoffset,"");
  leg_asym->Draw("SAME");


  c_save->SaveAs(savetext);

}
