// Based on fitPulser.C, from Brian Quinn
/* Comments: 
Here's a routine to fit Alexa's gr TGraphErrors.  You'll have to change 
plotPulser to make gr a global (just moving the declaration outside the 
routine). Then you run plotPulser as usual, .L fitPulser and 
fitPulser(n,m,s)  where n is number of adjustable parameters (in addition 
to delta), m is number of data points, and s is slope measured for 
dark-variable run (see below).  So I've been using fitPulser(2,150,0.).

   If you give it n=0, it does a fit with only delta varying (so a 
horizontal line fit to the data) but it also figures you would only do 
that for a dark delta run (which goes bananas, as it should, if you tell 
it to vary parameters based on delta=0) so it also does an LSF fit to the 
data and tells you the slope, s.  (I think it calls it LSF:m for slope, 
but it's the s value you enter for subsequent fits if you want to correct 
for the measured dark-Delta slope.)

   I had it working yesterday, but it would only work once and then you'd 
have to re-run plotPulser or fitPulser would give nonsense.  It turned out 
it was re-scaling the data in gr (I didn't even know it could do that) and 
then re-scaling it AGAIN by the same factor and making the data get very 
tiny.  It no longer messes with the data in gr, it copies it to separate 
arrays and messes with them. */


#include <iostream>
#include <fstream>
#include "TMath.h"

Double_t Delta; 
int params;

const static int nmax = 200;

TGraphErrors *grsc;
TGraph *gr2 = gROOT->FindObject("gr");


/* July 2015, changed to one additional parameter, rather than treating 
delta as a constant which is (somehow??) externally adjusted.  
The array par now contains the params-1 coefficients followed by delta.
eg for 4th order fit:
par(0)  coeff of x^2
par(1)  coeff of x^3
par(2)  coeff of x^4
par(3)  delta
*/


Double_t fcn(Double_t *x, Double_t *par){

  Double_t value=par[params-1];
  for(int i=2;i<params+1;i++){
    value += par[i-2]*pow(x[0]+Delta,i) - par[i-2]*pow(x[0],i);
  }
  return value;
}


Double_t draw_fd(Double_t del, Double_t alpha, Double_t beta, Double_t gamma){
/*
Calling this routine superimposes a finite-difference function using the 
parameters of your choice  del=delta, alpha,beta,gamma are coefficients of 
2nd,3rd and 4th order non linearity terms.
New kludge: gamma specifies the 3rd coeff, unless it's >100, then it specifies line color
 */
  int my_col=4; // blue by default (gamma-100 if gamma>100)
  if(gamma>100.){my_col=gamma-100.; gamma=0;}
  TF1 *my_fd = new TF1("my_fd",
" [0]+[1]*(pow(x+[0],2.)-x*x)+[2]*(pow(x+[0],3.)-pow(x,3))+[3]*(pow(x+[0],4.)-pow(x,4)) ", 0., 1.1);
  my_fd->SetParameter(0,del);
  my_fd->SetParameter(1,alpha);
  my_fd->SetParameter(2,beta);
  my_fd->SetParameter(3,gamma);
  my_fd->SetLineColor(my_col);
  my_fd->Draw("same");
  return 0;
}


Double_t responsefcn(Double_t x, Double_t *par){

  Double_t value=x;
  for(int i=2;i<params+1;i++){
    value += par[i-2]*pow(x,i); 
  }
  return value;
}


Double_t deriv(Double_t x, Double_t *par){

  Double_t value=1;
  for(int i=2;i<params+1;i++){
    value += i*par[i-2]*pow(x,i-1); 
  }
  if(abs(value)<1e-6) value=1;
  return value;
}


void fitPulser(int numparams, int range, Double_t cross_LSF){

// Fits 'gr' filled by plotPulser.C

/* If numparams=0 (only delta being adjusted) it does an LSF which is 
   useful for dark-delta runs.
   added cross_LSF argument which is actual slope (with sign) of LSF fit to
   dark delta cross-talk measurement.
*/

  char parameter[30];

  if(range>nmax){
    printf("Please increase array dimensions.\n\n");
      return();
  }

  Double_t ys[nmax], xs[nmax], ers[nmax];
  Double_t finiteDifferencesc[nmax],Responsesc[nmax],Yerrsc[nmax];

  Double_t *finiteDifference = gr->GetY();
  Double_t *Response = gr->GetX();
  Double_t *Yerr = gr->GetEY();

  for(int i=0;i<nmax;i++){
    finiteDifferencesc[i]=finiteDifference[i];
    Responsesc[i]=Response[i];
    Yerrsc[i]=Yerr[i];
  }


  if(Responsesc[0]<0){
    double offset=-Responsesc[0]+(Responsesc[range-1]-Responsesc[0])/100000.;
    printf("***********************************************************\n");
    printf("***********************************************************\n");
    printf("Tweaking data by %f to avoid negative intensity\n",offset);
    printf("***********************************************************\n");
    printf("***********************************************************\n");
    for(int i=0;i<range;i++){
	Responsesc[i]+=offset;
      }
  }

  Double_t Asum=0, Bsum=0,Csum=0,Dsum=0,Esum=0,det,xi,yi,w;

  for(int ipt=0;ipt<range;ipt++){
    finiteDifferencesc[ipt]=finiteDifferencesc[ipt]-cross_LSF*Responsesc[ipt];
    if(numparams==0){
      xi=Responsesc[ipt];
      yi=finiteDifferencesc[ipt];
      w=1/Yerrsc[ipt]/Yerrsc[ipt];
      Asum+=xi*xi*w;
      Bsum+=xi*w;
      Csum+=w;
      Dsum+=xi*yi*w;
      Esum+=yi*w;
    }
  }

  if(numparams==0){
    det=Asum*Csum-Bsum*Bsum;
    Double_t LSFm=(Csum*Dsum-Bsum*Esum)/det;
    Double_t LSFb=(Asum*Esum-Bsum*Dsum)/det;
    Double_t Sigmam=sqrt(Csum/det);
    Double_t Sigmab=sqrt(Asum/det);
    Double_t chi2lsf=0;
    for(int ipt=0;ipt<range;ipt++){
      chi2lsf+=pow((LSFm*Responsesc[ipt]+LSFb)-finiteDifferencesc[ipt],2)/
	pow(Yerrsc[ipt],2);
    }
    //  printf("LSF: M=%f+/-%f, b=%f+/-%f, chi^2=%f,  chi^2/nu=%f\n",
    //	 LSFm,Sigmam,LSFb,Sigmab,chi2,chi2/(range-2));
  }


  const int n = range;
  const int numpar = numparams+1;
  params = numparams+1;
  Double_t yerr[n] = 0;
  Double_t xerr[n] = 0;
  // Double_t ped = hsumP1->GetMean(1);
  // Double_t low = Responsesc[0];
  // Double_t high = Responsesc[n-2];

  Double_t high = Responsesc[0];
  Double_t DeltaInit = finiteDifferencesc[0]/high;
  Double_t LastDelta = 0;
  Delta = DeltaInit;

  Double_t newy0 = -1;
  Double_t y0 = 0;
  Double_t par[numpar]=0;
  Double_t chi2 = 0;
  Double_t oldchi2 = -1;

  Double_t eDeposited[n] = 0;
  Double_t eDepLast[n] = 0;
  for(int j=0;j<n;j++){
    eDeposited[j] = Responsesc[j]/high;
    Responsesc[j] = Responsesc[j]/high;
    finiteDifferencesc[j] = finiteDifferencesc[j]/high;
    //Double_t *Xerr= graphErrors->GetErrorX(j);
    yerr[j] = Yerrsc[j]/high;
    xerr[j] = 0;
  }  

  if(numparams==0){
    printf("LSF: M=%f+/-%f, b=%f+/-%f, chi^2=%f,  chi^2/nu=%f\n",
	 LSFm,Sigmam,LSFb/high,Sigmab/high,chi2lsf,chi2lsf/(range-2));
  }


  TF1 *fitFcn = new TF1("fitFcn", fcn, eDeposited[0], eDeposited[n-1], params);

  while(abs(chi2-oldchi2)>1e-6 || xchanged){
    int xchanged = 0;
    /*    cerr<<"Delta: "<<Delta<<endl;
    if(abs(LastDelta-Delta)>1e-6)xchanged = 1;
    */
    cerr<<"Delta: "<<par[params-1]<<endl;
    for(int j=0;j<n;j++){
      if(abs(eDepLast[j]-eDeposited[j])>1e-6) xchanged = 1;
      eDepLast[j] = eDeposited[j];
    }

    //TCanvas *can = new TCanvas();
    oldchi2 = chi2;
    grsc = new TGraphErrors(n, eDeposited, finiteDifferencesc, xerr, yerr);
    grsc->SetMarkerStyle(20);
    grsc->GetXaxis()->SetTitle("Scaled Light Deposited");
    grsc->GetYaxis()->SetTitle("Scaled Finite Difference");
    grsc->GetXaxis()->SetLabelSize(0.04);
    grsc->GetYaxis()->SetLabelSize(0.03);
    grsc->SetTitle("Graph to Fit");
    grsc->Draw("ap");
    //grsc->Fit("fitFcn","ERQ");    
    grsc->Fit("fitFcn");


    if(numparams==0){
      TLine *LSF = new TLine(0.,LSFb/high,1.,LSFm+LSFb/high);
      LSF->Draw();
    }

    fitFcn->GetParameters(par);
    chi2 = fitFcn->GetChisquare();

    for(int j=0;j<n;j++){
      y0=0.;
      newy0=-1.;
      while(abs(y0-newy0)>1e-6){
	y0 = responsefcn(eDeposited[j],par);
	eDeposited[j] = (Responsesc[j]-y0)/deriv(eDeposited[j],par) + eDeposited[j];
	newy0 = responsefcn(eDeposited[j],par);
	//cerr<<"j: "<<j<<" newy0 "<<newy0[j]<<" y0: "<<y0<<" eDeposited: "<<eDeposited[j]<<endl;
      }
    }

  }
  for(int i=0;i<n;i++){
    eDeposited[i]=high*eDeposited[i];
    Responsesc[i]=high*Responsesc[i];
  }

  //TCanvas *can = new TCanvas();
  TGraph *gr2 = new TGraph(n, eDeposited, Responsesc);
  gr2->GetXaxis()->SetTitle("Light Deposited");
  gr2->GetYaxis()->SetTitle("Response");
  gr2->GetXaxis()->SetLabelSize(0.04);
  gr2->GetYaxis()->SetLabelSize(0.03);
  gr2->SetTitle("Response Function");
  TLine *line = new TLine(0,0,high,high);

  //TCanvas *can = new TCanvas();
  TGraph *grflip = new TGraph(n, Responsesc, eDeposited);
  grflip->GetXaxis()->SetTitle("Response");
  grflip->GetYaxis()->SetTitle("Light Deposited");
  grflip->GetXaxis()->SetLabelSize(0.04);
  grflip->GetYaxis()->SetLabelSize(0.03);
  grflip->SetTitle("Response Function, Flipped");

  cout<<"Chi square: "<<chi2/(n-params)<<endl;
  cout<<"Scale factor: "<<high<<endl;

  printf("\nDone.\n\n");

}
