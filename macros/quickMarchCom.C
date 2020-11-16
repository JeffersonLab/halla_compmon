#include <TCut.h>
const bool kDoHelicityCut = kFALSE;
const bool kBCMNormalize = kFALSE;
const bool kBkSubtract = kTRUE;

const Int_t kAccumNum = 0+4;

double gONMinH = 0;
double gONMaxH = 0;
double gONMinD = -2.;
double gONMaxD = 2.;
double gONMinS = 11.0;
double gONMaxS = 20.0;
double gONMinA = -0.1;
double gONMaxA =  0.1;
double gONFitMinD = -0.3;
double gONFitMaxD = 1.0;
double gONFitMinA = -0.03;
double gONFitMaxA =  0.075;

double gOFFMinH = -0.5;
double gOFFMaxH =  2.0;
double gOFFMinD = -0.8;
double gOFFMaxD =  0.8;
double gOFFMinS =  0.0;
double gOFFMaxS =  2.2;


TChain *T;
TChain *triggerwise;
void quickMarchCom()
{
  T = new TChain("quartetwise");
  //T->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3839.root");
  //T->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3840.root");
  T->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3843.root");
  T->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3846.root");
  triggerwise = new TChain("triggerwise");
  triggerwise->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3843.root");
  triggerwise->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3846.root");

  if(kAccumNum==0) {
    if(kBCMNormalize) {
      T->SetAlias("HelP",Form("(PosHelAcc%d/PosHelNSamples%d)/bcm",kAccumNum,kAccumNum));
      T->SetAlias("HelN",Form("(NegHelAcc%d/NegHelNSamples%d)/bcm",kAccumNum,kAccumNum));
    } else {
      T->SetAlias("HelP",Form("(PosHelAcc%d/PosHelNSamples%d)",kAccumNum,kAccumNum));
      T->SetAlias("HelN",Form("(NegHelAcc%d/NegHelNSamples%d)",kAccumNum,kAccumNum));
    }
  } else {
    if(kBCMNormalize) {
      T->SetAlias("HelP",Form("(PosHelAcc%d)/bcm",kAccumNum));
      T->SetAlias("HelN",Form("(NegHelAcc%d)/bcm",kAccumNum));
    } else {
      T->SetAlias("HelP",Form("(PosHelAcc%d)",kAccumNum));
      T->SetAlias("HelN",Form("(NegHelAcc%d)",kAccumNum));
    }
  }

  TCut cutLON = "cavPowerCalibrated>1500";
  TCut cutLOFF = "cavPowerCalibrated<5";
  TCut helCut = "helicityState==helicityStateReported";
  TCut cuts = "bcm>10" ;//&& helCut;

  if(kAccumNum == 4) {
    gONMinD = -20*0;
    gONMaxD =  20.0*0;
    gONMinS = 0.;
    gONMaxS = 0.;

    gOFFMinH = 0;
    gOFFMaxH = 0;
    gOFFMinD = 0;
    gOFFMaxD = 0;
    gOFFMinS = 0.;
    gOFFMaxS = 0.;

    gONFitMinD = -2e6;
    gONFitMaxD =  12e6;
    gONFitMinA = -0.01;
    gONFitMaxA =  0.06;

  }

  // First plot the Laser OFF to get an idea of background
  TCanvas *canvOFF = new TCanvas("canvOFF","canvOFF",2*500,2*400);
  canvOFF->Divide(2,2);
  canvOFF->cd(1);
  TH1F *h1OFF = new TH1F("h1OFF","Laser OFF Pos Hel",100,gOFFMinH,gOFFMaxH);
  TH1F *h0OFF = new TH1F("h0OFF","Laser OFF Neg Hel",100,gOFFMinH,gOFFMaxH);
  TH1F *hDiffOFF = new TH1F("hDiffOFF","Laser OFF (Pos-Neg)",100,gOFFMinD,gOFFMaxD);
  TH1F *hSumOFF = new TH1F("hSumOFF","Laser OFF (Pos+Neg)",100,gOFFMinS,gOFFMaxS);
  T->Draw("HelP>>h1OFF",cuts&&cutLOFF);
  //h1OFF->Fit("gaus","","",0.3,1.03);
  canvOFF->cd(2);
  T->Draw("HelN>>h0OFF",cuts&&cutLOFF);
  //h0OFF->Fit("gaus","","",0.1,.8);
  canvOFF->cd(3);
  //T->Draw("(PosHelAcc0/PosHelNSamples0)-(NegHelAcc0/NegHelNSamples0)>>hDiffOFF","cavPowerCalibrated<5&&bcm>10&&helicityStateReported==helicityState");
  T->Draw("HelP-HelN>>hDiffOFF",cuts&&cutLOFF);
  canvOFF->cd(4);
  //T->Draw("(PosHelAcc0/PosHelNSamples0)+(NegHelAcc0/NegHelNSamples0)>>hSumOFF","cavPowerCalibrated<5&&bcm>10&&helicityState==helicityStateReported");
  T->Draw("HelP+HelN>>hSumOFF",cuts&&cutLOFF);

  Double_t bk0 = h0OFF->GetMean();
  Double_t bk1 = h1OFF->GetMean();


  //T->SetAlias("LONPos","(PosHelAcc0/PosHelNSamples0)-6.81167e-01");
  //T->SetAlias("LONNeg","(NegHelAcc0/NegHelNSamples0)-3.98470e-01");
  //T->SetAlias("LONPos","(PosHelAcc0/PosHelNSamples0)");
  //T->SetAlias("LONNeg","(NegHelAcc0/NegHelNSamples0)");
  TCanvas *canv = new TCanvas("canv","canv",2*500,3*400);
  canv->Divide(2,3);
  canv->cd(1);
  //T->SetAlias("LONPos","(PosHelAcc0/PosHelNSamples0)-6.81167e-01");
  //T->SetAlias("LONNeg","(NegHelAcc0/NegHelNSamples0)-3.98470e-01");
  if(kBkSubtract) {
    T->SetAlias("LONPos",Form("HelP-%g",bk1));
    T->SetAlias("LONNeg",Form("HelN-%g",bk0));
  } else {
    T->SetAlias("LONPos","HelP");
    T->SetAlias("LONNeg","HelN");
  }
    //T->SetAlias("LONPos","(PosHelAcc0/PosHelNSamples0)");
  //T->SetAlias("LONNeg","(NegHelAcc0/NegHelNSamples0)");
  TH1F *h1ON = new TH1F("h1ON","Laser ON Pos Hel",100,gONMinH,gONMaxH);
  TH1F *h0ON = new TH1F("h0ON","Laser ON Neg Hel",100,gONMinH,gONMaxH);
  TH1F *hDiffON = new TH1F("hDiffON","Laser ON (Pos-Neg)",100,gONMinD,gONMaxD);
  TH1F *hSumON = new TH1F("hSumON","Laser ON (Pos+Neg)",100,gONMinS,gONMaxS);
  T->Draw("LONPos>>h1ON",cuts&&cutLON);
  canv->cd(2);
  T->Draw("LONNeg>>h0ON",cuts&&cutLON);
  canv->cd(3);
  T->Draw("(LONPos)-(LONNeg)>>hDiffON",cuts&&cutLON);
  //hDiffON->Fit("gaus","","",0.,1.4);
  hDiffON->Fit("gaus","","",gONFitMinD,gONFitMaxD);
  canv->cd(4);
  T->Draw("(LONPos)+(LONNeg)>>hSumON",cuts&&cutLON);
  canv->cd(5);
  TH1F *hAsymON = new TH1F("hAsymON","Asym",100,gONMinA,gONMaxA);
  T->Draw("((LONPos)-(LONNeg))/((LONPos)+(LONNeg))>>hAsymON",cuts&&cutLON);
  hAsymON->Fit("gaus","","",gONFitMinA,gONFitMaxA);

  //canvOFF->cd(5);
  //T->Draw("((PosHelAcc0/PosHelNSamples0)-(NegHelAcc0/NegHelNSamples0))/((PosHelAcc0/PosHelNSamples0)+(NegHelAcc0/NegHelNSamples0))>>hAsymOFF","cavPowerCalibrated<5&&bcm>10&&helicityStateReported==helicityState");

 TCanvas *canvTrig = new TCanvas("canvtrig","canvtrig",600,500);
 triggerwise->Draw("sum","sum>-50e3&&sum<120e3&&!sumIsRandom&&cavPowerCalibrated>1500");
 T->Draw("cavPowerCalibrated",cuts&&cutLON);
}

