TCut cLaser[2] = {"laserState==3","laserState==1"};
TCut cCuts = "sumIsRandom==0&&abs(sumPre-3790)<10&&beamState==1";
TCut cHelicity[2] = {"helicityState==0","helicityState==1"};
const Int_t kNbins = 50;
const char* kDraw="sum-500";

void quickTrig3(Int_t run = 4294)
{
  TChain *T = new TChain("triggerwise");
  T->Add(Form("/data/cmuwork/rootfiles/prex/compmon_%d.root",run));
  T->Add(Form("/data/cmuwork/rootfiles/prex/compmon_%d.root",run-1));
  T->Add(Form("/data/cmuwork/rootfiles/prex/compmon_%d.root",run-2));
/*
  Float_t bcm;
  Float_t bcmON[2] = {0.,0.};
  Float_t bcmOFF[2] = {0.0, 0.0};
  Int_t mpsCoda;
  Int_t lastCoda = -1;
  Bool_t sumIsRandom;
  Int_t helicityState;
  Float_t cavPowerCalibrated;
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("mpsCoda",1);
  T->SetBranchStatus("bcm",1);
  T->SetBranchStatus("sumIsRandom",1);
  T->SetBranchStatus("helicityState",1);
  T->SetBranchStatus("cavPowerCalibrated",1);
  T->SetBranchAddress("mpsCoda",&mpsCoda);
  T->SetBranchAddress("bcm",&bcm);
  T->SetBranchAddress("sumIsRandom",&sumIsRandom);
  T->SetBranchAddress("cavPowerCalibrated",&cavPowerCalibrated);
  T->SetBranchAddress("helicityState",&helicityState);
  Int_t entries = T->GetEntries();
  for(Int_t entry = 2; entry < entries; entry++ ) {
    T->GetEntry(entry);
    //std::cout << "mpsCoda: " << mpsCoda << ", bcm: " << bcm
    //      << ", sumIsRandom:" << sumIsRandom << ", cavPowerCalibrated: " << cavPowerCalibrated << std::endl;
    if(mpsCoda != lastCoda && bcm > 25 && sumIsRandom==0) {
      if(cavPowerCalibrated>280) {
       for(Int_t h = 0; h < 2; h++) {
        if(helicityState==h) {
          bcmON[h] += bcm;
        }
       }
      } else if (cavPowerCalibrated<40) {
        for(Int_t h = 0; h < 2; h++) {
          if(helicityState==h) {
            bcmOFF[h] += bcm;
          }
        }
      }
    }
    lastCoda = mpsCoda;
  }
  std::cout
    << "bcmON[0]: " << bcmON[0] << std::endl
    << "bcmON[1]: " << bcmON[1] << std::endl
    << "bcmOFF[0]: " << bcmOFF[0] << std::endl
    << "bcmOFF[1]: " << bcmOFF[1] << std::endl;
  T->SetBranchStatus("*",1);
*/
  TH1F *hh[2];
  TH1F *hb[2];
  TH1F *hs[2];
  TH1F *hhn[2];
  TH1F *hbn[2];
  
  hh[0] = new TH1F("hh0","Neg Hel",kNbins,-1e3,40e3);
  hh[1] = new TH1F("hh1","Pos Hel",kNbins,-1e3,40e3);
  hb[0] = new TH1F("hb0","Neg Hel",kNbins,-1e3,40e3);
  hb[1] = new TH1F("hb1","Pos Hel",kNbins,-1e3,40e3);
  hhn[0] = new TH1F("hhn0","Neg Hel",kNbins,40e3,50e4);
  hhn[1] = new TH1F("hhn1","Pos Hel",kNbins,40e3,50e4);
  hbn[0] = new TH1F("hbn0","Neg Hel",kNbins,40e3,50e4);
  hbn[1] = new TH1F("hbn1","Pos Hel",kNbins,40e3,50e4);

  hs[0] = new TH1F("hs0","Neg Hel (bksub)",kNbins,-1e3,40e3);
  hs[1] = new TH1F("hs1","Pos Hel (bksub)",kNbins,-1e3,40e3);
  TH1F *hhd = new TH1F("hhd","Pos-Neg",kNbins,-1e3,40e3);
  TH1F *hsd = new TH1F("hsd","Pos-Neg (bksub)",kNbins,-1e3,40e3);
  TH1F *hhs = new TH1F("hhs","Pos-Neg",kNbins,-1e3,40e3);
  TH1F *hss = new TH1F("hss","Pos-Neg (bksub)",kNbins,-1e3,40e3);
  TH1F *hha = new TH1F("hha","Pos-Neg",kNbins,-1e3,40e3);
  TH1F *hsa = new TH1F("hsa","Pos-Neg (bksub)",kNbins,-1e3,40e3);
  TCanvas *canv = new TCanvas("canv","canv",2*500,3*400);
  for(Int_t h = 0; h < 2; h++) {
    hb[h]->SetLineColor(kRed+1);
    hb[h]->SetFillColor(kRed+1);
    hbn[h]->SetLineColor(kRed+1);
    hs[h]->SetLineColor(kMagenta+1);
    hs[h]->SetFillColor(kMagenta+1);
  }
  canv->Divide(4,2);
  Double_t norm[2];
  for(Int_t h = 0; h < 2; h++) {
   canv->cd(1+h);
   for(Int_t l = 1; l >= 0; l--) {
     T->Draw(Form("%s>>h%sn%d",kDraw,l==0?"b":"h",h),cCuts&&cHelicity[h]&&cLaser[l],l==1?"":"SAME");
   }
   norm[h]=hhn[h]->Integral()/hbn[h]->Integral();
   for(Int_t l = 1; l >= 0; l--) {
     T->Draw(Form("%s>>h%s%d",kDraw,l==0?"b":"h",h),cCuts&&cHelicity[h]&&cLaser[l],l==1?"":"SAME");
   }
   hb[h]->SaveAs(Form("/home/compton/cornejo/plots_prex_commissioning/trigsums_%d_bkg_h%d.C",run,h));
   hb[h]->Scale(norm[h]);
   hb[h]->Draw("HSAME");
   
   canv->cd(3+h);
   hs[h]->Add(hh[h],hb[h],1.0,-1.0);
   hh[h]->Draw();
   hs[h]->Draw("HSAME");
  }

  
  canv->cd(5);
  hhd->Add(hh[1],hh[0],1.0,-1.0);
  hhd->Draw();
  canv->cd(6);
  hhs->Add(hh[1],hh[0],1.0,1.0);
  hhs->Draw();

  canv->cd(7);
  hsd->Add(hs[1],hs[0],1.0,-1.0);
  //hsd->SetAxisRange(-0.0005,0.001,"Y");
  hsd->Draw("H");

  canv->cd(8);
  hss->Add(hs[1],hs[0],1.0,1.0);
  //hss->SetAxisRange(0.0,0.004,"Y");
  hss->Draw();


  TCanvas *canv2 = new TCanvas("canv2","canv",2*600,500*2);
  canv2->Divide(1,2);

  canv2->cd(1);
  hha->Add(hhd,1.0);
  hha->Divide(hhs);
  hha->SetAxisRange(-0.05,0.05,"Y");
  hha->Draw("E");

  canv2->cd(2);
  hsa->Add(hsd,1.0);
  hsa->Divide(hss);
  hsa->SetAxisRange(-0.05,0.05,"Y");
  hsa->Draw("E");
}
