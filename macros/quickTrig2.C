void quickTrig(Int_t run = 4268)
{
  TChain *T = new TChain("triggerwise");
  T->Add(Form("/data/cmuwork/rootfiles/prex/compmon_%d.root",run));
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
  Int_t entries = T->GetEntries()*0.05;
  for(Int_t entry = 2; entry < entries; entry++ ) {
    T->GetEntry(entry);
    //std::cout << "mpsCoda: " << mpsCoda << ", bcm: " << bcm
    //      << ", sumIsRandom:" << sumIsRandom << ", cavPowerCalibrated: " << cavPowerCalibrated << std::endl;
    if(mpsCoda != lastCoda && bcm > 90 && sumIsRandom==0) {
      if(cavPowerCalibrated>380) {
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
  TH1F *h0 = new TH1F("h0","Neg Hel",100,0,40e3);
  TH1F *h1 = new TH1F("h1","Pos Hel",100,0,40e3);
  TH1F *hb0 = new TH1F("hb0","Neg Hel",100,0,40e3);
  TH1F *hb1 = new TH1F("hb1","Pos Hel",100,0,40e3);
  TH1F *hs0 = new TH1F("hs0","Neg Hel (bksub)",100,0,40e3);
  TH1F *hs1 = new TH1F("hs1","Pos Hel (bksub)",100,0,40e3);
  TH1F *hd = new TH1F("hd","Pos-Neg",100,0,40e3);
  TH1F *hsd = new TH1F("hsd","Pos-Neg (bksub)",100,0,40e3);
  TH1F *hs = new TH1F("hs","Pos-Neg",100,0,40e3);
  TH1F *hss = new TH1F("hss","Pos-Neg (bksub)",100,0,40e3);
  TH1F *ha = new TH1F("ha","Pos-Neg",100,0,40e3);
  TH1F *hsa = new TH1F("hsa","Pos-Neg (bksub)",100,0,40e3);
  TCanvas *canv = new TCanvas("canv","canv",2*500,3*400);
  canv->Divide(4,3);
  canv->cd(1);
  T->Draw("sum>>h1",Form("(1./%f)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>380&&helicityState==1&&bcm>90)",bcmON[1]));
  hb1->SetLineColor(kRed+1);
  hb1->SetLineWidth(2);
  T->Draw("sum>>hb1",Form("(1./%f)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated<40&&helicityState==1&&bcm>90)",bcmOFF[1]),"SAME");
  canv->cd(2);
  T->Draw("sum>>h0",Form("(1./%f)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>380&&helicityState==0&&bcm>90)",bcmON[0]));
  //T->Draw("sum>>h0","(1./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>380&&helicityState==0&&bcm>90)");
  hb0->SetLineColor(kRed+1);
  hb0->SetLineWidth(2);
  T->Draw("sum>>hb0",Form("(1./%f)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated<40&&helicityState==0&&bcm>90)",bcmOFF[0]),"SAME");
  //T->Draw("sum>>hb0","(1./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated<38&&helicityState==0&&bcm>90)","SAME");

  canv->cd(3);
  hs1->Add(h1,hb1,1.0,-1.0);
  //hs0->GetYaxis()->SetMinimum(0.0);
  hs1->SetAxisRange(0.0,0.002,"Y");
  hs1->Draw();
  canv->cd(4);
  hs0->Add(h0,hb0,1.0,-1.0);
  hs0->SetAxisRange(0.0,0.002,"Y");
  hs0->Draw();
  
  canv->cd(5);
  hd->Add(h1,h0,1.0,-1.0);
  hd->Draw();
  canv->cd(6);
  hs->Add(h1,h0,1.0,1.0);
  hs->Draw();
  canv->cd(7);
  hsd->Add(hs1,hs0,1.0,-1.0);
  hsd->SetAxisRange(-0.0005,0.001,"Y");
  hsd->Draw();

  canv->cd(8);
  hss->Add(hs1,hs0,1.0,1.0);
  hss->SetAxisRange(0.0,0.004,"Y");
  hss->Draw();



  canv->cd(9);
  ha->Add(hd,1.0);
  ha->Divide(hs);
  ha->SetAxisRange(-0.05,0.05,"Y");
  ha->Draw("E");

  canv->cd(11);
  hsa->Add(hsd,1.0);
  hsa->Divide(hss);
  hsa->SetAxisRange(-0.05,0.05,"Y");
  hsa->Draw("E");
  
}
