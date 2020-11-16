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
  canv->Divide(4,2);
  canv->cd(1);
  T->Draw("sum>>h1",Form("(1./%f)*(520./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>300&&helicityState==1&&bcm>25)",bcmON[1]));
  hb1->SetLineColor(kRed+1);
  hb1->SetLineWidth(2);
  T->Draw("sum>>hb1",Form("(1./%f)*(520./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated<40&&helicityState==1&&bcm>25)",bcmOFF[1]),"SAME");
  gPad->SetLogy(1);
  h1->Scale(hb1->Integral()/h1->Integral());
  h1->Draw();
  hb1->Draw("SAME");
  TF1 *f1 = new TF1("f1","[0]*exp(-[1]*(x)-[2])",4600,15e3);//,3,4600,15e3);
  //hb1->Fit(f1,"","",4600,10e3);
  //hb1->SaveAs("test_out_jc2.C");
  canv->cd(2);
  T->Draw("sum>>h0",Form("(1./%f)*(520./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>300&&helicityState==0&&bcm>25)",bcmON[0]));
  //T->Draw("sum>>h0","(1./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>300&&helicityState==0&&bcm>25)");
  hb0->SetLineColor(kRed+1);
  hb0->SetLineWidth(2);
  T->Draw("sum>>hb0",Form("(1./%f)*(520./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated<40&&helicityState==0&&bcm>25)",bcmOFF[0]),"SAME");
  gPad->SetLogy(1);
  //T->Draw("sum>>hb0","(1./scaler_run7)*(!sumIsRandom&&sum>-50e3&&cavPowerCalibrated<38&&helicityState==0&&bcm>25)","SAME");

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


  TCanvas *canv2 = new TCanvas("canv2","canv",2*600,500);
  canv2->Divide(2,1);

  canv2->cd(0);
  ha->Add(hd,1.0);
  ha->Divide(hs);
  ha->SetAxisRange(-0.05,0.05,"Y");
  ha->Draw("E");

  /*canv2->cd(2);
  hsa->Add(hsd,1.0);
  hsa->Divide(hss);
  hsa->SetAxisRange(-0.05,0.05,"Y");
  hsa->Draw("E");
 */ 
}
