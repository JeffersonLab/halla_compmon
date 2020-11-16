void quickTrigDiff()
{
  TChain *T = new TChain("triggerwise");
  T->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3840.root");
  T->Add("/data/cmuwork/rootfiles/tests_2018/compmon_3843.root");
  TH1F *h0 = new TH1F("h0","Neg Hel",50,0,90e3);
  TH1F *h1 = new TH1F("h1","Neg Hel",50,0,90e3);
  TH1F *hd = new TH1F("hd","Pos-Neg",50,0,90e3);
  TH1F *hs = new TH1F("hs","Pos-Neg",50,0,90e3);
  TH1F *ha = new TH1F("ha","Pos-Neg",50,0,90e3);
  TCanvas *canv = new TCanvas("canv","canv",2*500,3*400);
  canv->Divide(2,3);
  canv->cd(1);
  T->Draw("sum>>h1","!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>1700&&helicityState==1");
  canv->cd(2);
  T->Draw("sum>>h0","!sumIsRandom&&sum>-50e3&&cavPowerCalibrated>1700&&helicityState==0");
  canv->cd(3);
  hd->Add(h1,h0,1.0,-1.0);
  hd->Draw();
  canv->cd(4);
  hs->Add(h1,h0,1.0,1.0);
  hs->Draw();
  canv->cd(5);
  ha->Add(hd,1.0);
  ha->Divide(hs);
  ha->Draw("E");
  
}
