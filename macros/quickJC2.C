TFitResultPtr Fit(TH1F *h, const char *msg)
{
  std::cout << "Fit result for " << msg << std::endl;
  Double_t mean = h->GetMean();
  Double_t sig  = h->GetRMS();
  Double_t mult = 1.5;
  Double_t min = mean-mult*sig;
  Double_t max = mean+mult*sig;
  return h->Fit("gaus","S","",min,max);
}
void PrintValErr(Double_t *v, const char *pre)
{
  std::cout << TString::Format("%20s: %+1.10f +/- %1.10f",pre,v[0],v[1]) << std::endl;
}

void quickJC2(Int_t run = 3803)
{
  gStyle->SetOptStat(1110);
  gStyle->SetOptFit(1011);
  TChain *quartetwise = new TChain("quartetwise");
  quartetwise->Add(TString::Format("/data/cmuwork/rootfiles/tests_2018/compmon_%d.root",run));
  quartetwise->SetAlias("MyHel0","0.5*PosHelAcc0/PosHelNSamples0");
  quartetwise->SetAlias("MyHel1","0.5*NegHelAcc1/NegHelNSamples1");
  Double_t xmin = -0.3;
  Double_t xmax =  0.3;
  Double_t dmin = -0.2;
  Double_t dmax =  0.2;
  TCanvas *canv = new TCanvas("canv","canv",500*2,400);
  canv->Divide(2,1);
  TH1F *h0 = new TH1F("h0","Reported Plus Helicity;Acc0/NAcc0[rau]",100,xmin,xmax);
  TH1F *h1 = new TH1F("h1","Reported Minus Helicity;Acc0/Nacc0[rau]",100,xmin,xmax);
  TH1F *hd = new TH1F("hd","Helicity Plus-Minus;[rau]",100,-0.2,0.2);
  h0->SetLineColor(kRed+1);
  h1->SetLineColor(kBlue+1);
  h0->SetLineWidth(2);
  h1->SetLineWidth(2);
  canv->cd(1);
  quartetwise->Draw("MyHel0>>h0");
  canv->cd(2);
  Fit(h0,"Hel0");
  quartetwise->Draw("MyHel1>>h1");
  //h0->Draw();
  //h1->Draw("SAME");
  Fit(h1,"Hel1");
  TCanvas *c1 = new TCanvas("c1","c1",500,400);
  c1->cd(0);
  quartetwise->SetAlias("MyDiff","MyHel0-MyHel1");
  quartetwise->Draw("MyDiff>>hd");
  TCanvas *c2 =new TCanvas("c2","c2",500*2,400*2);
  c2->Divide(2,2);
  TH1F *h[2][2];
  TH1F *hD[2];
  TCut cuts[3];
  cuts[0] = "quartetHelicityPattern==quartetReportedHelicityPattern";
  cuts[1] = "quartetHelicityPattern!=quartetReportedHelicityPattern";
  TString hname;
  Double_t rr[2][2];
  for (Int_t hr = 0; hr < 2; hr++) {
    c2->cd(1+hr*2);
    for (Int_t hs = 0; hs < 2; hs++) {
      hname = TString::Format("h%d%d",hr,hs);
      h[hr][hs] = new TH1F(hname.Data(),TString::Format("Acc0/NAcc0 (Source%sReported);[rau]",hr==0?"#neq":"="),100,xmin,xmax);
      h[hr][hs]->SetLineWidth(2);
      h[hr][hs]->SetLineColor(hs==0?kRed+1:kBlue+1);
      quartetwise->Draw(TString::Format("MyHel%d>>%s",hs,hname.Data()),cuts[hr],(hs==0?"":"SAME"));
    }
    hname = TString::Format("hD%d",hr);
    hD[hr] = new TH1F(hname,"Difference (H^{+}-H^{-});[rau]",100,dmin,dmax);
    hD[hr]->SetLineColor(kGreen-1);
    hD[hr]->SetLineWidth(2);
    c2->cd(2+hr*2);
    quartetwise->Draw(TString::Format("MyDiff>>%s",hname.Data()),cuts[hr]);
    TFitResultPtr r = Fit(hD[hr],hr==0?"R==S":"R!=S");
    rr[hr][0] = r->Value(1);
    rr[hr][1] = r->Error(1);
    if (hr == 0) {
      gPad->Update();
      TPaveStats *st = (TPaveStats*)hD[hr]->FindObject("stats");
      st->SetX1NDC(st->GetX1NDC()-0.5);
      st->SetX2NDC(st->GetX2NDC()-0.5);
      gPad->Update();
    }
  }
  Double_t dd[2];
  Double_t avg[2];
  dd[0] = rr[0][0]+rr[1][0];
  dd[1] = TMath::Sqrt(TMath::Power(rr[0][1],2.0)+TMath::Power(rr[1][1],2));
  avg[0] = 0.5*(rr[0][0]+rr[1][0]);
  avg[1] = dd[1];
  PrintValErr(rr[0],"Diff S==R");
  PrintValErr(rr[1],"Diff S!=R");
  PrintValErr(dd,"Sum of Diffs");
  PrintValErr(avg,"Sum of Diffs");
}
