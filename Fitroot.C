

void Fitroot(Int_t run = 4294)

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
  //std::cout
    // << "bcmON[0]: " << bcmON[0] << std::endl
    //  << "bcmON[1]: " << bcmON[1] << std::endl
    // << "bcmOFF[0]: " << bcmOFF[0] << std::endl
    // << "bcmOFF[1]: " << bcmOFF[1] << std::endl;
  T->SetBranchStatus("*",1);
  ;
  TH1F *h1 = new TH1F("h1","spectrum",100,0,60e3);
  TH1F *h2 = new TH1F("h2","Background",100,0,34);
  TH1F *h3 = new TH1F("h3","Background",100,0,34);
  TH1F *h4 = new TH1F("h4","Background",100,0,34);
  TCanvas *canv = new TCanvas("canv","canv",1000,1000);
  canv->Divide(1,3);
  canv->cd(1);
  T->Draw("sum>>h1",Form("(sumIsRandom!=1&&sum>-100e3&&sum<75e3&&cavPowerCalibrated>380 && bcm>25)"));

  canv->cd(2);
  T->Draw("sum/1000>>h2",Form("(sumIsRandom!=1&&sum>-100e3&&cavPowerCalibrated<40&&bcm>25 )"));

  canv->cd(3);
  TH1F *htemp__1 = new TH1F("htemp__1","simulation",100,0,34);
   htemp__1->SetBinContent(1,7144);
   htemp__1->SetBinContent(2,7041);
   htemp__1->SetBinContent(3,8239);
   htemp__1->SetBinContent(4,10846);
   htemp__1->SetBinContent(5,13594);
   htemp__1->SetBinContent(6,16252);
   htemp__1->SetBinContent(7,20274);
   htemp__1->SetBinContent(8,28907);
   htemp__1->SetBinContent(9,42744);
   htemp__1->SetBinContent(10,67845);
   htemp__1->SetBinContent(11,65883);
   htemp__1->SetBinContent(12,64608);
   htemp__1->SetBinContent(13,62746);
   htemp__1->SetBinContent(14,62194);
   htemp__1->SetBinContent(15,60675);
   htemp__1->SetBinContent(16,59621);
   htemp__1->SetBinContent(17,58499);
   htemp__1->SetBinContent(18,57532);
   htemp__1->SetBinContent(19,56798);
   htemp__1->SetBinContent(20,55713);
   htemp__1->SetBinContent(21,55256);
   htemp__1->SetBinContent(22,54315);
   htemp__1->SetBinContent(23,53145);
   htemp__1->SetBinContent(24,52609);
   htemp__1->SetBinContent(25,52354);
   htemp__1->SetBinContent(26,51279);
   htemp__1->SetBinContent(27,50699);
   htemp__1->SetBinContent(28,50334);
   htemp__1->SetBinContent(29,49703);
   htemp__1->SetBinContent(30,49174);
   htemp__1->SetBinContent(31,48599);
   htemp__1->SetBinContent(32,48344);
   htemp__1->SetBinContent(33,48074);
   htemp__1->SetBinContent(34,47472);
   htemp__1->SetBinContent(35,47315);
   htemp__1->SetBinContent(36,47063);
   htemp__1->SetBinContent(37,47072);
   htemp__1->SetBinContent(38,46701);
   htemp__1->SetBinContent(39,46861);
   htemp__1->SetBinContent(40,46625);
   htemp__1->SetBinContent(41,46540);
   htemp__1->SetBinContent(42,46543);
   htemp__1->SetBinContent(43,46505);
   htemp__1->SetBinContent(44,46868);
   htemp__1->SetBinContent(45,47292);
   htemp__1->SetBinContent(46,46633);
   htemp__1->SetBinContent(47,47362);
   htemp__1->SetBinContent(48,47825);
   htemp__1->SetBinContent(49,48067);
   htemp__1->SetBinContent(50,48041);
   htemp__1->SetBinContent(51,48749);
   htemp__1->SetBinContent(52,48716);
   htemp__1->SetBinContent(53,49419);
   htemp__1->SetBinContent(54,50091);
   htemp__1->SetBinContent(55,50215);
   htemp__1->SetBinContent(56,50711);
   htemp__1->SetBinContent(57,51851);
   htemp__1->SetBinContent(58,52298);
   htemp__1->SetBinContent(59,52745);
   htemp__1->SetBinContent(60,53318);
   htemp__1->SetBinContent(61,53236);
   htemp__1->SetBinContent(62,54149);
   htemp__1->SetBinContent(63,54492);
   htemp__1->SetBinContent(64,55006);
   htemp__1->SetBinContent(65,55618);
   htemp__1->SetBinContent(66,55915);
   htemp__1->SetBinContent(67,56556);
   htemp__1->SetBinContent(68,56573);
   htemp__1->SetBinContent(69,56853);
   htemp__1->SetBinContent(70,56791);
   htemp__1->SetBinContent(71,57375);
   htemp__1->SetBinContent(72,57147);
   htemp__1->SetBinContent(73,56709);
   htemp__1->SetBinContent(74,56988);
   htemp__1->SetBinContent(75,56340);
   htemp__1->SetBinContent(76,55607);
   htemp__1->SetBinContent(77,54969);
   htemp__1->SetBinContent(78,53916);
   htemp__1->SetBinContent(79,52650);
   htemp__1->SetBinContent(80,51122);
   htemp__1->SetBinContent(81,49429);
   htemp__1->SetBinContent(82,47729);
   htemp__1->SetBinContent(83,44970);
   htemp__1->SetBinContent(84,42311);
   htemp__1->SetBinContent(85,38763);
   htemp__1->SetBinContent(86,35097);
   htemp__1->SetBinContent(87,31106);
   htemp__1->SetBinContent(88,26404);
   htemp__1->SetBinContent(89,20821);
   htemp__1->SetBinContent(90,14984);
   htemp__1->SetBinContent(91,8465);
   htemp__1->SetBinContent(92,3813);
   htemp__1->SetEntries(4271842);
   htemp__1->SetDirectory(0);
   
  
   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   htemp__1->SetLineColor(ci);
   htemp__1->GetXaxis()->SetTitle("GSOCrystalPhysical_eDep");
   htemp__1->GetXaxis()->SetRange(1,100);
   htemp__1->GetXaxis()->SetLabelFont(42);
   htemp__1->GetXaxis()->SetLabelSize(0.035);
   htemp__1->GetXaxis()->SetTitleSize(0.035);
   htemp__1->GetXaxis()->SetTitleFont(42);
   htemp__1->GetYaxis()->SetLabelFont(42);
   htemp__1->GetYaxis()->SetLabelSize(0.035);
   htemp__1->GetYaxis()->SetTitleSize(0.035);
   htemp__1->GetYaxis()->SetTitleFont(42);
   htemp__1->GetZaxis()->SetLabelFont(42);
   htemp__1->GetZaxis()->SetLabelSize(0.035);
   htemp__1->GetZaxis()->SetTitleSize(0.035);
   htemp__1->GetZaxis()->SetTitleFont(42);
    htemp__1->Draw("");
   
   TPaveText *pt = new TPaveText(0.15,0.9342405,0.85,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   pt->Draw();

   T->Draw("sum/1000>>h3",Form("(sumIsRandom!=1&&sum>-100e3&&sum<75e3&&cavPowerCalibrated>380 && bcm>25)"),"same");
  h3->SetLineColor(kRed+1);
  h3->Scale(3);
  h3->DrawCopy("same");
  // htemp__1->Draw("same");


 canv->cd(2);
 h4->Add(h3,htemp__1,1,-1);
 h4->SetLineColor(kRed+1);
 h4->Scale(1/216.);
 h4->Draw("same");


  

// Double_t fitf(Double_t *x, Double_t *par)
//    {
     
//     Double_t fitval = par[0]*h1->GetBinContent() + par[1]*h2->GetBinContent();
      
//       return fitval;
// }

//  void myfit()
//  {
// TCanvas *c2 = new TCanvas("c2"," Fit",150,10,700,700);
//    //TGraphErrors* gS1 =new TGraphErrors(nlines,beta,s1,ebeta,es1);
//    //gStyle->SetOptFit(1111); gS1->SetTitle("Fit");
//    //gS1->SetName("Fit");
//    //gS1->SetMarkerStyle(20); gS1->SetMarkerSize(1.0); gS1->SetMarkerColor(1);
//       TH1F *h1 = (TH1F*)gDirectory->Get("h1");
//       TH1F *h2 = (TH1F*)gDirectory->Get("h2");
//       h1->Add(h2); 

//     TF1 *func = new TF1("fitf",fitf,-10,80000,2);
//     func->SetParNames("P0","P1");


//     h1->Fit("fitf","r");
//     h1->GetFunction("fitf")->SetLineColor(2);

// //Display data points and fit curves.
//    TCanvas *cSphere = new TCanvas("cSphere"," Fit",150,10,700,700);
//    h1->Draw("AP"); 
//    h1->GetXaxis()->SetTitle("Fit"); h1->GetYaxis()->SetTitle("fit");

//  }


 
}
