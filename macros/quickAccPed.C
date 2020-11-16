void quickAccPed(Int_t run = 3803)
{
  TChain *mpswise = new TChain("mpswise");
  mpswise->Add(TString::Format("/data/cmuwork/rootfiles/tests_2018/compmon_%d.root",run));
  TCut mpsCut = "helicityState==0";
  TCanvas *canv = new TCanvas("canv","canv",500*2,400);
  canv->Divide(2,1);
  canv->cd(1);
  //TH1F *hRaw = new TH1F("hRaw","Acc0/NAcc0 no Ped Correction",100,0,0);
  mpswise->Draw("Acc0/NAcc0>>hRaw",mpsCut);
  TH1F *hRaw = (TH1F*)gDirectory->Get("hRaw");
  Double_t meanRaw = hRaw->GetMean();
  Double_t sigmaRaw = hRaw->GetRMS();
  Double_t multSigma = 1.0;
  TFitResultPtr sRaw = hRaw->Fit("gaus","S","",meanRaw-multSigma*sigmaRaw,meanRaw+multSigma*sigmaRaw);
  TString pedStr = TString::Format("%20.10f",-sRaw->Value(1));
  std::cout << "Ped: " << pedStr.Data() << std::endl;
  mpswise->SetAlias("Acc0PedSub",TString::Format("Acc0+NAcc0*%s",pedStr.Data()));
  canv->cd(2);
  mpswise->Draw("Acc0PedSub>>hPed",mpsCut);
  TH1F *hPed = (TH1F*)gDirectory->Get("hPed");
  Double_t meanPed = hPed->GetMean();
  Double_t sigmaPed = hPed->GetRMS();
  TFitResultPtr sPed = hPed->Fit("gaus","S","",meanPed-multSigma*sigmaPed,meanPed+multSigma*sigmaPed);

  TString pedSubStr = TString::Format("%20.10f",-sPed->Value(1));
  std::cout << "After Ped Sub: " << pedSubStr.Data() << std::endl;
}
