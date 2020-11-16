const Int_t kNumTrigsPerMPS = 5;

void checkPrescale(Int_t run = 4299)
{
  TChain *chain = new TChain("mpswise");
  //chain->Add(TString::Format("/data/cmuwork/rootfiles/compmon_%d.root",run));
  chain->Add(TString::Format("$COMP_ROOTFILES/compmon_%d.root",run));

  chain->Draw(Form("scaler_run7/%d",kNumTrigsPerMPS),"laserState==1&&beamState==1");
} 
