void SetupColors(TTree *tree,Color_t color)
{
  tree->SetLineColor(color);
  tree->SetMarkerColor(color);
}

void inspectMiniRuns(Int_t run = 4299, Bool_t readFromFile = false)
{
  TChain *mpswise = new TChain("mpswise");
  mpswise->Add(TString::Format("$COMP_ROOTFILES/compmon_%d.root",run));

  if(readFromFile) {
    
  } else {
    TCanvas *canvasMpsTime = new TCanvas("canvasMpstime","canvasMpsTime",1600,600);
    //canvasMpsTime->Divide(1,2);
    SetupColors(mpswise,kBlack+1);
    mpswise->Draw("Acc0/NAcc0:mpsCoda","beamState==1&&Acc0/NAcc0>4");
    SetupColors(mpswise,kGreen+1);
    mpswise->Draw("Acc0/NAcc0:mpsCoda","laserState==1&&beamState==1","SAME");
    SetupColors(mpswise,kRed+1);
    mpswise->Draw("Acc0/NAcc0:mpsCoda","laserState==3&&beamState==1","SAME");
    gPad->SetGrid(1,1);
  }
} 
