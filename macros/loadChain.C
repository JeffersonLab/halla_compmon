TChain *mpswise = 0;
TChain *quartetwise = 0;
TChain *pulserwise = 0;
TChain *triggerwise = 0;

void resetChain(TChain *ch)
{
  if(ch) {
   delete ch;
  }
}

void loadChain(Int_t runnum)
{


  resetChain(mpswise);
  resetChain(quartetwise);
  resetChain(pulserwise);
  resetChain(triggerwise);

  mpswise = new TChain("mpswise");
  quartetwise = new TChain("quartetwise");
  pulserwise = new TChain("pulserwise");
  triggerwise = new TChain("triggerwise");
  std::vector<TChain*> chains;
  chains.push_back(mpswise);
  chains.push_back(quartetwise);
  chains.push_back(pulserwise);
  chains.push_back(triggerwise);

  TString filesPre = Form("%s/compmon_%d",getenv("COMP_ROOTFILES"),runnum);
  int nfiles = 0;
  for(size_t ch = 0; ch < chains.size(); ch++) {
    nfiles = chains[ch]->Add(filesPre+".root");
    nfiles += chains[ch]->Add(filesPre+"_*.root");

    if(nfiles<=0) {
	    std::cerr << "Looked for files under: " << filesPre+".root" << std::endl;
	    std::cerr << "Found no files to plot!" << std::endl;
	    return;
    }
  }

}
