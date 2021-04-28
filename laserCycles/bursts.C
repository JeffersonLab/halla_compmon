#include "laserUtils.h"
#include "../online/utils.h"

using namespace std;

void bursts(Int_t runNum, Int_t burstSize){
  vector<vector<int>> cycles = getCycleList(runNum);
  ofstream outfile; outfile.open(Form("%s/bursts_%i.dat", getenv("COMPMON_BURSTS"), runNum));

  vector<TChain *> runChains = loadChain(runNum);
  TTree *quartetwise = runChains[1];
  Int_t beamState, laserState, dithering, mpsNum;
  quartetwise->SetBranchAddress("laserState", &laserState);
  quartetwise->SetBranchAddress("beamState", &beamState);
  quartetwise->SetBranchAddress("dithering", &dithering);
  quartetwise->SetBranchAddress("firstMPSnumber", &mpsNum);

  for(Int_t i = 0; i < cycles.size(); i++){
    Int_t fOffStart = cycles[i][0]; Int_t fOffEnd = cycles[i][1];
    Int_t onStart = cycles[i][2]; Int_t onEnd = cycles[i][3];
    Int_t lOffStart = cycles[i][4]; Int_t lOffEnd = cycles[i][5];
    Int_t burstStart = -1; Int_t burstEnd = -1; Int_t burstCount = 0; Int_t burstState = -1;
    Int_t startEntry = 0;

    for(Int_t j = startEntry; j < quartetwise->GetEntries(); j++){
      quartetwise->GetEntry(j);
      if(mpsNum < fOffStart){continue;}
      if(mpsNum > lOffEnd) break;
      if(laserState < 4 && beamState == 1 && dithering == 0){
        if(mpsNum >= fOffStart && mpsNum <= fOffEnd){
          if(burstState != 0 && burstCount > 0){
            burstEnd = mpsNum - 4;
            if(burstCount >= 0.1*burstSize){
              outfile<<i+1<<","<<1<<","<<burstStart<<","<<burstEnd<<"\n";
            }
            burstStart = mpsNum; burstEnd = -1; burstCount = 1; burstState = 0;
          }
          else if(burstCount == 0){
            burstStart = mpsNum; burstCount++; burstState = 0;
          }
          else{
            burstCount++;
            if(burstCount >= burstSize){
              outfile<<i+1<<","<<1<<","<<burstStart<<","<<mpsNum<<"\n";
              burstCount = 0;
            }
          }
        }
        else if(mpsNum >= onStart && mpsNum <= onEnd){
          if(burstState != 1 && burstCount > 0){
            burstEnd = mpsNum - 4;
            if(burstCount >= 0.1*burstSize){
              outfile<<i+1<<","<<1<<","<<burstStart<<","<<burstEnd<<"\n";
            }
            burstStart = mpsNum; burstEnd = -1; burstCount = 1; burstState = 1;
          }
          else if(burstCount == 0){
            burstStart = mpsNum; burstCount++; burstState = 1;
          }
          else{
            burstCount++;
            if(burstCount >= burstSize){
              outfile<<i+1<<","<<2<<","<<burstStart<<","<<mpsNum<<"\n";
              burstCount = 0;
            }
          }
        }
        else if(mpsNum >= lOffStart && mpsNum <= lOffEnd){
          if(burstState != 0 && burstCount > 0){
            burstEnd = mpsNum - 4;
            if(burstCount >= 0.1*burstSize){
              outfile<<i+1<<","<<2<<","<<burstStart<<","<<burstEnd<<"\n";
            }
            burstStart = mpsNum; burstEnd = -1; burstCount = 1; burstState = 0;
          }
          else if(burstCount == 0){
            burstStart = mpsNum; burstCount++; burstState = 0;
          }
          else{
            burstCount++;
            if(burstCount >= burstSize){
              outfile<<i+1<<","<<3<<","<<burstStart<<","<<mpsNum<<"\n";
              burstCount = 0;
            }
          }
        }
      }
    }
    if(burstCount > 0){
      if(burstCount >= 0.1*burstSize){
        outfile<<i+1<<","<<3<<","<<burstStart<<","<<lOffEnd<<"\n";
      }
      burstCount = 0;
    }
    printf("  Finished with bursts for cycle %i!\n", i+1);
  }

  outfile.close();
}
