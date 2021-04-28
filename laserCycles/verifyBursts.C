#include "laserUtils.h"
#include "../online/utils.h"
#include "../online/runs.h"

using namespace std;

vector<vector<vector<vector<Int_t>>>> sortBurstsByCycle(vector<vector<int>> allBursts){
  vector<vector<vector<vector<Int_t>>>> sortedBursts;
  if(allBursts.size() == 0) return sortedBursts;
  
  Int_t curCycle = 0; Int_t curPeriod = 0;
  vector<vector<vector<Int_t>>> cycleBursts;
  vector<vector<Int_t>> periodBursts;
  for(Int_t i = 0; i < allBursts.size(); i++){
    vector<Int_t> burstLimits;
    if(i == 0){
      curCycle = allBursts[i][0]; curPeriod = allBursts[i][1];
    }
    if(curPeriod != allBursts[i][1]){
      cycleBursts.push_back(periodBursts);
      periodBursts.clear(); curPeriod = allBursts[i][1];
    }
    if(curCycle != allBursts[i][0]){
      sortedBursts.push_back(cycleBursts);
      cycleBursts.clear(); curCycle = allBursts[i][0];
    }

    burstLimits.push_back(allBursts[i][2]); burstLimits.push_back(allBursts[i][3]);
    periodBursts.push_back(burstLimits);
  }
  cycleBursts.push_back(periodBursts);
  sortedBursts.push_back(cycleBursts);

  return sortedBursts;
}

void verifyBursts(Int_t prexOrCrex){
  vector<vector<int>> runList = productionRunList(prexOrCrex);

  for(vector<int> snailRunList : runList){
    for(Int_t run : snailRunList){
      vector<vector<int>> burstList = getCycleList(run, true);
      vector<vector<vector<vector<Int_t>>>> sortedBursts = sortBurstsByCycle(burstList);
      for(Int_t cycle = 0; cycle < sortedBursts.size(); cycle++){
        size_t per1 = sortedBursts[cycle][0].size();
        size_t per2 = sortedBursts[cycle][1].size();
        size_t per3 = sortedBursts[cycle][2].size();

        //printf("Run %i, Cycle %i, Sizes of period bursts: %zu / %zu / %zu\n", run, cycle+1, per1, per2, per3);

        if(per1 <= 0)
          printf("Run %i, Cycle %i, Period 1 has zero bursts!\n", run, cycle+1);
        if(per2 <= 0)
          printf("Run %i, Cycle %i, Period 2 has zero bursts!\n", run, cycle+1);
        if(per3 <= 0)
          printf("Run %i, Cycle %i, Period 3 has zero bursts!\n", run, cycle+1);
        if(per1 > 20)
          printf("Run %i, Cycle %i, Period 1 has an excess of bursts! (%zu)\n", run, cycle+1, per1);
        if(per2 > 20)
          printf("Run %i, Cycle %i, Period 2 has an excess of bursts! (%zu)\n", run, cycle+1, per2);
        if(per3 > 20)
          printf("Run %i, Cycle %i, Period 3 has an excess of bursts! (%zu)\n", run, cycle+1, per3);
      }
    }
  }
}
