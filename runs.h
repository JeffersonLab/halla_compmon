#ifndef RUNS_H
#define RUNS_H

#include "TString.h"
#include "TRandom3.h"

#include <fstream>
#include <vector>
#include <string>

TString insults[15] = {"hippopotamus breath", "ya froghead", "weirdo", "oatmeal-for-brains", "cheesebrain",
                       "lunkhead", "poindexter", "ya weapons-grade buttface", "wiseguy", "commander witless", 
                       "ya clown", "sergeant major toe jam", "ya noodleloaf", "beetle-butt", "lieutenant 'I'm-so-good-at-following-directions'"};

const char* randInsult(){
  TRandom3 *r = new TRandom3(0);
  return insults[(int)(15*r->Rndm())].Data();
}

vector<vector<int>> productionRunList(int prexOrCrex){
  int snailMin = 1; int snailMax = 40;
  vector<vector<int>> runList;
  if(prexOrCrex == 2){
    snailMin = 101; snailMax = 220;
  }
  else if(prexOrCrex != 1 && prexOrCrex != 2){
    printf("Please enter correct parameter for experiment:\n");
    printf("    PREX == 1\n");
    printf("    CREX == 2\n\n");
    printf("Get it right this time, %s.\n", randInsult());
  }

  for(int snail = snailMin; snail <= snailMax; snail++){
    printf("Reading runs from file snail%i.list...\n", snail);
    string runNumStr;
    vector<int> snailList;
    ifstream infile(Form("%s/snail%i.list", getenv("COMPMON_SNAILS"), snail));
    if(infile.good()){
      while(getline(infile, runNumStr)){
        Int_t runNum = atoi(runNumStr.c_str());
        if(runNum == 0) continue;
        snailList.push_back(runNum);
      }
      runList.push_back(snailList);
      snailList.clear();
    }
  }
  return runList;
}



#endif
