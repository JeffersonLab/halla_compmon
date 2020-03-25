#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

vector<vector<int>> readEventCuts(int runNum){
  vector<vector<int>> eventCuts;
  ifstream mapfile(Form("%s/compmon_event_cuts_%i.map", getenv("COMPMON_MAPS"), runNum));
  if(!mapfile.good()){return eventCuts;}
  string readStr;
  while(getline(mapfile, readStr)){
    vector<int> cutRegion; stringstream ss(readStr);
    for(int i; ss >> i;){
      cutRegion.push_back(i);
      if(ss.peek() == ',')
        ss.ignore();
    }
  }
  return eventCuts;
}

bool isValidCycle(int *cycleInds, vector<int> section1, vector<int> section2, vector<int> section3, int runNum){
  int diff1 = section2[0] - section1[0];
  int diff2 = section3[0] - section2[0];
  bool temporallyLocked = diff1 == 1 && diff2 == 1;
  bool correctPattern = section1[1] == 0 && section2[1] == 1 && section3[1] == 0;
  vector<vector<int>> eventCuts = readEventCuts(runNum);
  bool noEventCuts = true;
  for(vector<int> cut : eventCuts){
    noEventCuts = noEventCuts && ((cut[0] <= section1[0] && cut[1] >= section1[0]) || 
                                  (cut[0] <= section3[1] && cut[1] >= section3[1]) || 
                                  (cut[0] <= section1[0] && cut[1] >= section3[1]));
  }
  if(!noEventCuts){
    printf("Applied event cuts on run %i\n", runNum);
  }
  //cout<<"Diff1: "<<diff1<<"; Diff2: "<<diff2<<"; section1: "<<section1[1]<<"; section2: "<<section2[1]<<
  //" section3: "<<section3[1]<<"; Final bool: "<<(temporallyLocked && correctPattern)<<endl;
  return temporallyLocked && correctPattern && noEventCuts;
}

bool indsAreValid(int *cycleInds){
  return cycleInds[0] >= 0 && cycleInds[1] > 0 && cycleInds[2] > 1;
}

//bool isValidSection(vector<int> section){
//  
//}

int numEntries(int *cycleInds){
  if(cycleInds[0] >= 0 && cycleInds[1] > 0 && cycleInds[2] > 1){return 3;}
  else if(cycleInds[0] >= 0 && cycleInds[1] > 0){return 2;}
  else if(cycleInds[0] >= 0){return 1;}
  else {return 0;}
}

vector<vector<int>> readCycles(int run_num){
  vector<vector<int>> cycles;

  string data_str;
  ifstream datafile(Form("%s/laserCycles_%i.dat", getenv("COMPMON_MINIRUNS"), run_num));
  while(getline(datafile, data_str)){
    vector<int> limits; stringstream ss(data_str);
    for(int i; ss >> i;){
      limits.push_back(i);
      if(ss.peek() == ' ')
        ss.ignore();
    }
    cycles.push_back(limits);
  }
  return cycles;
}

vector<vector<int>> refineCycles(vector<vector<int>> cycles){
  vector<vector<int>> shortCycles;

  for(int i = 0; i < cycles.size(); i++){
    printf("Cut vars for period %i - status: %i, beamOffFraction %.4f\n", cycles[i][0], cycles[i][9], cycles[i][6]*1.0/cycles[i][4]);
    if(cycles[i][9] == 1 && cycles[i][6]*1.0/cycles[i][4] < 0.2){
      shortCycles.push_back(cycles[i]);
    }
  }

  return shortCycles;
}

vector<vector<int>> defineCycles(vector<vector<int>> cycles, int runNum){
  vector<vector<int>> finalCycles;
  vector<int> currentCycle;

  int cycleInds[3] = {-1, -1, -1};
  for(int i = 0; i < cycles.size(); i++){
    int ents = numEntries(cycleInds);
    //cout<<"cycleInds: {"<<cycleInds[0]<<", "<<cycleInds[1]<<", "<<cycleInds[2]<<"}; Ents: "<<ents<<"; i: "<<i<<"; cycle id: "<<cycles[i][0]<<"; laser id: "<<cycles[i][1]<<endl;
    if(ents == 0){
      if(cycles[i][1] == 0){cycleInds[0] = i;}
    }
    else if(ents == 1){
      if(cycles[i][1] == 0){cycleInds[0] = i;}
      else if(cycles[i][0] - cycles[cycleInds[0]][0] != 1){cycleInds[0] = -1; cycleInds[1] = -1; cycleInds[2] = -1;}
      else{cycleInds[1] = i;}
    }
    else if(ents == 2){
      if(cycles[i][1] == 1){cycleInds[0] = -1; cycleInds[1] = -1; cycleInds[2] = -1;}
      else if(cycles[i][0] - cycles[cycleInds[1]][0] != 1){
        cycleInds[0] = i; cycleInds[1] = -1; cycleInds[2] = -1;
      }
      else{
        cycleInds[2] = i;
        if(indsAreValid(cycleInds) && isValidCycle(cycleInds, cycles[cycleInds[0]], cycles[cycleInds[1]], cycles[cycleInds[2]], runNum)){
          int startEvt = (cycles[cycleInds[0]][2] + cycles[cycleInds[0]][3])/2;
          int finalEvt = (cycles[cycleInds[2]][2] + cycles[cycleInds[2]][3])/2;
          currentCycle.push_back(startEvt); currentCycle.push_back(cycles[cycleInds[0]][3]);
          currentCycle.push_back(cycles[cycleInds[1]][2]); currentCycle.push_back(cycles[cycleInds[1]][3]);
          currentCycle.push_back(cycles[cycleInds[2]][2]); currentCycle.push_back(finalEvt);
          finalCycles.push_back(currentCycle);
          currentCycle.clear();
          cout<<"Cycle IDs: "<<cycles[cycleInds[0]][0]<<", "<<cycles[cycleInds[1]][0]<<", "<<cycles[cycleInds[2]][0]<<endl;
          cycleInds[0] = i; cycleInds[1] = -1; cycleInds[2] = -1;
        }
        else if(indsAreValid(cycleInds) && not isValidCycle(cycleInds, cycles[cycleInds[0]], cycles[cycleInds[1]], cycles[cycleInds[2]], runNum) && cycles[i][1] == 0){
          cycleInds[0] = i; cycleInds[1] = -1; cycleInds[2] = -1;
        }
        else{
          cycleInds[0] = -1; cycleInds[1] = -1; cycleInds[2] = -1;
        }
      }
    }
    else{
      cycleInds[0] = -1; cycleInds[1] = -1; cycleInds[2] = -1;
    }
  }

  return finalCycles;
}

void laserCycles(int snail_num){
  string run_num_str;
  ifstream infile(Form("%s/snail%i.list", getenv("COMPMON_SNAILS"), snail_num));
  while(getline(infile, run_num_str)){
    int run_num = atoi(run_num_str.c_str());
    vector<vector<int>> cycles = defineCycles(refineCycles(readCycles(run_num)), run_num);
    cout<<"Identified "<<cycles.size()<<" laser cycles."<<endl;

    ofstream output;
    output.open(Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), run_num));
    for(int cycle = 0; cycle < cycles.size(); cycle++){
      output<<cycles[cycle][0]<<","<<cycles[cycle][1]<<","<<cycles[cycle][2]<<","<<cycles[cycle][3]<<","<<cycles[cycle][4]<<","<<cycles[cycle][5]<<endl;
    }
    output.close();
    cout<<"Created file "<<Form("%s/cycles_%i.dat", getenv("COMPMON_MINIRUNS"), run_num)<<endl;
  }
  infile.close();
}
