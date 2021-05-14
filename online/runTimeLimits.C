#include "laserUtils.h"
#include "../online/utils.h"
#include "../online/runs.h"

using namespace std;


Int_t numDaysInMonth(Int_t month, Int_t year){
  if(month == 4 || month == 6 || month == 9 || month == 11){return 30;}
  else if(month == 2){
    if(year == 2012 || year == 2016 || year == 2020 || year == 2024){return 29;}
    else{return 28;}
  }
  else{return 31;}
}


vector<Int_t> dateNumbers(TTree *epicswise){
  vector<Int_t> dateData;
  string *date_ptr = 0;
  epicswise->SetBranchAddress("epics_datestring", &date_ptr);
  epicswise->GetEntry(0);

  TString date_str(*date_ptr); TObjArray *date_split = date_str.Tokenize(" ");
  TString month(      ( (TObjString *)(date_split->At(1)) )->String());
  TString date_num_str( ( (TObjString *)(date_split->At(2)) )->String());
  TString time_str    ( ( (TObjString *)(date_split->At(3)) )->String());
  TString year_str    ( ( (TObjString *)(date_split->At(5)) )->String());

  TObjArray *time_ptr = time_str.Tokenize(":");
  TString hour_str(   ( (TObjString *)(time_ptr->At(0)) )->String());
  TString minute_str( ( (TObjString *)(time_ptr->At(1)) )->String());
  TString second_str( ( (TObjString *)(time_ptr->At(2)) )->String());

  dateData.push_back(atoi(year_str.Data()));
  dateData.push_back(translateMonth(month)));
  dateData.push_back(atoi(date_num_str.Data()));
  dateData.push_back(atoi(hour_str.Data()));
  dateData.push_back(atoi(minute_str.Data()));
  dateData.push_back(atoi(second_str.Data()));

  return dateData;
}


void runTimeLimits(Int_t prexOrCrex){
  vector<vector<int>> runList = productionRunList(prexOrCrex);
  Int_t helFreq = 120;
  if(prexOrCrex == 1) helFreq = 240;
  
  for(Int_t i = 0; i < runList.size(); i++){
    for(Int_t j = 0; j < runList[i].size(); j++){
      TFile *f = TFile::Open(Form("%s/compmon_%i.root", getenv("COMP_ROOTFILES"), runList[i][j]));
      TTree *epicswise = f->Get("epicswise");
      TTree *mpswise = f->Get("mpswise");

      vector<Int_t> dateData = dateNumbers(epicswise);
      Int_t runTime = (Int_t)(mpswise->GetEntries()*1.0)/helFreq;
      Int_t endSeconds = dateData[5] + (runTime % 60);
      Int_t minuteAdd = 0;
      if(endSeconds > 59){endSeconds = endSeconds % 60; minuteAdd = 1;}
      Int_t endMinutes = dateData[4] + (((Int_t)(runTime/60)) % 60) + minuteAdd;
      Int_t hourAdd = 0;
      if(endMinutes > 59){endMinutes = endMinutes % 60; hourAdd = 1;}
      Int_t endHours = dateData[3] + ((Int_t)runTime/3600) + hourAdd;
      Int_t dateNumAdd = 0;
      if(endHours > 23){endHours = endHours % 24; dateNumAdd = 1;
      Int_t endDateNum = dateData[2] + dateNumAdd;
      Int_t monthAdd = 0;
      if(endDateNum > numDaysInMonth(dateData[1])){endDateNum = 1; monthAdd = 0;}
      Int_t endMonthNum = dateData[1] + monthAdd;
      Int_t endYear = dateData[0];

      ofstream outfile(Form("%s/Run%i_", getenv("COMPMON_WEB"), exptStr.Data()));
    }
  }
}
