#include "TROOT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"


class Root {
  
 private:
  TH1D  *ch1, *ch2, *ch3, *ch4, *ch5, *ch6, *ch7, *ch8;
  TH1D *ch9, *ch10, *ch11, *ch12, *ch13, *ch14, *ch15, *ch16;

  TH1D  *rch1, *rch2, *rch3, *rch4, *rch5, *rch6, *rch7, *rch8;
  TH1D *rch9, *rch10, *rch11, *rch12, *rch13, *rch14, *rch15, *rch16;

  TFile *rootfile;

  
  FILE *eventfile;
  char *eventname;
  
  double adc_thres[16];
  double ped_mean[16];



 public:
  int openrootfile(char*,int);
  int createhistos(int);
  int createADChistos(int);
  int deletehistos(int);
  int deleteADChistos(int);
  int createdirectory(char*,int);
  int gotodirectory(char*,int);
  int fillhistos(double event[][512],int);
  int fillADChistos(double event[][512],int);
  int writehistos(int);
  int closerootfile(int);
  int get_adc_thres(char*, int, int);
  int eventempty(double event[][512], int,int); 
};
