#ifndef vars_h
#define vars_h

#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"

#include <vector>
#include <string>
#include <fstream>
#include <istream>

typedef struct {Float_t mean, meanErr, rms, rmsErr;} DataVar;
typedef struct {Float_t mean, meanErr;} PolVar;
typedef struct {Float_t mean;} StdVar;
typedef struct {Float_t mean, meanErr, Chi2; Int_t NDF;} FitPolVar;

#endif 
