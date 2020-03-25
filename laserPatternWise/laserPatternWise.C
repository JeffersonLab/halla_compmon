#include "global.h"
#include <TMath.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <iostream>
#include <iomanip>
#include <fstream>

#define ACC 0
#define NUM_MEASUREMENTS 5 // H0, H1, Sum, Diff, Asym

//const bool kSaveIntHistos = false+1;
const bool kSaveIntHistos = false;
const bool kSaveIntHistos2D = false;
const bool kReadPedestals = false;
const double kNumSamplesAvg = 6601042;
const double kMinASub = -10;
const double kMaxASub =  10;
const bool kSetASubPlotLimits = true;
const TString kIntHistoPrefix="${COMPMON_WEB}/";
const Int_t kCanvSizeX = 1100; // US Letter size proportions (landscape)
const Int_t kCanvSizeY = 850; // US Letter size proportions (landscape)

const char *gMeasurementNames[5] = {"H0   ","H1   ","H1+H0","H1-H0","Asym "};
const char *gMeasureSafeNames[5] = {"H0","H1","Sum","Diff","Asym"};
const char *gUnitNames[5] = {"radc", "radc", "radc", "radc", "%"};
//const char *gStateNames[5] = {"OFF0","ON  ","OFF1","Bkg ","Sub "};
//const char *gStateSafeNames[5] = {"OFF0","ON","OFF1","Bkg","Sub"};

double gAnalyzingPowerOld;

std::vector<std::vector<std::vector<std::vector<double> > > > gResultsOld;
std::vector<ComptonVariableOld*> gVariablesOld;

static const char *kHelicityStructNames[5] = { "H0", "H1", "H0+H1", "H0-H1", "Asym" };
template<typename T>
struct HelicityStruct_t {
  T h0;
  T h1;
  T hSum;
  T hDiff;
  T hAsym;
  T null;
  T& operator[](int index) {
    switch(index) {
      case 0:
        return h0;
      case 1:
        return h1;
      case 2:
        return hSum;
      case 3:
        return hDiff;
      case 4:
        return hAsym;
      default:
        return null;
    }
  }

  T& sum() { return hSum; }
  T& diff() { return hDiff; }
};

struct VariableLimits_t {
  double low[NUM_MEASUREMENTS];
  double high[NUM_MEASUREMENTS];
  VariableLimits_t() {
    for(int m = 0; m < NUM_MEASUREMENTS; m++) {
      low[m] = 0;
      high[m] = 0;
    }
  };
};


std::map<TString, VariableLimits_t> gVariableLimits;

VariableLimits_t& getVariableLimits(std::string name)
{
  // Do we already have a map for this variable?
  std::pair< std::map<TString,VariableLimits_t>::iterator, bool> ret =
    gVariableLimits.insert(std::pair<TString,VariableLimits_t>(
          name,VariableLimits_t() ) );
  return ret.first->second;
}

void readVariableLimits()
{
  fstream in;
  in.open("compton_limits.dat",std::ios::in);
  int run;
  std::string svar,smeas,slow,shigh;

  std::string line;
  int ncomment = 0;
  int m = 0;
  double low,high;
  while(std::getline(in >> std::ws,line)) {
    m = -1;
    // Find any comments and get rid of them
    ncomment = line.find("#");
    if(ncomment != std::string::npos )
      line.erase(ncomment);
    ncomment = line.find(";");
    if(ncomment != std::string::npos )
      line.erase(ncomment);
    // If what follows has a run, acc and ath with its error, then continue
    std::istringstream iss(line);
    if( iss >> run && iss >> svar && iss >> smeas
        && iss >> slow && iss >> shigh) {
      // Is our run greater or equal to this run?
      if(run<=gRun) {
        TString tvar(svar);
        TString tmeas(smeas);
        TString tlow(slow);
        TString thigh(shigh);
        tvar.ToLower();
        tmeas.ToLower();
        VariableLimits_t &vp = getVariableLimits(tvar.Data());

        // Now that we have a pointer to the limits, let's find out
        // which limit we are setting
        if(tmeas.CompareTo("h")==0) {
          m = 0;
        } else if (tmeas.CompareTo("sum")==0) {
          m = 2;
        } else if (tmeas.CompareTo("diff") ==0) {
          m = 3;
        } else if (tmeas.CompareTo("asym") == 0) {
          m = 4;
        }

        if(m >= 0 && m <=4) {
          if (tlow.CompareTo("-") != 0) {
            low = tlow.Atof();
            vp.low[m] = low;
            if(m==0) {
              vp.low[1] = low; // Same as neg hel one
            }
          }
          if (thigh.CompareTo("-") != 0) {
            high = thigh.Atof();
            vp.high[m] = high;
            if(m==0) {
              vp.high[1] = high; // Same as neg hel one
            }
          }
        }
      }
    }
  }


  in.close();
}



struct Variables {
  std::vector<VComptonVariable*> vars;
  std::vector<LaserPatternObjectList<HelicityStruct_t<MyHistoData> > > results;
  std::vector<MyData> aths; // If this variable has a polarization
  std::vector<std::vector<MyData> > pols; // If this variable has a polarization
  std::vector<MyErrorWeightedData> avg_pol; // Average polarization for this variable (over all cycles)
  std::vector<LaserPatternObject<HelicityStruct_t<MyErrorWeightedData> > > avg_results;
  std::vector<VariableLimits_t> limits;
  void AddVariable(VComptonVariable* var, MyData *ath = 0) {
    vars.push_back(var);
    TString vname(var->GetName());
    vname.ToLower();
    limits.push_back(getVariableLimits(vname.Data()));
    aths.push_back(*ath);
    std::vector<MyData> tpol;
    tpol.resize(gLaserPatterns.size());
    pols.push_back(tpol);
    avg_pol.push_back(MyErrorWeightedData());
    results.push_back(LaserPatternObjectList<HelicityStruct_t<MyHistoData> >(&gLaserPatterns));
    avg_results.push_back(LaserPatternObject<HelicityStruct_t<MyErrorWeightedData> >());
    std::cout << "Added variable: " << var->GetName() << ", with limits: ";
    for(int m = 0; m < NUM_MEASUREMENTS; m++) {
      std::cout << gMeasureSafeNames[m] << "["
        << limits[limits.size()-1].low[m] << ", "
        << limits[limits.size()-1].high[m] << "]  ";
    }
    std::cout << std::endl;
  }
  size_t size() { return vars.size(); }

};

Variables gVariables;

std::vector<std::vector<std::vector<double> > > gAccResults;
std::vector<std::vector<std::vector<double> > > gBCMResults;
double gAccNoCycleONResults[10];
double gAccNoCycleOFFResults[10];
double gAccNoCycleSubResults[4];



template<typename T>
struct VarPatternList_t {
  std::vector<LaserPatternObjectList<HelicityStruct_t<T> > > list;
  VarPatternList_t() {
    Init();
  }
  void Init() {
    list.clear();
    for(size_t v = 0; v < gVariables.size(); v++) {
      list.push_back(LaserPatternObjectList<HelicityStruct_t<T> >(&gLaserPatterns));
    }
  };
  LaserPatternObjectList<HelicityStruct_t<T> >& operator[](int index) {
    return list[index];
  }
  size_t size() { return list.size(); }
};

void setupHisto(TH1 *h, int s) {
  h->SetLineWidth(2);
  h->SetLineColor(kGraphColors[s]);
  h->SetMarkerColor(kGraphColors[s]);
}

void backgroundSubtract( double *yOn, double *yOff, double *aOn, double *aOff,
    double *ped, double *results)
{
  double frac = yOff[0]/yOn[0];
  double fracE = TMath::Abs(frac)*TMath::Sqrt(
      TMath::Power(yOn[1]/yOn[0],2)+
      TMath::Power(yOff[1]/yOff[0],2)
      );
  double pedCorr = 0;
  if(kReadPedestals) {
    pedCorr = ped[0]/yOn[0];
  }
  double asym = (aOn[0]-frac*aOff[0])/(1-frac+pedCorr);
  double e0 = 1./(1.-frac);
  double asymE = TMath::Sqrt(
      TMath::Power(e0*aOn[1],2)+
      TMath::Power(e0*e0*(aOn[0]-aOff[0])*fracE,2)+
      TMath::Power(frac*e0*aOff[1],2)
      );
  results[0] = frac;
  results[1] = fracE;
  results[2] = asym;
  results[3] = asymE;
}

//void backgroundSubtract(std::vector<std::vector<double> > &results)
void backgroundSubtract(size_t v, size_t c)
{
  int on = 1;
  int bk = 3;
  int sb = 4;
  int yl = 4;
  int as = 8;
  double tmp[4];
  double yOn[2];
  double yBk[2];
  double aOn[2];
  double aBk[2];
  yOn[0] = gResultsOld[v][c][on][yl];
  yOn[1] = gResultsOld[v][c][on][yl+1];
  yBk[0] = gResultsOld[v][c][bk][yl];
  yBk[1] = gResultsOld[v][c][bk][yl+1];
  aOn[0] = gResultsOld[v][c][on][as];
  aOn[1] = gResultsOld[v][c][on][as+1];
  aBk[0] = gResultsOld[v][c][bk][as];
  aBk[1] = gResultsOld[v][c][bk][as+1];
  double ped[2] = {0.0};
  if(kReadPedestals&&v==1) {
    ped[0] = (gPedestals[1][2*c]-gPedestals[2][2*c])*kNumSamplesAvg;
    ped[1] = TMath::Sqrt(
        TMath::Power(gPedestals[1][2*c+1],2)+
        TMath::Power(gPedestals[2][2*c+1],2)
        );
  }

  backgroundSubtract(yOn,yBk,aOn,aBk,ped,tmp);
  gResultsOld[v][c][sb][0] = tmp[0];
  gResultsOld[v][c][sb][1] = tmp[1];
  gResultsOld[v][c][sb][2] = tmp[2];
  gResultsOld[v][c][sb][3] = tmp[3];
}

void printVariablePatternResults(size_t v)
{
  const char *prefix = "[%02d] %s %s %s ";
  for(int p = 0; p < int(gVariables.results[v].size()); p++) {
    printLine();
    for(size_t s = 0; s < LaserPatternObject<TH1F*>::size(); s++) {
      for(size_t m = 0; m < NUM_MEASUREMENTS; m++) {
        printResult(TString::Format("[%02d] %s %s %s ",p,
              gVariables.vars[v]->GetNameLeft(),gStateNames[s],
              gMeasurementNames[m]),
            gVariables.results[v][p][s][m],
            gVariables.vars[v]->GetUnits());
      }
      std::cout << std::endl;
    }
    // Print out an empty line
    if(gVariables.pols[v].size()>0) { // Has a polarization
      printResult(TString::Format("[%02d] %s Sub  Pol   ",p,
            gVariables.vars[v]->GetNameLeft()),
          gVariables.pols[v][p],"%");
      std::cout << std::endl;
    }

    /*
    printResult(TString::Format("[%02d] %s Bk / Sig   ",c,
          gVariables.vars[v]->GetNameLeft()), gVariables.results[v][c][4][0],
        gVariables.results[v][c][4][1],"");
    printResult(TString::Format("[%02d] %s Sub Asym   ",c,
          gVariables.vars[v]->GetNameLeft()), gVariables.results[v][c][4][2],
        gVariables.results[v][c][4][3],"%");
        */
  }

}

void printVariableCycleResults(size_t v)
{
  const char *prefix = "[%02d] %s %s %s ";
  for(int c = 0; c < int(gResultsOld[v].size()); c++) {
    printLine();
    for(size_t s = 0; s < NUM_STATES-1; s++) {
      for(size_t m = 0; m < NUM_MEASUREMENTS; m++) {
        printResult(TString::Format("[%02d] %s %s %s ",c,
              gVariablesOld[v]->GetNameLeft(),gStateNames[s],
              gMeasurementNames[m]),
            gResultsOld[v][c][s][2*m],gResultsOld[v][c][s][2*m+1],
            gVariablesOld[v]->GetUnits());
      }
      // Print out an empty line
      std::cout << std::endl;
    }
    if(gVariablesOld[v]->BackgroundSubtract()) {
      printResult(TString::Format("[%02d] %s Bk / Sig   ",c,
            gVariablesOld[v]->GetNameLeft()), gResultsOld[v][c][4][0],
          gResultsOld[v][c][4][1],"");
      printResult(TString::Format("[%02d] %s Sub Asym   ",c,
            gVariablesOld[v]->GetNameLeft()), gResultsOld[v][c][4][2],
          gResultsOld[v][c][4][3],"%");
    }
  }
}

void computeErrorWeighted(size_t v, size_t s, size_t m, double &mean, double &error,
    bool clear = true, bool finish = true)
{
  if(clear) {
    mean = 0;
    error = 0;
  }
  double tmp = 0;
  for(size_t c = 0; c < gResultsOld[v].size(); c++) {
    if(gLaserCyclesStatus[c]) {
      if(gResultsOld[v][c][s][m+1] != 0) {
        tmp    = TMath::Power(1./gResultsOld[v][c][s][m+1],2);
        mean  += gResultsOld[v][c][s][m]*tmp;
        error += tmp;
      }
    }
  }
  if(finish) {
    mean  = mean/error;
    error = TMath::Sqrt(1.0/error);
  }
}

void makeVariableGraph(size_t v)
{
  int printType = gVariablesOld[v]->GetPrintType();
  if(printType == 1) {
    TString titles[4];
    const char *name = gVariablesOld[v]->GetName();
    titles[0] = TString::Format("%s Laser-OFF Asymmetry;Entry;Asymmetry (%%)",name);
    titles[1] = TString::Format("%s Laser-ON  Asymmetry;Entry;Asymmetry (%%)",name);
    titles[2] = TString::Format("%s Bkg/Sig Fraction;Entry;Fraction",name);
    titles[3] = TString::Format("%s Bkg Sub Asymmetry;Entry;Asymmetry (%%)",name);
    std::vector<double > vecData[4][4];
    int states[4] = { 3, 1, 4, 4};
    int measurements[4] = { 8, 8, 0, 2};
    for(size_t c = 0; c < gResultsOld[v].size(); c++) {
      if(!gLaserCyclesStatus[c])
        continue;
      for(size_t k = 0; k < 4; k++) {
        vecData[k][0].push_back( 0.5*(gLaserCyclesStart[0][c] +
              gLaserCyclesEnd[2][c] ) );
        vecData[k][1].push_back( 0.0);
        vecData[k][2].push_back( gResultsOld[v][c][states[k]][measurements[k]] );
        vecData[k][3].push_back( gResultsOld[v][c][states[k]][measurements[k]+1] );
      }
    }
    TGraphErrors *graphs[4];
    TCanvas *canvas = new TCanvas("canvas","canvas",400*2,400*2);
    canvas->Divide(2,2);
    TH2F *hist = 0;
    for(size_t g =  0; g < 4; g++) {
      canvas->cd(g+1);
      graphs[g] = new TGraphErrors(vecData[g][0].size(),
          vecData[g][0].data(),vecData[g][2].data(),
          vecData[g][1].data(),vecData[g][3].data());
      graphs[g]->SetLineWidth(2);
      graphs[g]->SetLineColor(kGraphColors[g]);
      graphs[g]->SetMarkerColor(kGraphColors[g]);
      graphs[g]->SetMarkerStyle(kGraphStyles[g]);
      gPad->SetGrid(1,1);
      if(g<3||!kSetASubPlotLimits) {
        graphs[g]->Draw("AP 0");
      } else {
        hist = new TH2F(TString::Format("hTemp%s",name),
            titles[3],
            10,graphs[2]->GetXaxis()->GetXmin(),graphs[2]->GetXaxis()->GetXmax(),
            10,kMinASub,kMaxASub);
        hist->SetStats(0);
        hist->Draw();
        graphs[g]->Draw("SAME P0");
      }
      graphs[g]->SetTitle(titles[g]);
    }
    canvas->SaveAs(TString::Format("results/g%s_asymmetries_%d.png",name,gRun));
    delete canvas;
    if(hist)
      delete hist;
  }
}

void makeGraphs()
{
  for(size_t v = 0; v < gResultsOld.size(); v++) {
    makeVariableGraph(v);
  }
}

void printVariableSummary(size_t v, std::fstream *out = 0)
{
  TString prefix = TString::Format("Avg  %s ",gVariables.vars[v]->GetNameLeft());
  const char *units = gVariables.vars[v]->GetUnits();
  // First, identify what types should be printed
  int printType = gVariables.vars[v]->GetPrintType();
  if(printType == 1) {
    if(out) {
      printLineF(*out);
    } else {
      printLine();
    }
    for(size_t s = 0; s < LaserPatternObject<TH1F*>::size(); s++) {
      for(size_t m = 0; m < NUM_MEASUREMENTS; m++) {
        printResult(prefix + TString::Format("%s %s ",gStateNames[s],
              gMeasurementNames[m]),gVariables.avg_results[v][s][m],
            gVariables.vars[v]->GetUnits(),out);
      }
      if(out) {
        *out << std::endl;
      } else {
         std::cout << std::endl;
      }

    }
    if(out) {
      *out << std::endl;
    } else {
      std::cout << std::endl;
    }
    printResult(prefix+"Sub  Pol   ",
        gVariables.avg_pol[v],"%",out);
    // Print out an empty line
    if(out) {
      *out << std::endl;
    } else {
      std::cout << std::endl;
    }

  }/* else if (printType == 2) {
    double yield[2] = {0};
    double diff[2] = {0};
    double asym[2] = {0};
    computeErrorWeighted(v,1,4,yield[0],yield[1],false,false); // Sum L-ON
    computeErrorWeighted(v,3,4,yield[0],yield[1],false); // Sum Bk
    printResult(prefix+TString("H0+H1      "),yield[0],yield[1],units);

    computeErrorWeighted(v,1,6,diff[0],diff[1],false,false); // Diff L-ON
    computeErrorWeighted(v,3,6,diff[0],diff[1],false); // Diff Bk
    printResult(prefix+TString("H0-H1      "),diff[0],diff[1],units);

    computeErrorWeighted(v,1,8,asym[0],asym[1],false,false); // Asym L-ON
    computeErrorWeighted(v,3,8,asym[0],asym[1],false); // Asym Bk
    printResult(prefix+TString("Asym       "),asym[0],asym[1],"%");
  }
  */

}

void printVariableSummaryOld(size_t v)
{
  TString prefix = TString::Format("Avg  %s ",gVariablesOld[v]->GetNameLeft());
  const char *units = gVariablesOld[v]->GetUnits();
  // First, identify what types should be printed
  int printType = gVariablesOld[v]->GetPrintType();
  if(printType == 1) {
    // Accumulator printouts
    // Print out the average asymmetry
    double asymOn[2];
    double asymOff[2];
    double yieldOn[2];
    computeErrorWeighted(v,1,8*0+4,asymOn[0],asymOn[1]);
    computeErrorWeighted(v,3,8*0+4,asymOff[0],asymOff[1]);
    computeErrorWeighted(v,1,4,yieldOn[0],yieldOn[1],true,true);
    printResult(prefix+TString("L-ON  Asym "),asymOn[0],asymOn[1],"%");
    printResult(prefix+TString("L-OFF Asym "),asymOff[0],asymOff[1],"%");
    printResult(prefix+TString("L-ON  Sum  "),yieldOn[0],yieldOn[1],"radc");

    if(gVariablesOld[v]->BackgroundSubtract()) {
      double frac[2];
      double asymSub[2];
      computeErrorWeighted(v,4,0,frac[0],frac[1]);
      computeErrorWeighted(v,4,2,asymSub[0],asymSub[1]);
      printResult(prefix+TString("Bk / Sig   "),frac[0],frac[1],"%");
      printResult(prefix+TString("BkSub Asym "),asymSub[0],asymSub[1],"%");

      if(gVariablesOld[v]->ComputePol()) {
        double pol[2];
        pol[0] = 100*asymSub[0]/gAnalyzingPowerOld;
        pol[1] = 100*asymSub[1]/gAnalyzingPowerOld;
        printResult(prefix+TString("BkSub Pol  "),pol[0],pol[1],"%");
      }
    }
  } else if (printType == 2) {
    double yield[2] = {0};
    double diff[2] = {0};
    double asym[2] = {0};
    computeErrorWeighted(v,1,4,yield[0],yield[1],false,false); // Sum L-ON
    computeErrorWeighted(v,3,4,yield[0],yield[1],false); // Sum Bk
    printResult(prefix+TString("H0+H1      "),yield[0],yield[1],units);

    computeErrorWeighted(v,1,6,diff[0],diff[1],false,false); // Diff L-ON
    computeErrorWeighted(v,3,6,diff[0],diff[1],false); // Diff Bk
    printResult(prefix+TString("H0-H1      "),diff[0],diff[1],units);

    computeErrorWeighted(v,1,8,asym[0],asym[1],false,false); // Asym L-ON
    computeErrorWeighted(v,3,8,asym[0],asym[1],false); // Asym Bk
    printResult(prefix+TString("Asym       "),asym[0],asym[1],"%");
  }
}

void printSummary(std::fstream *out = 0)
{
  if(out) {
    printLineF(*out);
  } else {
    printLine();
  }
  if(out) {
    (*out) << "Run " << gRun << " analysis summary:" << std::endl;
  } else {
    std::cout << "Run " << gRun << " analysis summary:" << std::endl;
  }
  for(size_t v = 0; v < gVariables.vars.size(); v++) {
    if(v>0) {
      if(out) {
        *out << std::endl;
      } else {
        std::cout << std::endl;
      }
     }
    printVariableSummary(v,out);
  }
}

void printSummaryOld()
{
  printLine();
  std::cout << "Run " << gRun << " analysis summary:" << std::endl;
  for(size_t v = 0; v < gVariablesOld.size(); v++) {
    if(v>0)
      std::cout << std::endl;
    printVariableSummary(v);
  }
  printLine();
}

struct MyHistoVarInfo_t {
  TString varname;
  int pat;
  TString measname;
  float x0;
  float y0;
  float x;
  float y;
  float w;
  float h;


  TPaveLabel *pVar;
  TPaveLabel *pPat;
  TPaveLabel *pMea;
  MyHistoVarInfo_t(TString vname = "", int p = 0,TString mname = "") :
    varname(vname),pat(p),measname(mname), pVar(0), pPat(0), pMea(0) {
      Init();
  }
  void Init() {
    x0 = gStyle->GetPadLeftMargin();
    y0 = 0.95;
    w = (1.0-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())/3.;
    h = 0.05;
    x = y = 0;
    pVar = new TPaveLabel();
    pPat = new TPaveLabel();
    pMea = new TPaveLabel();
    Customize(pVar);
    Customize(pPat);
    Customize(pMea);
  }
  void Customize(TPaveLabel *p) {
    p->SetBorderSize(1);
    p->SetFillColor(kWhite);
    p->SetFillStyle(1001);
  }
  void SetPat(int p) {
    pat = p;
  }
  void SetMeasuremennt(TString mname) {
    measname = mname;
  }
  void SetVariable(TString vname) {
    varname = vname;
  }
  void Draw() {
    if(!pVar||!pPat||!pMea)
      return;
    x = x0;
    y = y0;
    pVar->DrawPaveLabel(x,y,x+w,y+h,"Variable: "+varname,"NDC NB");
    x+=w;
    pMea->DrawPaveLabel(x,y,x+w,y+h,"Type: "+measname,"NDC NB");
    x+=w;
    pPat->DrawPaveLabel(x,y,x+w,y+h,TString::Format("Pattern: %02d",pat),"NDC NB");
  }
};

struct MyHistoTitleBox_t {
  TPaveText *pbox;
  TPaveLabel *ptext;
  TString pname;
  bool init;
  MyHistoTitleBox_t(TString name = "") : pbox(0), ptext(0), pname(name),
    init(0) {
  }
  void Init(int color, float x, float y, float w, float h,
      float xlw,float ylw) {
    ptext = new TPaveLabel(x+xlw,y+ylw,x+w-xlw,y+h+ylw,pname,"NDC NB");
    ptext->SetBorderSize(1);
    ptext->SetFillColor(kWhite);
    pbox = new TPaveText(x,y,x+w,y+h,"NDC");
    pbox->SetFillColor(color);
    pbox->SetFillStyle(1001);
    pbox->SetBorderSize(0);
    init = true;
  }
  void Draw()
  {
    if(!ptext||!pbox)
      return;

    pbox->Draw();
    ptext->Draw();
  }
};

struct HistoDrawHelper_t {
  float pvw; // PaveStat width
  float pvh; // PaveStat height
  float pvx; // Current PaveStat x position
  float pvy; // Current PaveStat y position
  bool init;
  bool reuse;
  int dim;
  MyHistoVarInfo_t *measInfo;

  std::vector<TH1*> histos; // List of histograms for easy drawing
  std::vector<MyHistoTitleBox_t*> btitles; // Helper titles to use for these histos

  // Simple constructor
  HistoDrawHelper_t(int d = 1, MyHistoVarInfo_t *info = 0) :
    dim(d), measInfo(info) {
    Clear();
    reuse = false;
  }
  void Add(TH1 *h, TString name = "") {
    histos.push_back(h);
    if(!reuse) {
      btitles.push_back(new MyHistoTitleBox_t(name));
    }
  }
  void NewRepeatList() {
    histos.clear(); // Clear histos only, but preserve box titles
    ResetPavePositions();
    reuse = true;
  }
  void NewList() {
    NewRepeatList();
    histos.clear();
    btitles.clear();
  }
  void Clear()
  {
    pvw = pvh = pvx = pvy = 0.0;
    ResetPavePositions();
    init = reuse = false;
    NewList();
  }
  void DrawHistos()
  {
    // Loop over all histograms and draw them
    if (histos.size() == 0)
      return; // Nothing to do

    // If it's not initialized, do so now
    if(!init)
      Init();

    // First, make sure that the histogram is properly padded so that
    // the stat boxes won't overlap anything.
    // First, find the largest value of all histograms
    float maxy = -1e6;
    float tmp = 0;
    for(size_t i = 0; i < histos.size(); i++ ) {
      tmp = histos[i]->GetBinContent(histos[i]->GetMaximumBin());
      if(tmp>maxy)
        maxy=tmp;
    }

    // Now, set the max for the first histogram
    histos[0]->SetMaximum(maxy*1.1);

    std::vector<TPaveStats*> pvec;
    TH1 *hh = 0;
    TPaveStats *pv = 0;
    for(size_t i = 0; i < histos.size(); i++ ) {
      hh = histos[i];
      hh->SetTitle("");
      if(i==0) {
        hh->Draw("");
      } else {
        hh->Draw("SAMES");
      }
      gPad->Update(); // To ensure PaveStat has been placed

      // Now configure the PaveStats
      if(!reuse) { // The title boxes should be configured
        btitles[i]->Init(hh->GetLineColor(),pvx,pvy,pvw,0.05,0.00,0.01);
      }
      pv = (TPaveStats*)hh->FindObject("stats");
      if(pv) { // Found the pave stats, good!
        pv->SetX1NDC(pvx);
        pv->SetY1NDC(pvy);
        // Increment x position
        pvx += pvw;
        pv->SetX2NDC(pvx);
        pv->SetY2NDC(pvy-pvh);
        //pv->SetLineColor(hh->GetLineColor());
        //pv->SetFillColor(kWhite);
        pv->SetFillStyle(1001); // Solid style
        pvec.push_back(pv);
      } else { // Increment x position only
        pvx += pvw;
      }
    }
    // Draw the title boxes (even if no pv)
    for(size_t i = 0; i < btitles.size(); i++) {
      btitles[i]->Draw();
    }

    // Now, loop over the PaveStats to re-draw them over the histograms
    for(std::vector<TPaveStats*>::iterator it = pvec.begin();
        it != pvec.end(); it++) {
      (*it)->Draw();
    }

    // Draw the measurement info
    if(measInfo)
      measInfo->Draw();
    gPad->Update(); // Force it to re-draw all pave stats
  }
  void ResetPavePositions()
  {
    // Get the top right corner of the histogram in NDC units
    pvx = gStyle->GetPadLeftMargin();
    // Additional 0.05 lower to leave room for the stat box title
    //pvy = 1.0-gStyle->GetPadTopMargin()-0.06;
    pvy = 0.885;
    if(dim==2) {
      //pvy -= 0.0075;
    }
  }
  void Init() {
    init = true;
    ResetPavePositions();
    pvw = (1.0-gStyle->GetPadLeftMargin()-gStyle->GetPadRightMargin())/
      float(histos.size());
    pvh = 0.15;
    if(dim==2) {
      pvh = 0.20;
    }
  }
};

void walkTreeAndAnalyze()
{
  VarPatternList_t<TH1F*> histos;
  VarPatternList_t<TH2F*> histos2;
  std::vector<LaserPatternObjectList<HelicityStruct_t<MyHistoData> > >
    &results = gVariables.results;

  for(size_t v = 0; v < gVariables.size(); v++) {
    for(size_t p = 0; p < histos[v].size(); p++) {
      for(int s = 0; s < LaserPatternObject<TH1F*>::size()-1; s++ ) {
        for(int m = 0; m < NUM_MEASUREMENTS; m++) {
          const char *vname = gVariables.vars[v]->GetName();
          const char *sname = kLaserPatternObjectNames[s];
          const char *mname = kHelicityStructNames[m];
          histos[v][p][s][m] = new TH1F(TString::Format("hV%s_P%03d_S%s_M%s",
                vname,int(p),sname,mname),
              TString::Format("%s %s %s P[%02d];%s (%s)",vname,
                gMeasurementNames[m],gStateNames[s],int(p),vname,
                gVariables.vars[v]->GetUnits()),
              100,gVariables.limits[v].low[m],gVariables.limits[v].high[m]);
          int mpslow = gLaserPatterns[p].offLeft->start;
          int mpshigh = gLaserPatterns[p].offRight->end;
          int mpsrange = (mpshigh-mpslow)/10;
          if (mpsrange>1e2)
            mpsrange=1e2;
          mpslow -= mpsrange*0.05;
          mpshigh += mpsrange*0.05;
          histos2[v][p][s][m] = new TH2F(TString::Format("h2V%s_P%03d_S%s_M%s",
                vname,int(p),sname,mname),
              TString::Format("%s %s %s P[%02d];mps;%s (%s)",vname,
                gMeasurementNames[m],gStateNames[s],int(p),vname,
                gVariables.vars[v]->GetUnits()),
              mpsrange,mpslow,mpshigh,
              100,gVariables.limits[v].low[m],gVariables.limits[v].high[m]);
          setupHisto(histos[v][p][s][m],s);
          setupHisto(histos2[v][p][s][m],s);
        }
      }
    }
  }

  bool goodPattern;
  std::vector<int> patternNumbers;
  std::vector<int> patternStates;
  std::vector<int> cycleNumbers;
  std::vector<int> cycleStates;
  std::vector<std::vector<double> > variablePatValues;

  double tmpVal[5];
  // Now step through the tree and extract the values
  for(int entry = gStartEntry; entry < gEntries; entry++) {
    // Get the entry from the multiplet tree
    gChain->GetEntry(entry);


    // Reset variables
    goodPattern = true; // The default is to mark it good

    // First, check that this entry falls inside a valid laser pattern
    patternNumbers.clear();
    patternStates.clear();
    for(size_t p = 0; p < gLaserPatterns.size(); p++) {
      int s = gLaserPatterns[p].mpsIsInState(mpsCoda,1);
      if(s != -1) {
        patternNumbers.push_back(p);
        patternStates.push_back(s);
      }
    }

    if(patternNumbers.size() == 0) {
      goodPattern = false;
    }

    // Check to make sure timing information is correct
    //if(TMath::Abs(double(runScaler2)/40e6-kAvgMPSTime)>kAvgMPSTimeDeviation)
    //  goodPattern = false;
    float bcm = (posBCM + negBCM)/2.0;
    if(beamState==1) {
      goodPattern = false;
    }

    // Check that the helicity is valid
    if(helicityState != 0 && helicityState != 1) {
      goodPattern = false;
    }

    // If all is well, then we can use this pattern
    if(goodPattern) {
      for(size_t j = 0; j < patternNumbers.size(); j++) {
        int p = patternNumbers[j];
        int s = patternStates[j];
        for(size_t v = 0; v < gVariables.size(); v++) {
          for(int h = 0; h < 2; h++) {
            tmpVal[h] = gVariables.vars[v]->val(h);
          }
          tmpVal[2] = tmpVal[1]+tmpVal[0];
          tmpVal[3] = tmpVal[1]-tmpVal[0];
          tmpVal[4] = tmpVal[3]/tmpVal[2];
          for(int m = 0; m < 5; m++) {
            histos[v][p][s][m]->Fill(tmpVal[m]);
            histos2[v][p][s][m]->Fill(mpsCoda,tmpVal[m]);
            if(s==0||s==2) { // bkg measurement
              histos[v][p][3][m]->Fill(tmpVal[m]);
              //histos2[v][p][3][m]->Fill(mpsCoda,tmpVal[m]);
            }
          }
        }
      }
    }

  }

  // At the very end, get the results from the histogram
  for(size_t v = 0; v < gVariables.size(); v++) {
    for(size_t p = 0; p < histos[v].size(); p++) {
      for(int s = 0; s < LaserPatternObject<TH1F*>::size()-1; s++ ) {
        for(int m = 0; m < 5; m++) {
          results[v][p][s][m].Set(histos[v][p][s][m]);
        }
      }

      // Now, compute the bkg subtracted values
      for(int m = 0; m < 4; m++) {
        MyData dat = results[v][p].getOn()[m] - results[v][p].getBkg()[m];
        results[v][p][4][m].Set(dat);
      }
      // Compute the background subtracted asymmetry
      results[v][p][4][4].Set(results[v][p][4][3]/results[v][p][4][2]);
      // If this variable has an analyzing power, then compute polarization
      // for this laser pattern
      if(gVariables.aths[v].val() != -1.0 ) {
        gVariables.pols[v][p].Set(100.*results[v][p][4][4]/gVariables.aths[v]);
      }
    }
    // Print out the pattern results for this variable
    //printVariablePatternResults(v);
    // Finally close any open PDF files we have
  }

  // Now compute the average weighted results
  for(size_t v = 0; v < gVariables.size(); v++) {
    for(size_t p = 0; p < histos[v].size(); p++) {
      for(int s = 0; s < LaserPatternObject<TH1F*>::size(); s++ ) {
        for(int m = 0; m < 5; m++) {
          gVariables.avg_results[v][s][m].add(results[v][p][s][m]);
        }
      }
      if(gVariables.aths[v].val() != -1.0 ) {
        gVariables.avg_pol[v].add(gVariables.pols[v][p]);
      }
    }
    for(int s = 0; s < LaserPatternObject<TH1F*>::size(); s++ ) {
      for(int m = 0; m < 5; m++) {
        gVariables.avg_results[v][s][m].finish();
      }
    }
    if(gVariables.aths[v].val() != -1.0 ) {
      gVariables.avg_pol[v].finish();
    }
  }

  // Save histograms?
  if(kSaveIntHistos) {
    TCanvas *cTemp = new TCanvas("cTemp","cTemp",kCanvSizeX,kCanvSizeY);
    cTemp->Divide(2,3,10./kCanvSizeX,10./kCanvSizeY);
    TCanvas *cTemp2 = new TCanvas("cTemp2","cTemp2",kCanvSizeX,kCanvSizeY);
    cTemp2->Divide(2,3,10./kCanvSizeX,10./kCanvSizeY);

    HistoDrawHelper_t *drawHelper1[NUM_MEASUREMENTS];
    HistoDrawHelper_t *drawHelper2[NUM_MEASUREMENTS];
    MyHistoVarInfo_t *measInfo[NUM_MEASUREMENTS];
    for(int m = 0; m < NUM_MEASUREMENTS; m++) {
      measInfo[m] = new MyHistoVarInfo_t("",0,gMeasureSafeNames[m]);
      drawHelper1[m] = new HistoDrawHelper_t(1,measInfo[m]);
      drawHelper2[m] = new HistoDrawHelper_t(2,measInfo[m]);
    }

    for(size_t v = 0; v < gVariables.size(); v++) {
      TString cTempfname = gComptonOutPath+TString::Format("/laserwise_1d_%s.pdf",
          gVariables.vars[v]->GetName());
      TString cTemp2fname = gComptonOutPath+TString::Format("/laserwise_2d_%s.pdf",
          gVariables.vars[v]->GetName());

      // Open the PDF for writing (nothing gets saved yet)
      cTemp->Print(cTempfname+"[");
      if(kSaveIntHistos2D) {
        cTemp2->Print(cTemp2fname+"[");
      }
      for(int m = 0; m < NUM_MEASUREMENTS; m++) {
        measInfo[m]->SetVariable(gVariables.vars[v]->GetName());
      }

      // Now loop through the patterns and fill the PDF pages with
      // each laser pattern
      for(size_t p = 0; p < histos[v].size(); p++) {
        // This time loop over measurements first
        for(int m = 0; m <NUM_MEASUREMENTS; m++) {
          measInfo[m]->SetPat(p);
          cTemp->cd(1+m);
          for(int s = 0; s < 3; s++) {
            drawHelper1[m]->Add(histos[v][p][s][m],
                gStateSafeNamesPretty[s]);
          }
          drawHelper1[m]->DrawHistos();
          drawHelper1[m]->NewRepeatList(); // Start a new list, but preserve old settings

          // Now do the same for the 2D histos
          cTemp2->cd(1+m);
          for(int s = 0; s < 3; s++) {
            drawHelper2[m]->Add(histos2[v][p][s][m],
                gStateSafeNamesPretty[s]);
          }
          drawHelper2[m]->DrawHistos();
          drawHelper2[m]->NewRepeatList(); // Start a new list, but preserve old settings

          /*
            if(s!=3) {
              histos[v][p][*si][m]->Draw(firstState?"":"SAME");
            if(s < LaserPatternObject<TH1F*>::size()-2) {
              cTemp2->cd(1+m);
              histos2[v][p][s][m]->Draw(firstState?"":"SAME");
            }
          }
          */
        }
        // Now append this page to the PDF
        cTemp->Print(cTempfname);
        if(kSaveIntHistos2D) {
          cTemp2->Print(cTemp2fname);
        }
      }
      // Now close the PDFs for this variable
      cTemp->Print(cTempfname+"]");
      if(kSaveIntHistos2D) {
        cTemp2->Print(cTemp2fname+"]");
      }
    }
  }


  // Cleanup
  for(size_t v = 0; v < gVariables.size(); v++) {
    for(size_t p = 0; p < histos[v].size(); p++) {
      for(int s = 0; s < LaserPatternObject<TH1F*>::size()-1; s++ ) {
        for(int m = 0; m < NUM_MEASUREMENTS; m++) {
          delete histos[v][p][s][m];
        }
      }
    }
  }
}

void walkTreeAndAnalyzeOld()
{
  std::vector<std::vector<std::vector<std::vector<TH1F*> > > > histograms;
  gResultsOld.resize(gVariablesOld.size());
  histograms.resize(gVariablesOld.size());
  std::vector<std::vector<double> > variablePatValues;
  std::vector<std::vector<double> > variableMeasurements;
  variablePatValues.resize(gVariablesOld.size());
  variableMeasurements.resize(gVariablesOld.size());
  for(size_t v = 0; v < gVariablesOld.size(); v++) {
    gResultsOld[v].resize(gLaserCyclesStart[0].size());
    histograms[v].resize(gLaserCyclesStart[0].size());
    variablePatValues[v].resize(kHelStructure);
    variableMeasurements[v].resize(2*NUM_MEASUREMENTS);

    for(size_t c = 0; c < gLaserCyclesStart[0].size(); c++) {
      gResultsOld[v][c].resize(NUM_STATES);
      histograms[v][c].resize(NUM_STATES-1); // No bkg sub histograms!
      for(size_t s = 0; s < NUM_STATES; s++) {
        gResultsOld[v][c][s].resize(2*NUM_MEASUREMENTS);
        if(s<NUM_STATES-1) {
          for(size_t m = 0; m < NUM_MEASUREMENTS; m++) {
            histograms[v][c][s].push_back(new TH1F(
                  TString::Format("h%sL%zu%s%s",
                  gVariablesOld[v]->GetName(),c,gStateSafeNames[s],
                  gMeasureSafeNames[m]),"",100,-5,5));
          }
        }
      }
    }
  }

  std::vector<double> rAccAsyms[6][4];

  int patternTypeActual;
  int patternTypeReported;
  int *patternType = &patternTypeActual;

  double result[5];
  double bcmResult[5];
  double adc[kHelStructure] = {0.};
  double bcmPat[kHelStructure] = {0.};
  bool goodPattern;
  int patternActual[kHelStructure] = {-1};
  int patternReported[kHelStructure] = {-1};
  std::vector<int> cycleNumbers;
  std::vector<int> cycleStates;

  // Now we are ready to walk the tree and make measurements
  for(int entry = gStartEntry; entry < gEntries-kHelStructure; entry+=kHelStructure) {

    // First, check that this entry falls inside a valid laser pattern;
    cycleNumbers.clear();
    cycleStates.clear();
    for(size_t i = 0; i < gLaserCyclesStart[0].size(); i++) {
      for( int j = 0; j < 3; j++) {
        if(mpsCoda >= gLaserCyclesStart[j][i] &&
            mpsCoda+kHelStructure <= gLaserCyclesEnd[j][i]) {
          cycleNumbers.push_back(i);
          cycleStates.push_back(j);
        }
      }
    }

    // If the cycle is not found, then don't bother with this helicity pattern
    //if(cycleNumbers.size() ==0)
    //  continue;

    // Reset variables
    goodPattern = true;

    // Loop through the four entries in the pattern
    for(int hel = 0; hel < kHelStructure; hel++) {
      // Get the next entry from the tree
      gChain->GetEntry(entry+hel);
      patternActual[hel] = helicityState;
      patternReported[hel] = helicityState;
      float bcm = (posBCM + negBCM)/2.0;      

      // Fill variable values
      for(size_t v = 0; v < gVariablesOld.size(); v++) {
        variablePatValues[v][hel] = gVariablesOld[v]->GetValue();
        //debugging
        if(v==1)
          variablePatValues[v][hel] /= NAcc0;

        if(kBCMNormalizeMPSLevel) {
          if(gVariablesOld[v]->BCMNormalize()) {
            if(bcm != 0) {
              variablePatValues[v][hel] /= bcm;
            } else {
              goodPattern = false;
            }
          }
        }
      }

      // Check to make sure timing information is correct
      //if(TMath::Abs(double(runScaler2)/40e6-kAvgMPSTime)>kAvgMPSTimeDeviation)
      //  goodPattern = false;

      if((kBeamOn&&bcm < kMinCurrent)||(!kBeamOn&&bcm>kMaxCurrent)) {
        goodPattern = false;
      }

      // Check that the helicity is valid
      if(helicityState != 0 && helicityState != 1) {
        goodPattern = false;
      }
    }

    patternTypeActual = GetPatternType(patternActual,kHelStructure);
    patternTypeReported = GetPatternType(patternReported,kHelStructure);

    // First, that both the actual and reported patterns are valid
    if(patternTypeActual == PATTERN_UNKNOWN
        || patternTypeReported == PATTERN_UNKNOWN ) {
      goodPattern = false;
      std::cout << "ERROR: This is an unknown pattern! "
       << PatternToString(patternActual,kHelStructure).c_str() << " vs "
       << PatternToString(patternReported,kHelStructure).c_str() << std::endl;
    }

    // If all is well, then we can use this pattern
    if(goodPattern) {
      for(size_t v = 0; v < gVariablesOld.size(); v++) {
        processPattern(*patternType,variablePatValues[v],
            variableMeasurements[v],kHelStructure);
        for(int m = 0; m < NUM_MEASUREMENTS; m++ ) {
          for(size_t c = 0; c < cycleNumbers.size(); c++) {
            histograms[v][cycleNumbers[c]][cycleStates[c]][m]->Fill(
                variableMeasurements[v][m]);
            // Also create one of the average background
            if(cycleStates[c] == 0 || cycleStates[c] == 2) {
              histograms[v][cycleNumbers[c]][3][m]->Fill(
                variableMeasurements[v][m]);
            }
          }
        }
      }
    }
  }

  // Once we finished walking the tree, now extract the results
  // from the histograms
  TH1F *hTemp = 0;
  TCanvas *cTemp = new TCanvas("cTemp","cTemp",400,400);
  for(size_t v = 0; v < gVariablesOld.size(); v++ ) {
    for(size_t c = 0; c < gLaserCyclesStart[0].size(); c++ ) { // Cycles
      for(size_t s = 0; s < NUM_STATES-1; s++ ) { // States
        for(size_t m = 0; m < NUM_MEASUREMENTS; m++) {
          hTemp = histograms[v][c][s][m];
          gResultsOld[v][c][s][2*m] = hTemp->GetMean();
          gResultsOld[v][c][s][2*m+1] = hTemp->GetRMS()/
            TMath::Sqrt(hTemp->GetEntries());
          if(kSaveIntHistos) {
            hTemp->Draw();
            cTemp->SaveAs(TString::Format("results/h%s_%d_c%02d_%d_m%d.png",
                  gVariablesOld[v]->GetName(),gRun,int(c),int(s),int(m)));
          }
          /*
             if(m==4) {
             if(s<3) {
             rAccAsyms[s][0].push_back(
             (gLaserCyclesStart[j][i]+gLaserCyclesEnd[j][i])/2);
             } else {
             rAccAsyms[j][0].push_back(
             (gLaserCyclesStart[0][i]+gLaserCyclesEnd[2][i])/2);
             }
             rAccAsyms[j][1].push_back(0.0);
             rAccAsyms[j][2].push_back(hAcc[j][m][i]->GetMean());
             rAccAsyms[j][3].push_back(hAcc[j][m][i]->GetRMS()/
             TMath::Sqrt(hAcc[j][m][i]->GetEntries()));
             }*/
        }
      }
      if(gVariablesOld[v]->BackgroundSubtract())
        backgroundSubtract(v,c);
      /*
         if(gAccResults[i][4][3]<15 && gAccResults[i][4][2] == gAccResults[i][4][2]) {
         for(int l = 0; l < 2; l++) {
         rAccAsyms[4+l][0].push_back((gLaserCyclesStart[0][i]+gLaserCyclesEnd[2][i])/2);
         rAccAsyms[4+l][1].push_back(0.0);
         rAccAsyms[4+l][2].push_back(gAccResults[i][4][2*l]);
         rAccAsyms[4+l][3].push_back(gAccResults[i][4][2*l+1]);
         }
         }
         */
    }
    printVariableCycleResults(v);
    //printCycleResults(gAccResults,TString::Format("Acc%d",ACC));
  }

  if(kSaveGraphs)
    makeGraphs();
    /*
    TGraphErrors *graphs[6];
    TF1 *fits[6];
    for(int i = 0; i < 6; i++ ) {
      graphs[i] = new TGraphErrors(rAccAsyms[i][0].size(),
          rAccAsyms[i][0].data(),rAccAsyms[i][2].data(),
          rAccAsyms[i][1].data(),rAccAsyms[i][3].data());
      graphs[i]->SetLineWidth(2);
      graphs[i]->SetLineColor(kGraphColors[i]);
      graphs[i]->SetMarkerColor(kGraphColors[i]);
      graphs[i]->SetMarkerStyle(kGraphStyles[i]);
      fits[i] = new TF1(TString::Format("fitFn%d",i),"pol0");
      fits[i]->SetLineColor(kGraphColors[i]);
      graphs[i]->Fit(TString::Format("fitFn%d",i),"Q");
    }
    TCanvas *canvas2 = new TCanvas("canvas2","canvas2",400*2,400*2);
    canvas2->Divide(2,2);
    canvas2->cd(1);
    graphs[0]->Draw("AP");
    graphs[0]->SetTitle("Laser OFF Asymmetry;Entry Number;Asymmetry (%)");
    graphs[2]->Draw("SAME P");
    graphs[3]->Draw("SAME P");
    canvas2->cd(2);
    graphs[1]->Draw("AP");
    graphs[1]->SetTitle("Laser ON Asymmetry;Entry Number;Asymmetry (%)");
    canvas2->cd(3);
    graphs[4]->Draw("AP");
    graphs[4]->SetTitle("Bkg / Signal Fraction;Entry Number;Bk/Signal");
    canvas2->cd(4);
    TH2F *h2 = new TH2F("h2","BkSub Asymmetry;Entry Number;Asymmetry (%)",10,
        graphs[0]->GetXaxis()->GetXmin(),graphs[0]->GetXaxis()->GetXmax(),10,-15,15);
    h2->SetStats(kFALSE);
    h2->Draw();
    graphs[5]->Draw("SAME P");
    //graphs[5]->GetYaxis()->SetLimits(-15,15);
    graphs[5]->Draw("AP");
    graphs[5]->SetTitle("BkSub Asymmetry;Entry Number;Asymmetry (%)");
    saveHisto(canvas2,TString::Format("gAcc%d_%4d_asymmetries.png",ACC,gRun));
    delete canvas2;
  }

  // Do the same for the non-cycle based analyzis
  for(int m  = 0; m < 5; m++) {
    gAccNoCycleONResults[2*m] = hAccON[m]->GetMean();
    gAccNoCycleONResults[2*m+1] = hAccON[m]->GetRMS()/
      TMath::Sqrt(hAccON[m]->GetEntries());
    gAccNoCycleOFFResults[2*m] = hAccOFF[m]->GetMean();
    gAccNoCycleOFFResults[2*m+1] = hAccOFF[m]->GetRMS()/
      TMath::Sqrt(hAccOFF[m]->GetEntries());
  }
  backgroundSubtract(&gAccNoCycleONResults[4],&gAccNoCycleOFFResults[4],
      &gAccNoCycleONResults[8],&gAccNoCycleOFFResults[8],
      gAccNoCycleSubResults);
      */
  // Finally, clear out the histograms
  for(size_t v = 0; v < gVariablesOld.size(); v++) {
    for(size_t c = 0; c < gLaserCyclesStart[0].size(); c++) {
      for(size_t s = 0; s < NUM_STATES-1; s++) {
        for(size_t m = 0; m < NUM_MEASUREMENTS; m++) {
          delete histograms[v][c][s][m];
        }
        histograms[v][c][s].clear();
      }
      histograms[v][c].clear();
    }
    histograms[v].clear();
  }
  histograms.clear();
}

void configVariable(TString treeName, void *ptr)
{
  gChain->SetBranchStatus(treeName,1);
  gChain->SetBranchAddress(treeName,ptr);
}

void addVariable(TString treeName, TString printName, TString units, int print,
    int *ptr, bool norm, bool sub, bool pol = false)
{
  gVariablesOld.push_back(new ComptonVariableOld(printName,ptr,units,print,norm,sub,pol));
  configVariable(treeName,ptr);
}

void addVariable(TString treeName, TString printName, TString units, int print,
    float *ptr, bool norm, bool sub, bool pol = false)
{
  gVariablesOld.push_back(new ComptonVariableOld(printName,ptr,units,print,norm,sub,pol));
  configVariable(treeName,ptr);
}

void addVariable(TString treeName, TString printName, TString units, int print,
    double *ptr, bool norm, bool sub, bool pol = false)
{
  gVariablesOld.push_back(new ComptonVariableOld(printName,ptr,units,print,norm,sub,pol));
  configVariable(treeName,ptr);
}

void laserPatternWise(int run = 4430, bool readCyclesFromFile = false)
{
  gStyle->SetPaperSize(TStyle::kUSLetter);
  gStyle->SetOptStat("emr"); // remove 'n' which is the histogram name
  gStyle->SetPadTopMargin(0.25);
  gStyle->SetPadLeftMargin(150./kCanvSizeX);
  gStyle->SetPadRightMargin(0.1);
  gStyle->SetLabelSize(0.045,"XYZ");
  gStyle->SetTitleSize(0.05,"XYZ");
  gErrorIgnoreLevel = kWarning;

  if(!init(run,readCyclesFromFile))
    return;

  readVariableLimits();
  gAnalyzingPowerOld = getAnalyzingPower(gRun);
  // Now configure the variables we want to perform an analysis on
  //addVariable("bcm","BCM","uA",2,&bcm,false,false,false);
  //addVariable("Acc0","Acc0","rau",2,&Acc0,false,true,true);
  //addVariable("Acc0","Acc0BN","radc",1,&Acc0,true,true,true);
  //addVariable("Acc4","Acc4","rau",1,&Acc4,false,true,true);
  //addVariable("Acc4","Acc4BN","radc",1,&Acc4,true,true,true);
  //addVariable("runScaler0","DSBg1","",2,&runScaler0,false,false,false);
  //addVariable("runScaler1","DSBg2","",2,&runScaler1,false,false,false);
  //addVariable("runScaler4","USBg1","",2,&runScaler4,false,false,false);
  //addVariable("runScaler0","DSBg1BN","",2,&runScaler0,true,false,false);
  //addVariable("runScaler1","DSBg2BN","",2,&runScaler1,true,false,false);
  //addVariable("runScaler4","USBg1BN","",2,&runScaler4,true,false,false);
  ComptonVariable<double> varNAcc0("NAcc0",NAcc[0],"",0,0);
  ComptonVariable<double> varAcc0("Acc0",Acc[0],"srau",1,&varNAcc0);
  ComptonVariable<double> varAcc4("Acc4",Acc[4],"srau",1,0);
  gVariables.AddVariable(&varAcc0,&(gAnalyzingPower[0]));
  gVariables.AddVariable(&varAcc4,&(gAnalyzingPower[4]));

  //if(kSaveIntHistos) {
    //gSystem->mkdir(TString::Format("$COMPTON_WEB/Run%d",gRun),kTRUE);
  //}

  // Read pedestals from file
  //if(kReadPedestals)
  //  readPedestals();

  // Now walk the tree and start analysing
  //walkTreeAndAnalyze();
  walkTreeAndAnalyze();

  // Print out variable pattern polarizations

  // Now compute the overal laser based asymmetry
  std::fstream summaryOutFile;
  summaryOutFile.open(TString::Format("%s/laserwise_summary.txt",gComptonOutPath.Data()),std::ios::out);
  printSummary();
  // Now save the same summary to file
  printSummary(&(summaryOutFile));
  summaryOutFile.close();
  return;

  // Delete the ComptonVariableOlds
  for(size_t v = 0; v < gVariablesOld.size(); v++) {
    delete gVariablesOld[v];
  }
  gVariablesOld.clear();


  printLine();
  printResult("NO Cycle  BkSub Asym ",gAccNoCycleSubResults[2],
      gAccNoCycleSubResults[3],"%");
  printResult("NO Cycle  BkSub Pol  ",100*gAccNoCycleSubResults[2]/gAnalyzingPowerOld,
      100*gAccNoCycleSubResults[3]/gAnalyzingPowerOld,"%");
  printLine();
  // Now print out the results
  //printAsymmetries();
}
