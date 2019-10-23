#include "comptonHelTree.h"
#include <iostream>

#include <TTree.h>
#define IMP_MAKECOMPVARTYPE(MYTYPE)\
  template void comptonHelTree::AddVariable<MYTYPE>(comptonVariable<MYTYPE>*, const char*vname,\
     const char *vmultname, bool avg);

// Status for multiplet
const int comptonHelTree::kMultOK = 0;
const int comptonHelTree::kMultUnbalanced = 1;
const int comptonHelTree::kMultBadCount = 1<<2;

template<typename T>
void DefineComptonVarBranches(comptonHelTree *tree, comptonVariable<T> *var,
  const char *vname, const char *vmultname)
{
  tree->TreeMPS()->Branch(vname,&(var->mpsval));
  tree->TreeMult()->Branch(TString::Format("PosHel%s",vmultname),&(var->multval[HELPLUS]));
  tree->TreeMult()->Branch(TString::Format("NegHel%s",vmultname),&(var->multval[HELMINUS]));
  (*var).mpsval = 0;
  (*var).multval[0] = 0;
  (*var).multval[1] = 0;
}

IMP_MAKECOMPVARTYPE(int)
IMP_MAKECOMPVARTYPE(float)
IMP_MAKECOMPVARTYPE(double)

comptonHelTree::comptonHelTree(TTree *mps, TTree *mult) : fTreeMPS(mps), fTreeMult(mult) {
}

void comptonHelTree::ClearMultiplet()
{
  for(std::vector<vcomptonVariable*>::iterator it = fVars.begin();
    it != fVars.end(); it++) {
    (*it)->ClearMultiplet();
  }
  fCountHel[0] = fCountHel[1] = 0;
}

void comptonHelTree::ProcessHelicity(int hel)
{
  for(std::vector<vcomptonVariable*>::iterator it = fVars.begin();
    it != fVars.end(); it++) {
    (*it)->ProcessHelicity(hel);
  }
  fCountHel[hel]++;
}

template<typename T>
void comptonHelTree::AddVariable(comptonVariable<T> *var, const char *vname,
    const char *vmultname, bool avg)
{
  fVars.push_back(var);
  fDoAvg.push_back(avg);
  DefineComptonVarBranches<T>(this,var,vname,vmultname);
}

int comptonHelTree::MultStatus()
{
  int status = kMultOK;
  if(fCountHel[0] != fCountHel[1]) // Same number of helicities
    status |= kMultUnbalanced;
  if(fCountHel[0]+fCountHel[1] != fHelStructure)  // The expected number of helicitie states in this multiplet
    status |= kMultBadCount;

  return status;
}


void comptonHelTree::SetHelStructure(int hel)
{
  fHelStructure = hel;
  fNorm = 2./double(fHelStructure);
}

void comptonHelTree::ProcessFullMult()
{
  // Normalize the helicity structure for those variables that need it
  for(size_t v = 0; v < fVars.size(); v++) {
    if(fDoAvg[v]) {
      fVars[v]->Norm(fNorm);
    }
  }
}
