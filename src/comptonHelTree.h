/*
 * comptonVariable Class
 * 
 * A single place to keep all the helicity correlated variables which can be pushed
 * to either the mpstree or the multiplettree
*/

#ifndef comptonHelTree_h
#define comptonHelTree_h

#include <vector>

#define HELPLUS 1
#define HELMINUS 0


class TTree;

struct vcomptonVariable {
  virtual void ClearMultiplet() = 0;
  virtual void ProcessHelicity(int hel) = 0;
  virtual void Norm(double) = 0;
};

template<typename T>
struct comptonVariable : public vcomptonVariable {
  T mpsval; ///< MPS wise variable
  T multval[2]; ///< Multiplet wise variable
  virtual void ClearMultiplet() {
    multval[0] = multval[1] = 0;
  }
  virtual void ProcessHelicity(int hel) {
    multval[hel] += mpsval;
  }
  virtual void Norm(double norm) {
    multval[0] *= norm;
    multval[1] *= norm;
  }
};

// This class will contain the helicity correlated variables that we are storing
// in the mpswise and quartetwise trees.
class comptonHelTree {
public:
  comptonHelTree(TTree *mps = 0, TTree *mult = 0);
  ~comptonHelTree();

  template<typename T>
  void AddVariable(comptonVariable<T> *var, const char *vname, const char *vnamemult = 0, bool davg = false);

  TTree* TreeMPS() { return fTreeMPS; }
  TTree* TreeMult() { return fTreeMult; }
  void ClearMultiplet();
  void ProcessHelicity(int hel);

  int MultStatus();

  int CountHel(int h) {
    return fCountHel[h];
  }
  void SetHelStructure(int hel);
  void ProcessFullMult();

  // Status for multiplet
  static const int kMultOK;
  static const int kMultUnbalanced;
  static const int kMultBadCount;

private:
  TTree *fTreeMPS;
  TTree *fTreeMult;
  double fNorm;
  std::vector<vcomptonVariable*> fVars;
  std::vector<bool> fDoAvg;
  int fCountHel[2];
  int fHelStructure;
};


template<typename T>
void DefineComptonVarMultipletBranch(comptonHelTree *tree, comptonVariable<T> *var,
    const char *vname);

template<typename T>
void DefineComptonVarMPSBranch(comptonHelTree *tree, comptonVariable<T> *var,
    const char *vname);


#endif


