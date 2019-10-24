void jc2LinOld(int nparams = 2,int nsteps=147, double dark_delta = 0.003161)
{
  gROOT->ProcessLine(".x plotPulser.C");
  //gROOT->ProcessLine(Form(".x fitPulser.C(%d,150-3,0.003161)",nparams));
  gROOT->ProcessLine(Form(".x fitPulser.C(%d,%d,%g)",nparams,nsteps,dark_delta));
}


// 2019-10-17: New version that allows us to start with index > 0
//void jc2Lin(int nparams = 2, int firstIndex = 3, int lastIndex=150, double dark_delta = 0.003161, bool runPlot = true)
//void jc2Lin(int run, int halfGain = 0, int nparams = 2, int firstIndex = 0, int lastIndex=150, double dark_delta = 0.002539, bool runPlot = true)
void jc2Lin(int run, int halfGain = 0, int nparams = 2, int firstIndex = 0, int lastIndex=150, double dark_delta = 0.002612, bool runPlot = true)
{
  if(runPlot) {
    gROOT->ProcessLine(Form(".x plotPulser.C(%d)",run));
  }
  //gROOT->ProcessLine(Form(".x fitPulser.C(%d,150-3,0.003161)",nparams));
  gROOT->ProcessLine(Form(".x fitPulser.C(%d,%d,%d,%g,%d)",nparams,firstIndex,lastIndex,dark_delta,halfGain!=0?1:0));
}
