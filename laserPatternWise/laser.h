#ifndef LASER_H
#define LASER_H

#include <vector>
#include <iostream>

typedef struct {
  double coordinateValue;
  double coordinateError;
  double result;
  double error;
} LaserResult_t;

typedef struct {
  int laserOn;
  int start;
  int end;
  int mpsCount;
  int beamOnCount;
  int beamOffCount;
  int beamOnStatus;
  int beamOffStatus;
  int cycleStatus;
  int start_entry;
  int end_entry;
} LaserCycle_t;

// Define a single instance of Laser Cycles
std::vector<LaserCycle_t> gLaserCycles;
//LaserCycles gLaserCycles;

class LaserCycleResults {
public:
  LaserCycleResults() : fInitialized(false) {};
  virtual ~LaserCycleResults() {};

  void init() {
    for(int s = 0; s < 4; s++) {
      for(size_t c = 0; c < gLaserCycles.size(); c++) {
        fResults[s].push_back(0.0);
      }
    }
    fInitialized = true;
  }
  std::vector<double>& operator[](int var) {
    if(!fInitialized)
      init();
    return fResults[var];
  }

  size_t size() { return fResults[0].size(); }
private:
  std::vector<double> fResults[4];
  bool fInitialized;
};

/*
template<class T>
class LaserCycleObjects {
public:
  LaserCycleObjects() {
    fObjects.resize(gLaserCycles.size());
  };
  virtual ~LaserCycleObjects() {};
private:
  std::vector<T> fObjects;
};
*/

// Structure defining laser pattern
typedef struct {
  LaserCycle_t *offLeft;
  LaserCycle_t *on;
  LaserCycle_t *offRight;
  int offLeftNum;
  int onNum;
  int offRightNum;
  int status;
  int patNum;
  int mpsIsInState(int mpsCoda, int half = 0) {
    if( mpsCoda>= (offLeft->start + half*offLeft->mpsCount/2)
        && mpsCoda <= offLeft->end) {
      return 0;
    } else if( mpsCoda>= on->start && mpsCoda <= on->end) {
      return 1;
    } else if( mpsCoda>= offRight->start && mpsCoda <=
        (offRight->end - half*offRight->mpsCount/2.) ) {
      return 2;
    } else {
      return -1;
    }
  }
  int entryIsInState(int entry) {
    if( entry>= offLeft->start_entry && entry <= offLeft->end_entry) {
      return 0;
    } else if( entry>= on->start_entry && entry <= on->end_entry) {
      return 1;
    } else if( entry>= offRight->start_entry && entry <= offRight->end_entry) {
      return 2;
    } else {
      return -1;
    }
  }

} LaserPattern_t;

// Single instance of all laser patterns for beam on and off
std::vector<LaserPattern_t> gLaserPatterns;
std::vector<LaserPattern_t> gLaserPatternsBeamOff;

const char *kLaserPatternObjectNames[5] = { "OffLeft", "On", "OffRight", "Background",
  "BackgroundSubtracted" };

template<class T>
class LaserPatternObject {
public:
  LaserPatternObject() {};
  virtual ~LaserPatternObject() {};

  T& operator[](size_t state) {
    switch(state) {
      case 0:
        return fOffLeft;
      case 1:
        return fOn;
      case 2:
        return fOffRight;
      case 3:
        return fBackground;
      case 4:
        return fBkgSubtracted;
      default:
        return fNull;
    }
  }

  T& getBkg() { return fBackground; }
  T& getOn() { return fOn; }
  T& getSub() { return fBkgSubtracted; }
  static int size() { return 5; }
private:
  T fOffLeft;
  T fOn;
  T fOffRight;
  T fBackground;
  T fBkgSubtracted;
  T fNull;
};

template<class T>
class LaserPatternObjectList {
public:
  LaserPatternObjectList(std::vector<LaserPattern_t> *patterns = 0) :
    fInitialized(0) {
    fPatterns = patterns;
    init();
  };
  virtual ~LaserPatternObjectList() {};

  LaserPatternObject<T>& operator[] (size_t pat) {
    return fPatObjects[pat];
  }

  std::vector<T*> byCycle(int cycle) {
    if(cycle<0) {
      std::vector<T*> tmp;
      return tmp;
    }
    return fCycleToPatObjects[cycle];
  }

  size_t size() {
    return fPatterns->size();
  }
private:
  std::vector<LaserPattern_t> *fPatterns;
  std::vector<LaserPatternObject<T> > fPatObjects;
  std::vector<std::vector<T*> > fCycleToPatObjects;
  bool fInitialized;
  void init()
  {
    if(!fPatterns)
      return;
    fPatObjects.resize(fPatterns->size());
    fCycleToPatObjects.resize(gLaserCycles.size());
    for(size_t c = 0; c < gLaserCycles.size(); c++) {
      for(size_t p = 0; p < fPatterns->size(); p++) {
        if(fPatterns->at(p).offLeftNum == int(c)) {
          fCycleToPatObjects[c].push_back(&fPatObjects[p][0]);
          fCycleToPatObjects[c].push_back(&fPatObjects[p][3]);
        } else if (fPatterns->at(p).onNum == int(c)) {
          fCycleToPatObjects[c].push_back(&fPatObjects[p][1]);
        } else if (fPatterns->at(p).offRightNum == int(c)) {
          fCycleToPatObjects[c].push_back(&fPatObjects[p][2]);
          fCycleToPatObjects[c].push_back(&fPatObjects[p][3]);
        }
      }
    }
    fInitialized = 1;
  }
};

typedef struct {
  double patNum;
  double value;
  double error;
  double rms;
  double entries;
} LaserPatternResult_t;

class LaserPatternResults :
  public LaserPatternObjectList<LaserPatternResult_t> {
public:
  LaserPatternResults(std::vector<LaserPattern_t> *pattern) :
    LaserPatternObjectList<LaserPatternResult_t>(pattern) {
  };
  virtual ~LaserPatternResults() {};
private:
};

#endif // LASER_H
