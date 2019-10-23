#ifndef MYMATH_H
#define MYMATH_H

#include <TMath.h>
#include <TH1.h>
#include <iostream>
#include <iomanip>

const double kStandardErrorValue = -4e-6;

class MyData {
public:
  MyData(double val = 0., double err = kStandardErrorValue) :
    fValue(val), fError(err) {};

  MyData& operator=(const MyData &rhs) {
    this->fValue = rhs.fValue;
    this->fError = rhs.fError;
    return *this;
  }

  MyData& operator+=(const MyData &rhs) {
    fValue += rhs.fValue;
    fError = TMath::Sqrt( fError*fError + rhs.fError*rhs.fError);
    return *this;
  }

  MyData& operator/=(const MyData &rhs) {
    this->fError = TMath::Sqrt( TMath::Power(this->fError/this->fValue,2.0) +
        TMath::Power(rhs.fError/rhs.fValue,2) );
    this->fValue /= rhs.fValue;
    this->fError *= this->fValue;
    return *this;
  }
  MyData& operator*=(double mult ) {
    this->fValue *= mult;
    this->fError *= mult;
    return *this;
  }

  friend MyData operator+(MyData lhs, const MyData &rhs) {
    lhs += rhs;
    return lhs;
  }

  MyData& operator-=(const MyData &rhs) {
    fValue -= rhs.fValue;
    fError = TMath::Sqrt( fError*fError + rhs.fError*rhs.fError);
    return *this;
  }

  friend MyData operator-(MyData lhs, const MyData &rhs) {
    lhs -= rhs;
    return lhs;
  }

  friend MyData operator/(MyData lhs, const MyData &rhs) {
    lhs /= rhs;
    return lhs;
  }

  friend MyData operator*(MyData lhs, double mult) {
    lhs *= mult;
    return lhs;
  }

  friend MyData operator*(double mult, MyData lhs) {
    lhs *= mult;
    return lhs;
  }



  friend std::ostream& operator<<(std::ostream& out, const MyData& d) {
    return out << std::fixed
      << std::setprecision(5) << std::setw(8) << d.fValue
      << " +/- "
      << std::setprecision(5) << std::setw(8) << d.fError;
  }

  void Set(double val) { fValue = val; }
  void SetErr(double err) { fError = err; }
  void Set(double val, double err) {
    Set(val);
    SetErr(err);
  }

  virtual void Set(MyData dat) {
    MyData::Set(dat.val(),dat.err());
  }

  double val() { return fValue; }
  double err() { return fError; }

private:
  double fValue;
  double fError;
};

class MyHistoData : public MyData {
public:
  MyHistoData() : MyData(0.,0.), fRMS(0), fEntries(0) {
  }
  void Set(TH1* h) {
    MyData::Set(h->GetMean(),h->GetMeanError());
    fRMS = h->GetRMS();
    fEntries = h->GetEntries();
  };

  void Set(MyData dat) {
    MyData::Set(dat);
    fRMS = fEntries = 0;
  }

private:
  double fRMS;
  int fEntries;
};

class MyErrorWeightedData {
public:
  MyErrorWeightedData() : fVal(0.), fErr(0.) {}
  void add(MyData dat) {
    fVal += dat.val()/(dat.err()*dat.err());
    fErr += 1./(dat.err()*dat.err());
  }
  void finish() {
    fVal /= fErr;
    fErr = 1./TMath::Sqrt(fErr);
  }
  double val() { return fVal; }
  double err() { return fErr; }
private:
  double fVal;
  double fErr;
};

#endif // MYMATH_H
