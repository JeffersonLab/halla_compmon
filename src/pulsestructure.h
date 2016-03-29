#ifndef pulsestructure_h
#define pulsestructure_h

struct pulsestructure
{
  bool Valid;       //true if valid pulse
  int NumberSamples;  //# fadc samples in waveform
  int PulseNumber;   //number of thispulse withing coda event
                     //0 to NumberSamples-1
  int Index;      //index of word zero in fadc data
  int Clock;           //timing info
  int Sum;
  int RandomTimes;   //=1 if data is random samples 9/14/2004
  int* Data;          //pointer to fadc data
  int* UserBits;      //pointer to fadc user bit data
};

#endif
