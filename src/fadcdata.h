// fadcdata.h
// used to unpack and access FADC data
//
#include <stdlib.h>
#include <stdint.h>
#include "pulsestructure.h"
#include "textParams.h"
#include "vmeauxdata.h"
#include "THaCodaFile.h"

#ifndef fadcdata_h
#define fadcdata_h

#define MAX_FADC_CHANNELS 8
#define MAX_FADC_EVENTS 5000
#define MAX_UNPACKED_SUMS 5000
class fadcdata {
  private:
  textParams* theParams; //pointer to instance of parameter Class
  //    double PedValue[MAX_FADC_CHANNELS];
    // fadc sample info
    bool SamplesValid[MAX_FADC_CHANNELS];
    bool AccumValid[MAX_FADC_CHANNELS];
    bool SumsValid[MAX_FADC_CHANNELS];
    int* RawSamples[MAX_FADC_CHANNELS];
    int* UserBits[MAX_FADC_CHANNELS];
    int NumberSamples[MAX_FADC_CHANNELS];
    int triggermode;
    int RandomTimes[MAX_FADC_CHANNELS];;
    int sampleformat;
    int NumberEvents[MAX_FADC_CHANNELS];
    int SamplesPerEvent[MAX_FADC_CHANNELS];
    int samplepointer[MAX_FADC_CHANNELS];
    int  UnpackedClock[MAX_FADC_CHANNELS][MAX_FADC_EVENTS];
    //
    // stuff from  "sums" subbank
    int Sums_InputRegister;
    int Sums_NumberFADCChannels;
    int Sums_PulserIndex;
    int Sums_NumberFADCTriggers[MAX_FADC_CHANNELS];
    int Sums_NumberTriggersSummed[MAX_FADC_CHANNELS];
    int Sums_NumberSamplesSummed[MAX_FADC_CHANNELS];
    int Sums_PedestalSubtracted[MAX_FADC_CHANNELS];
    int Sums_Data[MAX_FADC_CHANNELS][MAX_UNPACKED_SUMS];
    int Sums_Clock[MAX_UNPACKED_SUMS];
    int Sums_PulserSynch[MAX_UNPACKED_SUMS];
    //
    // stuff from accumulators subbank

    int64_t AccumValue[MAX_FADC_CHANNELS][6];
    int AccumNumSamples[MAX_FADC_CHANNELS][6];
    //parameters
    int DacSetting[MAX_FADC_CHANNELS];
    int AccumThresh[MAX_FADC_CHANNELS][2];
    int N4before[MAX_FADC_CHANNELS];
    int N4after[MAX_FADC_CHANNELS];
    int N5before[MAX_FADC_CHANNELS];
    int N5after[MAX_FADC_CHANNELS];
    int N6before[MAX_FADC_CHANNELS];
    int N6after[MAX_FADC_CHANNELS];
    pulsestructure CurrentPulse[MAX_FADC_CHANNELS];

    int PulseIndex[MAX_FADC_CHANNELS][MAX_FADC_EVENTS];
    int TimerPar1;  //HAPPEX timer board parameter 1
    int TimerPar2;  //parameter 2
    int TimerDac[2];  // DAC1 and DAC2 settings

  public:
    fadcdata();
    fadcdata(textParams* theParamsIn);
    int newMPS();   //called at start of MPS to clear data
    int UnpackAllSubbanks(THaCodaFile* codaData, vmeauxdata* theVMEauxdata);
    int UnpackSamples(bankstructure bank);
    int UnpackAccumulators(bankstructure bank);
    int UnpackTimer(bankstructure bank);
    int UnpackSums(bankstructure bank, int verbose, int abortOnError);
    void DumpBank(bankstructure bank);
    //
    //sums info
    int GetSumsNumberADCChannels(){
      return Sums_NumberFADCChannels;}; //# FADC chans
    int GetSumsPulserIndex(){
      return Sums_PulserIndex;};
    int GetSumsNumberTriggers(int channel){
      return Sums_NumberFADCTriggers[channel];};
    int GetSumsNumberTriggersSummed(int channel){
      return Sums_NumberTriggersSummed[channel];}; //# triggers actually stored
    int GetNumberSamplesSummed(int channel) {
      return Sums_NumberSamplesSummed[channel];};
    int GetSumsPedestalSubtracted(int channel){
      return Sums_PedestalSubtracted[channel];};
    int GetSums(int channel,int trigger){
      return Sums_Data[channel][trigger];}; //pedestal corrected sum
    int GetSumsClock(int channel, int trigger){
      return Sums_Clock[trigger];}; //return clock time
    int GetSumsPulserSynch(int channel, int trigger){
      return Sums_PulserSynch[trigger];}; //return clock time
    //
    //sample info
    int* GetSamples(int channel){
      return RawSamples[channel]; };  //return pointer to array of FADC samples
    int* GetUserBits(int channel){
      return UserBits[channel]; };    //return pointer to array of user bit input

    pulsestructure GetPulse(int chan,int pulsenumber);

    int GetNumberSamples(int channel) {
      return NumberSamples[channel]; };  //Return number of samples
    int GetNumberEvents(int channel) {
      return NumberEvents[channel]; }; //return number of events for multevent mode
    int GetSamplesPerEvent(int channel){
      return SamplesPerEvent[channel]; };//return number samples per event
    int GetSamplePointer(int channel) {
      return samplepointer[channel] ; } ; //return sample index from FADC
    //accumulator info
    int GetAccumNumSamples(int chan, int accum){
      return AccumNumSamples[chan][accum]; };
    uint64_t GetAccumValue(int chan, int accum){
      return AccumValue[chan][accum]; };
    //Get Settings
    int GetDac(int chan){
      return DacSetting[chan]; };
    int GetThresh(int chan, int which) {		// near = 1, far = 2
      return AccumThresh[chan][which] ; } ;
    int GetN4before(int chan){
      return N4before[chan]; };
    int GetN4after(int chan){
      return N4after[chan]; };
    int GetN5before(int chan){
      return N5before[chan]; };
    int GetN5after(int chan){
      return N5after[chan]; };

    bool IsSamplesValid(int channel){
      return SamplesValid[channel]; };   //true if Samples data have been unpacked
    bool IsAccumValid(int channel) {
      return AccumValid[channel]; };    //true if Accumulator data unpacked
    bool IsSumsValid(int channel) {
      return SumsValid[channel]; };    //true if Sums data unpacked
    bool IsRandomTimes(int channel) {
      return RandomTimes[channel]; };  //true if random time sample

    int GetTimerPar1() {
      return TimerPar1; } ;
    int GetTimerPar2() {
      return TimerPar2; } ;
    int GetPulserSetting(int flag) {
      return TimerDac[flag]; } ;
};

#endif
