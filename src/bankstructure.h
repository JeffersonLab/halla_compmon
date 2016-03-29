//  Creates structure used for CODA subbank info
//   pointer to Data explicitly set to uint32_t 2/11/2015  gbf
#include <stdint.h>
#ifndef bankstructure_h
#define bankstructure_h

struct bankstructure
{
  int Tag;     //bank tag
  int Header;  //bank header word
  int DataWords;  //# datawords following header
  int BankIndex;  //index of word zero of bank within event buffer
  uint32_t* Data;      //pointer to data
};

#endif
