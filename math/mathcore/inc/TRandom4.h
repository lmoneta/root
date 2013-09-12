// @(#)root/mathcore:$Id: TRandom4.h $
// Author: Pierre Schweitzer 20/09/2012

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRandom4
#define ROOT_TRandom4

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRandom4                                                             //
//                                                                      //
// random number generator class (periodicity > 2**128)                 //
//  based on Random123                                                  //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TRandom
#include "TRandom.h"
#endif

namespace ROOT {
   namespace Math { 
      class Random123; 
   }
}

class TRandom4 : public TRandom {

private:

   ROOT::Math::Random123 * fRandom123;

public:
   TRandom4(UInt_t seed=4357);
   virtual ~TRandom4();
   // get the current seed (only first element of the seed table)
   virtual  UInt_t    GetSeed() const;
   virtual  Double_t  Rndm(Int_t i=0);
   virtual  void      RndmArray(Int_t n, Float_t *array);
   virtual  void      RndmArray(Int_t n, Double_t *array);
   virtual  void      SetSeed(UInt_t seed=0);

   ClassDef(TRandom4,1)  //Random number generator: Threefry
};

R__EXTERN TRandom *gRandom;

#endif
