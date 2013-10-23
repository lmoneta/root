// @(#)root/mathcore:$Id: TRandom4.h $
// Author: Lorenzo Moneta Sep/2013

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TRandom5
#define ROOT_TRandom5

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TRandom5                                                             //
//                                                                      //
// random number generator class                                        //
//  based on MixMax                                                     //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TRandom
#include "TRandom.h"
#endif

namespace ROOT {
   namespace Math { 
      class RandomMixMax; 
   }
}

class TRandom5 : public TRandom {

private:

   ROOT::Math::RandomMixMax * fRng;

public:
   TRandom5(UInt_t seed=1);
   virtual ~TRandom5();
   // get the current seed (only first element of the seed table)
   virtual  UInt_t    GetSeed() const;
   virtual  Double_t  Rndm(Int_t i=0);
   virtual  void      RndmArray(Int_t n, Float_t *array);
   virtual  void      RndmArray(Int_t n, Double_t *array);
   virtual  void      SetSeed(UInt_t   seed=1);
   virtual  void      SetSeed64(ULong64_t seed=1);
   virtual	void 	  SeedUniqueStream(UInt_t clusterID, UInt_t machineID, UInt_t runID, UInt_t  streamID);

   ClassDef(TRandom5,1)  //Random number generator: Threefry
};

R__EXTERN TRandom *gRandom;

#endif
