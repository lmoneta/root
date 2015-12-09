// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Aug 4 2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// random engines based on ROOT 

#ifndef ROOT_Math_RanLuxEngine
#define ROOT_Math_RanLuxEngine

#ifndef ROOT_Math_TRandomEngine
#include "Math/TRandomEngine.h"
#endif

#include "TRandom1.h"

#include <stdint.h>

namespace ROOT {

   namespace Math {

      /**
         Random number generator class based on RanLux
         
         @ingroup Random
      */
      
      class RanLuxEngine : public TRandomEngine {


      public:

         typedef  TRandomEngine BaseType; 
         
         RanLuxEngine(unsigned int seed=111, int lux = 3) : fRndm(seed,lux)  {
            //SetSeed(seed);
         }

         virtual ~RanLuxEngine() {}

         void SetSeed(unsigned int seed)  { fRndm.SetSeed(seed); }

         static int Size() { return 24; }

         int Counter() const { return fRndm.fCount24; }

         // set a state from a vector of 32 bit floats
         void SetState(const std::vector<float> & state) {
            unsigned int n = Size(); 
            for (unsigned int i = 0; i < n; ++i)
               fRndm.fFloatSeedTable[i] = state[i];
            fRndm.fCount24 = n; // to make sure we re-iterate on the new state
         }
         // set state from a vectior of 32 bits integers
         void SetState(const std::vector<unsigned int> & state) {
            unsigned int n = Size(); 
            for (unsigned int i = 0; i < n; ++i)
               fRndm.fFloatSeedTable[i] = *((float*) &state[i]);
            fRndm.fCount24 = n; // to make sure we re-iterate on the new state
         }
         void GetState(std::vector<float> & state) const {
            state.resize(Size());
            for (unsigned int i = 0; i < state.size(); ++i)
               state[i] = fRndm.fFloatSeedTable[i];
         }
         void GetState(std::vector<unsigned int> & state) const {
            state.resize(Size());
            for (unsigned int i = 0; i < state.size(); ++i)
               state[i] = *( (unsigned int*) &fRndm.fFloatSeedTable[i]);
         }


         virtual double Rndm() {
            return Rndm_impl();
         }
         inline double operator() () { return Rndm_impl(); }

         static unsigned int MaxInt() { return 0xffffff; } // 2^24 - 1 
         
         unsigned int IntRndm() {
            const int imax = std::pow(2,23); 
            return fRndm.Integer( imax );
         }

      private:

         double Rndm_impl() { return fRndm.Rndm(); }

         
         TRandom1 fRndm;

      };
      

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_TRandomEngines */
