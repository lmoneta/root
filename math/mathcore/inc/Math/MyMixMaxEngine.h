// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Aug 4 2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// random engines based on ROOT 

#ifndef ROOT_Math_MyMixMaxEngine
#define ROOT_Math_MyMixMaxEngine

#include <cstdint>
#include <vector>

#ifndef ROOT_Math_TRandomEngine
#include "Math/TRandomEngine.h"
#endif




namespace ROOT {

   namespace Math {
      

      /**
         My impelmentation of MixMax
         
         @ingroup Random
      */

      class MyMixMaxEngine : public TRandomEngine {


      public:

         typedef  TRandomEngine BaseType; 
         
         MyMixMaxEngine(uint64_t seed=1, int n = 0);

         virtual ~MyMixMaxEngine() {}

         /// get the state of the generator
         void GetState(std::vector<uint64_t> & state) const {
            state = fState; 
         }

         /// Get the counter (between 0 and Size-1)
         int Counter() const { return fCounter; }

         /// Get the size of the generator
         static int Size();
         static void SetSize(int n);
         static void SetSpecial(int s);

         /// set the generator seed 
         void  SetSeed(unsigned int ) {}

         /// set the generator seed using a 64 bits integer
         void SetSeed64(uint64_t ) {}

         ///set the full initial generator state and warm up generator by doing some iterations
         void SetState(const std::vector<uint64_t> & state, bool warmup = true) {
            fState = state;
            fCounter = fState.size(); 
         }

         // generate a random number (virtual interface)
         virtual double Rndm() { return Rndm_impl(); }

         /// generate a double random number (faster interface)
         inline double operator() () { return Rndm_impl(); }

         /// generate an array of random numbers 
         void RndmArray (int n, double * array) {}

         /// maximum integer that can be generated. For MIXMAX is 2^61-1         
         static uint64_t MaxInt() { return  0x1fffffffffffffff; } //  2^61 -1 

         /// generate a 64  bit integer number
         uint64_t IntRndm();

      private:

         /// implementation function to generrate the random number
         double Rndm_impl() {
            //#define INV_MERSBASE (0x1p-61)
            return ( (double)IntRndm() ) *(0x1p-61) ; 
         }

         int fCounter; 
         std::vector<uint64_t> fState;  // mix-max generator state
         std::vector<uint64_t> fMatrix;              // Mix Max matrix 
         uint64_t              fSumTot;              // total sum
         
      };


   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_TRandomEngines */
