// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Aug 4 2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// random engines based on ROOT 

#ifndef ROOT_Math_RanLuxPPEngine
#define ROOT_Math_RanLuxPPEngine

#ifdef USE_RANLUXPP

#include "Math/TRandomEngine.h"

#include <cstdint>
#include <vector>
#include <string>

class ranluxpp;

namespace ROOT {

   namespace Math {

      /**

         Random engine baded on version 3.3 of RanLux
         see ..


         @ingroup Random
      */

      


      /**
         RanLuxPP engine
       */
     class RanLuxPPEngine : public TRandomEngine {


      public:

         typedef  TRandomEngine BaseType;
         typedef  uint64_t Result_t;
         typedef  uint64_t StateInt_t;


	 RanLuxPPEngine(uint64_t seed=4357, int level=-1);

         virtual ~RanLuxPPEngine();

         void SetSeed(Result_t seed) {
            Init(seed);
         }

         virtual double Rndm() {
	   return Rndm_impl();
         }
         inline double operator() () { return Rndm_impl(); }

         uint32_t IntRndm() {
            return IntRndm_impl();
         }
       
         /// minimum integer taht can be generated
         static unsigned int MinInt() { return 0; }
         /// maximum integer taht can be generated
         static uint64_t MaxInt() { return 0xffffffffffff; }  //  2^24 -1

         static int Size() { return 9; }

         
         static std::string Name() {
            return "RanLuxPPEngine";
         }

      protected:
         // functions used for testing
         
         void SetState(const std::vector<uint64_t> & state);

         void GetState(std::vector<uint64_t> & state);

         void Init(uint64_t seed);

         double Rndm_impl();
	 
         float  Rndm_impl_float();

      private:

         uint32_t IntRndm_impl() {return (int) Rndm_impl_float()*MaxInt(); }

         /* int fLevel; */
         /* double fX[1]; */

	 ranluxpp * fRlxpp;


      };

     class RanLuxPPEngineFloat : public RanLuxPPEngine {

     public:
       
       RanLuxPPEngineFloat (uint64_t seed=4357, int level=-1) :
       RanLuxPPEngine(seed,level)
	 {}

       virtual ~RanLuxPPEngineFloat() {}

       virtual double Rndm() {
            return Rndm_impl_float();
         }
       inline double operator() () { return Rndm_impl_float(); }

     };


   } // end namespace Math

} // end namespace ROOT

#endif  /* USE_RANLUXPP */
#endif /* ROOT_Math_RanluxPPEngine */
