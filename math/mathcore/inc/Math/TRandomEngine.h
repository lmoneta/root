// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Aug 4 2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// random engines based on ROOT 

#ifndef ROOT_Math_TRandomEngine
#define ROOT_Math_TRandomEngine

#include <cstdint>
#include <vector>
#include <cassert>

namespace ROOT {

   namespace Math {

      class RandomBaseEngine {
      public:
         virtual double Rndm() = 0;
         virtual ~RandomBaseEngine() {}
      };


      class TRandomEngine : public RandomBaseEngine {
      public:
         virtual ~TRandomEngine() {}
      };

      class LCGEngine : public TRandomEngine {


      public:

         typedef  TRandomEngine BaseType;
         typedef  uint32_t result_t;
         typedef  uint32_t StateInt_t;

         LCGEngine() : fSeed(65539) { }

         virtual ~LCGEngine() {}

         void SetSeed(uint32_t seed) { fSeed = seed; }

         virtual double Rndm() {
            //double Rndm() {
            return Rndm_impl();
         }
         inline double operator() () { return Rndm_impl(); }

         uint32_t IntRndm() {
            fSeed = (1103515245 * fSeed + 12345) & 0x7fffffffUL;
            return fSeed; 
         }

         /// minimum integer taht can be generated
         static unsigned int MinInt() { return 0; }
         /// maximum integer taht can be generated
         static unsigned int MaxInt() { return 0xffffffff; }  //  2^32 -1
         /// Size of the generator state
         static int Size() { return 1; }

      protected:
         // for testing all generators
         void SetState(const std::vector<uint32_t> & state) {
            assert(!state.empty());
            fSeed = state[0]; 
         }

         void GetState(std::vector<uint32_t> & state) {
            state.resize(1);
            state[0] = fSeed;
         }
         int Counter() const { return 0; }
      private:

         double Rndm_impl() {
            const double kCONS = 4.6566128730774E-10; // (1/pow(2,31)
            unsigned int rndm = IntRndm(); // generate integer number 
            if (rndm != 0) return  kCONS*rndm;
            return Rndm_impl();
         }
         
         uint32_t fSeed;
      };
      

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_TRandomEngine */
