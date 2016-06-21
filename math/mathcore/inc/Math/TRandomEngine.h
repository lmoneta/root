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

#include <vector>



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
         
         LCGEngine() : fSeed(65539) { }

         virtual ~LCGEngine() {}

         static int Size() { return 1; }
         int Counter() const { return 0; } // not used

         void SetState(const std::vector<unsigned int> & state) { SetSeed(state[0]); }

         void GetState(std::vector<unsigned int> & state) const {state[0] = fSeed;  }

         void SetSeed(unsigned int seed) { fSeed = seed; }
         
         virtual double Rndm() {
            //double Rndm() {
            return Rndm_impl();
         }
         double Rndm_impl() { 
            const double kCONS = 4.6566128730774E-10; // (1/pow(2,31)
            unsigned int rndm = IntRndm(); // generate integer number 
            if (rndm != 0) return  kCONS*rndm;
            return Rndm_impl();
         }
         inline double operator() () { return Rndm_impl(); }

         /// maximum integer taht can be generated
         static unsigned int MaxInt() {
            return 0xffffffff; // 2^32 -1
         }
         
         unsigned int IntRndm() {
            fSeed = (1103515245 * fSeed + 12345) & 0x7fffffffUL;
            return fSeed; 
         }

      private:
         unsigned int fSeed; 
      };
      

   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_TRandomEngine */
