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

#include "Math/TRandomEngine.h"

#include <cstdint>
#include <vector>
#include <string>

class RanluxS;
class RanluxD;
namespace ROOT {

   namespace Math {

      /**

         Random engine baded on version 3.3 of RanLux
         see ..


         @ingroup Random
      */

      

      class RanLuxSEngine : public TRandomEngine {


      public:

         typedef  TRandomEngine BaseType;
         typedef  uint32_t Result_t;
         typedef  uint32_t StateInt_t;


         RanLuxSEngine(uint32_t level, uint32_t seed=4357) :
            fRlxs(0)
         {
            fLevel = level;
            Init(level,seed);
         }

         virtual ~RanLuxSEngine();

         void SetSeed(Result_t seed) {
            Init(fLevel,seed);
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
         static unsigned int MaxInt() { return 0xffffff; }  //  2^24 -1

         static int Size();

         
         static std::string Name() {
            return "RanLuxSEngine";
         }

      protected:
         // functions used for testing
         
         void SetState(const std::vector<uint32_t> & state);

         void GetState(std::vector<uint32_t> & state);

         void Init(int level, int seed);


      private:

         double Rndm_impl();
         uint32_t IntRndm_impl() {return (int) Rndm_impl()*MaxInt(); }

         int fLevel;
         float fX[1];

         RanluxS * fRlxs;

      };

      /**
         48 buts engine
       */
     class RanLuxDEngine : public TRandomEngine {


      public:

         typedef  TRandomEngine BaseType;
         typedef  uint32_t Result_t;
         typedef  uint32_t StateInt_t;


        RanLuxDEngine(uint32_t level, uint32_t seed=4357) :
           fRlxd(0)
        {
           fLevel = level;
           Init(level,seed);
         }

         virtual ~RanLuxDEngine();

         void SetSeed(Result_t seed) {
            Init(fLevel,seed);
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

         static int Size();

         
         static std::string Name() {
            return "RanLuxDEngine";
         }

      protected:
         // functions used for testing
         
         void SetState(const std::vector<uint32_t> & state);

         void GetState(std::vector<uint32_t> & state);

         void Init(int level, int seed);


      private:

         double Rndm_impl();
         uint32_t IntRndm_impl() {return (int) Rndm_impl()*MaxInt(); }

         int fLevel;
         double fX[1];

        RanluxD * fRlxd;


      };


   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_TRandomEngines */
