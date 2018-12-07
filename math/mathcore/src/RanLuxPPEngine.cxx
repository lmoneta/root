// @(#)root/mathcore:$Id$
// Authors: L. Moneta    8/2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015 , ROOT MathLib Team                             *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// implementation file of Mersenne-Twister engine
//
//
// Created by: Lorenzo Moneta  : Tue 4 Aug 201
//
//
#include "Math/RanLuxPPEngine.h"
#include "ranluxpp.h"
#include <cassert>



namespace ROOT {
namespace Math {

   ///////RanluxPP  engine

  uint64_t skipValues[] = {24,48,97,223,389,1024,2048};

  RanLuxPPEngine::RanLuxPPEngine(uint64_t seed, int luxlevel) {
    uint64_t p = 2048;
    if (luxlevel >=0 &&  luxlevel < 7) p = skipValues[luxlevel];    
    fRlxpp = new ranluxpp(seed, p); 
   }


  RanLuxPPEngine::~RanLuxPPEngine() {
      if (fRlxpp) delete fRlxpp; 
   }

  /// init the generator
  void RanLuxPPEngine::Init(uint64_t seed) {
    if (fRlxpp) return;
    fRlxpp->init(seed);
      //rlxd_init(level,seed);
   }

   /// generate a random double number 
   double RanLuxPPEngine::Rndm_impl() {
      // generate a double random number
     return (*fRlxpp)(static_cast<double>(1.)); 
   }
   /// generate a double random double number 
   float RanLuxPPEngine::Rndm_impl_float() {
      // could use a counter to be faster
     return (*fRlxpp)(static_cast<float>(1.)); 
   }

   void RanLuxPPEngine::SetState(const std::vector<uint64_t> & state) {
      assert( (int) state.size() == Size() );
      uint64_t * s = fRlxpp->getstate();
      std::copy(state.begin(), state.end(), s); 
   }

   void RanLuxPPEngine::GetState(std::vector<uint64_t> & state) {
      size_t size = Size();
      uint64_t * s = fRlxpp->getstate( );
      state = std::vector<uint64_t>( s, s+size); 
   }

   // int RanLuxPPEngine::Size() {
   //    return  9; 
   // }
      
  
   } // namespace Math
} // namespace ROOT
