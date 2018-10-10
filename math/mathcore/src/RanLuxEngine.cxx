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
#include "Math/RanLuxEngine.h"
#include "ranlxs.h"
#include "ranlxd.h"
#include <cassert>


namespace ROOT {
namespace Math {

   /// init the generator
   void RanLuxSEngine::Init(int level, int seed) {
      rlxs_init(level,seed);
   }

   /// generate a random double number 
   double RanLuxSEngine::Rndm_impl() {
      // could use a counter to be faster
      ranlxs(fX,1);
      return fX[0]; 
   }

   void RanLuxSEngine::SetState(const std::vector<uint32_t> & state) {
      rlxs_reset((int*) state.data() );
   }

   void RanLuxSEngine::GetState(std::vector<uint32_t> & state) {
      rlxs_get((int*) state.data() );
   }

   int RanLuxSEngine::Size() {
      return rlxs_size(); 
   }

   ///////Double engine
      
     /// init the generator
   void RanLuxDEngine::Init(int level, int seed) {
      rlxd_init(level,seed);
   }

   /// generate a random double number 
   double RanLuxDEngine::Rndm_impl() {
      // could use a counter to be faster
      ranlxd(fX,1);
      return fX[0]; 
   }

   void RanLuxDEngine::SetState(const std::vector<uint32_t> & state) {
      assert( (int) state.size() == Size() );
      rlxd_reset((int*) state.data() );
   }

   void RanLuxDEngine::GetState(std::vector<uint32_t> & state) {
      size_t size = Size(); 
      if (state.size() != size) state.resize(size); 
      rlxd_get((int*) state.data() );
   }

   int RanLuxDEngine::Size() {
      return rlxd_size(); 
   }
      
  
   } // namespace Math
} // namespace ROOT
