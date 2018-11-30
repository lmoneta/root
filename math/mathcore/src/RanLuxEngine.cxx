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
//#include "ranlxs.h"
#include "ranlxd.h"
#include <cassert>

#include "RanluxS.h"
#include "RanluxD.h"


namespace ROOT {
namespace Math {

   RanLuxSEngine::~RanLuxSEngine() {
      if (fRlxs) delete fRlxs; 
   }

   /// init the generator
   void RanLuxSEngine::Init(int level, int seed) {
//      rlxs_init(level,seed);

      if (fRlxs) delete fRlxs;
      fRlxs = new RanluxS(); 
      fRlxs->rlxs_init(level,seed);
   }

   /// generate a random double number 
   double RanLuxSEngine::Rndm_impl() {
      // could use a counter to be faster
      //ranlxs(fX,1);
      fRlxs->ranlxs(fX,1);
      return fX[0]; 
   }

   void RanLuxSEngine::SetState(const std::vector<uint32_t> & state) {
      fRlxs->rlxs_reset((int*) state.data() );
   }

   // void RanLuxSEngine::GetState(std::vector<uint32_t> & state) {
   //    fRlxs->rlxs_reset((int*) state.data() );
   // }

   void RanLuxSEngine::GetState(std::vector<uint32_t> & state) {
      fRlxs->rlxs_get((int*) state.data() );
   }

   int RanLuxSEngine::Size() {
      return RanluxS::rlxs_size(); 
   }

   ///////Double engine

    RanLuxDEngine::~RanLuxDEngine() {
      if (fRlxd) delete fRlxd; 
   }

     /// init the generator
   void RanLuxDEngine::Init(int level, int seed) {
      if (fRlxd) delete fRlxd;
      fRlxd = new RanluxD(); 
      fRlxd->rlxd_init(level,seed);
      //rlxd_init(level,seed);
   }

   /// generate a random double number 
   double RanLuxDEngine::Rndm_impl() {
      // could use a counter to be faster
      fRlxd->ranlxd(fX,1);
      return fX[0]; 
   }

   void RanLuxDEngine::SetState(const std::vector<uint32_t> & state) {
      assert( (int) state.size() == Size() );
      fRlxd->rlxd_reset((int*) state.data() );
   }

   void RanLuxDEngine::GetState(std::vector<uint32_t> & state) {
      size_t size = Size(); 
      if (state.size() != size) state.resize(size); 
      fRlxd->rlxd_get((int*) state.data() );
   }

   int RanLuxDEngine::Size() {
      return  RanluxD::rlxd_size(); 
   }
      
  
   } // namespace Math
} // namespace ROOT
