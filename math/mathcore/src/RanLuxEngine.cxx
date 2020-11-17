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

#include <iostream>

namespace ROOT {
namespace Math {

   RanLuxSEngine::~RanLuxSEngine() {
      if (fRlxs) delete fRlxs; 
   }

   /// init the generator
   void RanLuxSEngine::Init(int level, int seed) {
//      rlxs_init(level,seed);

      uint32_t maxval = 2147483648; // 2^31
      seed = seed % maxval;
      std::cout << "init with seed " << seed << std::endl;
      if (seed == 0) seed =1 ; 
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
      // state contains extra bits (e.g. carry bits , size etc..)
      std::vector<int> actState( RanluxS::rlxs_size() );
      fRlxs->rlxs_get(actState.data() );
      // copy state in elements from 1 to 96
      for (int i = 0; i < 96; ++i)
         actState[i+1] = state[i] % (1 << 24);

      for (int i = 97; i < 100; ++i) actState[i] = 0; 
      fRlxs->rlxs_reset(actState.data() );
   }

   // void RanLuxSEngine::GetState(std::vector<uint32_t> & state) {
   //    fRlxs->rlxs_reset((int*) state.data() );
   // }

   void RanLuxSEngine::GetState(std::vector<uint32_t> & state) {
      std::vector<int> actState( RanluxS::rlxs_size() );
      fRlxs->rlxs_get(actState.data() );
      std::copy(actState.begin()+1, actState.begin()+97, state.begin() ); 
   }

   int RanLuxSEngine::Size() {
      return 96; 
   }

   ///////Double engine

    RanLuxDEngine::~RanLuxDEngine() {
      if (fRlxd) delete fRlxd; 
   }

     /// init the generator
   void RanLuxDEngine::Init(int level, int seed) {
      uint32_t maxval = 2147483648; // 2^31
      seed = seed % maxval;
      //std::cout << "init with seed " << seed << std::endl;
      if (seed == 0) seed =1 ; 

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
      // state contains extra bits (e.g. carry bits , size etc..)
      std::vector<int> actState( RanluxD::rlxd_size() );
       fRlxd->rlxd_get(actState.data() );
      // copy state in elements from 1 to 96
      for (int i = 0; i < 96; ++i)
         actState[i+1] = state[i] % (1 << 24);
      
      //assert( (int) state.size() == Size() );
      fRlxd->rlxd_reset( actState.data() );
   }

   void RanLuxDEngine::GetState(std::vector<uint32_t> & state) {
      std::vector<int> actState( RanluxD::rlxd_size() );
      fRlxd->rlxd_get(actState.data() );
      std::copy(actState.begin()+1, actState.begin()+97, state.begin() ); 
      // size_t size = Size(); 
      // if (state.size() != size) state.resize(size); 
      // fRlxd->rlxd_get((int*) state.data() );
   }

   int RanLuxDEngine::Size() {
      return  96; //RanluxD::rlxd_size(); 
   }
      
  
   } // namespace Math
} // namespace ROOT
