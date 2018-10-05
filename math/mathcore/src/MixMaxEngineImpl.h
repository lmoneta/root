#ifndef ROOT_Math_MixMaxEngineImpl
#define ROOT_Math_MixMaxEngineImpl


#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <vector>



#if (ROOT_MM_N==17)
namespace mixmax_17 {
#elif (ROOT_MM_N==240)
namespace mixmax_240 {
#elif (ROOT_MM_N==256)
#define USE_MIXMAX_OLD   
namespace mixmax_256 {
#elif (ROOT_MM_N==8)
namespace mixmax_8 {
#elif (ROOT_MM_N==40)
namespace mixmax_40 { 
#elif (ROOT_MM_N==44)
#define USE_MIXMAX_OLD   
namespace mixmax_44 { 
#elif (ROOT_MM_N==88)
#define USE_MIXMAX_OLD   
namespace mixmax_88 { 
#elif (ROOT_MM_N==10)
#define USE_MIXMAX_OLD   
namespace mixmax_10 { 
#elif
error "value of ROOT_MM_N is not supported"
#endif

#ifdef WIN32
#define __thread __declspec(thread)
#endif

#include "mixmax.icc"

#undef N
}  // end namespace mixmax_xx

#include "Math/MixMaxEngine.h"

#include <iostream>

#if (ROOT_MM_N==17)
using namespace mixmax_17;
#elif (ROOT_MM_N==240)
using namespace mixmax_240;
#elif (ROOT_MM_N==256)
using namespace mixmax_256;
#elif (ROOT_MM_N==8)
using namespace mixmax_8;
#elif (ROOT_MM_N==40)
using namespace mixmax_40;
#elif (ROOT_MM_N==44)
using namespace mixmax_44;
#elif (ROOT_MM_N==88)
using namespace mixmax_88;
#elif (ROOT_MM_N==10)
using namespace mixmax_10;
#endif


namespace ROOT {
   namespace Math {


         // dummy implementation
   template<int N>
   class MixMaxEngineImpl {
   public:
      MixMaxEngineImpl(uint64_t) {
         std::cerr << "MixMaxEngineImpl - These template parameters are not supported for MixMaxEngine" << std::endl;
      }
      ~MixMaxEngineImpl() {}
      void SetSeed(uint64_t) { }
      double Rndm() { return -1; }
      double IntRndm() { return 0; }
      void SetState(const std::vector<uint64_t> &) { }
      void GetState(std::vector<uint64_t> &) { }
      int Counter() { return -1; }
      void SetCounter(int) {}
      void Iterate() {}
   };



template<> 
class MixMaxEngineImpl<ROOT_MM_N> {
   rng_state_t * fRngState;
public:

   typedef MixMaxEngine<ROOT_MM_N,0>::StateInt_t StateInt_t; 
   typedef MixMaxEngine<ROOT_MM_N,0>::Result_t Result_t;
   
   MixMaxEngineImpl(uint64_t seed) {
      fRngState = rng_alloc();
      SetSeed(seed); 
   }
   ~MixMaxEngineImpl() {
      rng_free(fRngState);
   }
   void SetSeedFast(Result_t seed) {
      seed_spbox(fRngState, seed);
   }
   void SetSeed(Result_t seed) {
#if (ROOT_MM_N==10 || ROOT_MM_N==40 || ROOT_MM_N==44 || ROOT_MM_N==88)
      SetSeedFast(seed);
#else
      seed_uniquestream(fRngState, 0, 0, (uint32_t)(seed>>32), (uint32_t)seed );
#endif
   }
   double Rndm() {
       return get_next_float(fRngState);
   }
   // generate one integer number 
   Result_t IntRndm() {
      return get_next(fRngState);
   }
   void SetState(const std::vector<StateInt_t> & state) {
      if (fRngState) rng_free(fRngState);
      fRngState = rng_copy(const_cast<StateInt_t*>(state.data()) );
   }
   void GetState(std::vector<StateInt_t> & state) const {
      int n =  rng_get_N(); 
      state.resize(n);
      for (int i = 0; i < n; ++i)
         state[i] = fRngState->V[i];
   }
   void Iterate() {
      iterate(fRngState); 
   }
   int Counter() const {
      return fRngState->counter; 
   }
   void SetCounter(int val) {
      fRngState->counter = val; 
   }
   static int Size()  {
      return rng_get_N(); 
   }

   // to silent some warning
   void RndmArray(int n, double * array) {
      fill_array(fRngState, n, array); 
   }
   void ReadState(const char filename[] ) {
      read_state(fRngState, filename);
   }
   // branch generator given a vector of seed (at least 4 32 bit values)
   void Branch(uint32_t * seedvec) {
      branch_inplace(fRngState, seedvec); 
   }
     
   
};

      
   } // end namesapce Math
} // end namespace ROOT

#endif
