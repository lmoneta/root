/*
 * mixmax.h version 0.1
 *
 * C++11 implementation of the MIXMAX random number generator.
 *
 *  Created by Konstantin Savvidy.
 *
 *  The code is released under GNU Lesser General Public License v3
 *
 *	G.K.Savvidy and N.G.Ter-Arutyunian,
 *  On the Monte Carlo simulation of physical systems,
 *	J.Comput.Phys. 97, 566 (1991);
 *  Preprint EPI-865-16-86, Yerevan, Jan. 1986
 *
 *  K.Savvidy
 *  The MIXMAX random number generator
 *  Comp. Phys. Commun. 196 (2015), pp 161–165
 *  http://dx.doi.org/10.1016/j.cpc.2015.06.003
 *
 *  K.Savvidy and G.Savvidy
 *  Spectrum and Entropy of C-systems. MIXMAX random number generator
 *  Chaos, Solitons & Fractals, Volume 91, (2016) pp. 33–38
 *  http://dx.doi.org/10.1016/j.chaos.2016.05.003
 *
 */

#ifndef __MIXMAX_H
#define __MIXMAX_H

#include <array>
#include <vector>
#include <cstdint>

namespace mixmax { 

template <typename T, T __min, T __max> class _Generator
// Boilerplate code, from Andrzej, it is required to be compatible with libstdc++ interfaces, see example.cpp for how to use.
{
public:
   using result_type = T;								
   static constexpr T min() {return __min;}
   static constexpr T max() {return __max;}
   void seed (result_type val = 1);
   T operator()();
};

typedef uint32_t myID_t;
typedef uint64_t myuint;


constexpr int Ndim = 240;

constexpr int BITS=61;
constexpr myuint M61=2305843009213693951ULL;
constexpr myuint MERSBASE=M61;
constexpr double INV_MERSBASE=(0.43368086899420177360298E-18);


/*
  Table of suggested parameters for MIXMAX
 
  Vector size |                                                                                    period q
  N           |    SPECIAL                |   SPECIALMUL   |           MOD_MULSPEC               | log10(q)  |   entropy  | skip
  -------------------------------------------------------------------------------------------------------------------------------
  17          |         0                 |     36         |                none                 |    294    |    374.3   | none
  240         | 487013230256099140        |     51         |   fmodmulM61( 0, SPECIAL , (k) )    |   4389    |   8679.2   | none
  256         |         -1                |     0          |     (MERSBASE - (k))                |   4682    |    164.5   |  2
  256         | 487013230256099064        |     0          |   fmodmulM61( 0, SPECIAL , (k) )    |   4682    |    193.6   | none

  Figure of merit is entropy: best generator overall is N=240
 
*/
// Interface C++11 std::random

namespace ROOT {
   namespace Math {
      template<int N, int S> 
      class MixMaxEngine;
   }
}


template <int Ndim>
class mixmax_engine: public _Generator<std::uint64_t, 0, 0x1FFFFFFFFFFFFFFF> // does not work with any other values
{

   template<int N1,int N2>
   friend class ::ROOT::Math::MixMaxEngine;
   
   static const int N = Ndim;
   static constexpr long long int SPECIAL   = ((N==17)? 0 : ((N==240)? 487013230256099140ULL:-1) ); // etc...
   static constexpr long long int SPECIALMUL= ((N==17)? 36: ((N==240)? 51                   : 0) ); // etc...
   // Note the potential for confusion...
   
   struct rng_state_st
   {
      //myuint V[N];      //
      std::array<myuint, N> V;
      //myuint *V=A.data();
      myuint sumtot;
      int counter;
   };
   
   typedef struct rng_state_st rng_state_t; // C struct alias

   // state of the generator 
   rng_state_t S;

    
public:
   using T = result_type;								  // should it be double?
   static constexpr int rng_get_N() {return N;}
   static constexpr long long int rng_get_SPECIAL()    {return SPECIAL;}
   static constexpr int rng_get_SPECIALMUL() {return SPECIALMUL;}
   void seed_uniquestream( rng_state_t* Xin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
   void print_state();
   void read_state(const char filename[] );
   myuint get_next() ;
   double get_next_float();
   
   int iterate();
   mixmax_engine Branch();
   void BranchDaughter(int b=5); // valid values are between b=5 and b=60
   void BranchMother();
   
   mixmax_engine(myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );	   // Constructor with four 32-bit seeds
   void seed(uint64_t seedval){seed_uniquestream( &S, 0, 0, (myID_t)(seedval>>32), (myID_t)seedval );} // seed with one 64-bit seed
   mixmax_engine(); // Constructor, no seeds
   
   mixmax_engine& operator=(const mixmax_engine& other );
   
   inline T operator()()
   {
      return get_next();
   }
    
private:
   myuint MOD_MULSPEC(myuint k);
   void seed_vielbein(rng_state_t* X, unsigned int i); // seeds with the i-th unit vector, i = 0..N-1,  for testing only
   myuint iterate_raw_vec(myuint* Y, myuint sumtotOld);
   myuint apply_bigskip(myuint* Vout, myuint* Vin, myID_t clusterID, myID_t machineID, myID_t runID, myID_t  streamID );
   myuint modadd(myuint foo, myuint bar);
   myuint fmodmulM61(myuint cum, myuint s, myuint a);
#if defined(__x86_64__)
   inline myuint mod128(__uint128_t s);
#endif
};


#define ARRAY_INDEX_OUT_OF_BOUNDS   0xFF01
#define SEED_WAS_ZERO               0xFF02
#define ERROR_READING_STATE_FILE    0xFF03
#define ERROR_READING_STATE_COUNTER       0xFF04
#define ERROR_READING_STATE_CHECKSUM      0xFF05

}  // end namespace mixmax   

#include "Math/mixmax/mixmax.icc"


#endif		// __MIXMAX_H

//class template mixmax_engine<256>;


