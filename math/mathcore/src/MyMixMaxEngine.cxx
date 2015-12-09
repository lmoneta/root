// @(#)root/mathcore:$Id$
// Authors: L. Moneta    8/2015

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2015 , ROOT MathLib Team                             *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// implementation file of MixMax engine
//
//
// Created by: Lorenzo Moneta  : Tue 4 Aug 2015
//
//
#include "Math/MyMixMaxEngine.h"
#include "Math/MersenneTwisterEngine.h"

#include <iostream>



int state_size = 256;
int special = 0;
const uint64_t M61 = 2305843009213693951ULL; // 2^61-1

bool debug = false; 


namespace ROOT {
namespace Math {

   void MyMixMaxEngine::SetSize(int n)  {
      state_size = n;
   }
   void MyMixMaxEngine::SetSpecial(int s)  {
      special = s;
   }

   MyMixMaxEngine::MyMixMaxEngine(uint64_t seed, int n) { 
   // fill the state
      if (n <= 3) n = state_size;
      else    state_size = n;
   fState.resize(n);
   fMatrix.resize(n * n);
   fCounter = state_size;
      MersenneTwisterEngine rndm(seed);
      //const int imax32 = std::numeric_limits<uint32_t>::max(); 
      for (int i = 0; i < n; ++i) {
         uint64_t s  =  (uint64_t) rndm.IntRndm() << 32 | rndm.IntRndm();
         fState[i] = (s & M61) + (s >> 61);  // to make a 61 bit integer
      }
      // fill the matrix 
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
            fMatrix[i*n  + j ] = 1;
         }
         for (int j = 1; j <= i ; ++j) {
            fMatrix[i*n + j ] += i - (j-1);
         }
      }
      // add special
      fMatrix[2 *n + 1] += special;

      // print
      if (debug) { 
      std::cout << "Matrix for n = " << state_size << std::endl;
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j) {
            std::cout << fMatrix[i*n+j] << "    "; 
         }
         std::cout << std::endl;
      }
      std::cout << std::endl;
   
      std::cout << "initial state " << std::endl;
      for   (int i = 0; i < n; ++i)
         std   ::cout << fState[i] << "   ";
      std::cout << std::endl;
      }
}



   // void MixMaxEngine::SeedUniqueStream(unsigned int clusterID, unsigned int machineID, unsigned int runID, unsigned int  streamID) { 
   //    seed_uniquestream(fRngState, clusterID,  machineID,  runID,   streamID);
   // }


   int MyMixMaxEngine::Size()  {
      return state_size; 
   }


   // unsigned int MyMixMaxEngine::GetSeed() const { 
   //    return get_next(fRngState);
   // }
         

   // // generate one random number in interval ]0,1]
   // double MyMixMaxEngine::Rndm_impl()  { 
      
   // }

   // generate one integer number 
   uint64_t MyMixMaxEngine::IntRndm() { 
      // perform the matrix multiplications
      if (fCounter >= state_size) {
         // iterate
         int n = state_size; 
         std::vector<uint64_t> newstate(n);

#ifdef MATRIX_IMPL         
         for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
               newstate[i] += fMatrix[i*n+j]*fState[j];
            }
            newstate[i] = (newstate[i] & M61) + (newstate[i] >> 61);  // to make a 61 bit integer
         }
#else
         // efficient impelmentation as described in the paper
         newstate[0] = fSumTot;
         uint64_t b_i = 0;
         uint64_t sumtot = newstate[0];
         // loop starting from 1 
         for (int i = 1; i < n; ++i) {
            b_i += fState[i];
            newstate[i] = newstate[i-1] + b_i;
            sumtot += newstate[i]; 
         }
         // special correction
         newstate[2] += special * fState[1];
         sumtot += special * fState[1];
         fSumTot = sumtot; 
#endif         
         fState = newstate;
         fCounter = 0; 
      }
      return fState[fCounter++];
   }

                  
   // void  MyMixMaxEngine::RndmArray(int n, double *array){
   //    // Return an array of n random numbers uniformly distributed in ]0,1]
   //    fill_array(fRngState, n,  array);
   // }



   
  
   } // namespace Math
} // namespace ROOT
