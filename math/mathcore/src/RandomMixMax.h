/*  MixMax random number  implementation
*	Generator described in 
*	N.Z.Akopov, G.K.Savvidy and N.G.Ter-Arutyunian, Matrix Generator of Pseudorandom Numbers, 
*	J.Comput.Phys. 97, 573 (1991); 
*	Preprint EPI-867(18)-86, Yerevan Jun.1986;
*/

#ifndef ROOT_Math_RandomMixMax
#define ROOT_Math_RandomMixMax

#include "mixmax.h"


namespace ROOT { 
   namespace Math { 


      
      class RandomMixMax { 

      public:

         RandomMixMax() { 
            fRngState = rng_alloc();
         }

         ~RandomMixMax() { 
            rng_free(fRngState);
         }


		  void SeedUniqueStream(UInt_t clusterID, UInt_t machineID, UInt_t runID, UInt_t  streamID) { 
			  seed_uniquestream(fRngState, clusterID,  machineID,  runID,   streamID);
		  }

		  void SetSeed(UInt_t seed) { 
			  seed_spbox(fRngState, seed); iterate(fRngState);			  
            //seed_lcg(fRngState, seed); 
            //seed_vielbein(fRngState, seed);
         }

		  void SetSeed64(ULong64_t seed) { 
			  seed_spbox(fRngState, seed); iterate(fRngState);			  
		 }

		 unsigned int GetSeed() const { 
            return get_next(fRngState);
         }
         

         // generate one random number in interval ]0,1]
         double operator() () { 
            return get_next_float(fRngState);
         }

		  
		virtual  void  RndmArray(Int_t n, Double_t *array){
			 // Return an array of n random numbers uniformly distributed in ]0,1]
			  fill_array(fRngState, n,  array);
		  }

      private: 

         rng_state_t * fRngState;

      };

    }  // end namespace Math
}  // end namespace ROOT

#endif
