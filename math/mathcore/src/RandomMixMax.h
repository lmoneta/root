// MixMax random number  implementation
//

#ifndef ROOT_Math_RandomMixMax
#define ROOT_Math_RandomMixMax

#include "mixmax.h"
//include "mixmax.c"


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


         void SetSeed(unsigned int seed) { 
            //seed_lcg(fRngState, seed); 
            seedvielbein(fRngState, seed);
         }

         unsigned int GetSeed() const { 
            return get_next(fRngState);
         }
         

         // generate one random number in interval ]0,1]
         double operator() () { 
            return get_next_float(fRngState);
         }


      private: 

         rng_state_t * fRngState;

      };

    }  // end namespace Math
}  // end namespace ROOT

#endif
