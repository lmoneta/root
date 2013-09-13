// MixMax random number  implementation
//

#ifndef ROOT_Math_RandomMixMax
#define ROOT_Math_RandomMixMax

#include "random123/include/Random123/threefry.h"
#include "random123/include/Random123/u01.h"


namespace ROOT { 
   namespace Math { 


// this is needed for CINT ? 
// #ifdef __CINT
//       struct r123array2x32 {
//          UInt_t v[2];
//          UInt_t & operator[](UInt_t i){return v[i];}
//       };


      
//       typedef r123array2x32 threefry2x32_ctr_t;
//       typedef r123array2x32 threefry2x32_key_t;
//       typedef r123array2x32 threefry2x32_ukey_t;
// #endif
      
      class Random123 { 

      public:

         // generate one random number in interval ]0,1]
         double operator() () { 
               fCtr.v[0]++;
               fCtr = threefry2x32(fCtr, fKey);
               return fCtr.v[0] * R123_0x1p_32f;
         }

         unsigned int GetSeed() const { 
            return  fKey.v[0];
         }

         void SetSeed(unsigned int seed) { 
            fKey = (r123array2x32){{seed, 0xdecafbad}};
            fCtr = (r123array2x32){{0, 0xf00dcafe}};
         }

      private: 
         
         threefry2x32_key_t fKey;
         threefry2x32_ctr_t fCtr;

      };

    }  // end namespace Math
}  // end namespace ROOT

#endif
