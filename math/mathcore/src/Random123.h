// MixMax random number  implementation
//

#ifndef ROOT_Math_Random123
#define ROOT_Math_Random123

#include "random123/include/Random123/threefry.h"
#include "random123/include/Random123/philox.h"
#include "random123/include/Random123/u01.h"


namespace ROOT { 
   namespace Math { 



      class Random123 { 

      public:

         // default for Threefry is 20 for Philox is 10
         //typedef r123::Threefry4x32_R<THREEFRY2x32_DEFAULT_ROUNDS> RNG;
         typedef r123::Philox4x64_R<7> RNG;
         typedef RNG::ctr_type::value_type result_type; 
         typedef RNG::ctr_type ctr_type; 
         typedef RNG::key_type key_type; 
         typedef RNG::ukey_type ukey_type; 
         // typedef r123array4x64 ARRAY;
         // typedef r123array2x64 ARRAY2;

         static const int BITS = 32;


         Random123(int seed = 1): fNcount(0), fLastElem(0) {           
            for (unsigned int i = 0; i < fCtr.size(); ++i) fCtr.v[i] = 0;       
            for (unsigned int i = 0; i < fKey.size(); ++i) fKey.v[i] = 0; 
            fKey.v[0] = seed;

            const size_t W = std::numeric_limits<result_type>::digits;
            assert (W >= BITS); // should be a static assert (t.b.d. in C++11)
            
         }

         // ~Random123() { 
         //    // printf("vmax value = %f  - %f\n ",(double) fMax, double(fMax)* R123_0x1p_64); 
         // }


         inline void FillBuffer() { 
            // jam n into the high bits of c
            // use same procedure as in random123/MicroURNG.hpp
            // could move to use that when C++11 is available
            const size_t W = std::numeric_limits<result_type>::digits;
            size_t N = fCtr.size(); 
            ctr_type c = fCtr; 
            c.v[N-1] |= fNcount<< ( W - BITS);
            fBuffer = generate(c, fKey);
            fNcount++;
            fLastElem = fBuffer.size(); 
         }

         // generate one random number in interval ]0,1]
         inline double operator() () { 
            if (fLastElem == 0 ) {
               FillBuffer();
               //printf("Size %d \n",(int) fCtr.size() );
            }
            //return (fBuffer.v[fBufCount++] >> 11) * R123_0x1p_53;
            //  return (fBuffer.v[fBufCount++] >> 11) * R123_0x1p_53;
          
            return double(fBuffer.v[--fLastElem]) * R123_0x1p_64;

         }

         unsigned int GetSeed() const { 
            return  fKey.v[0];
         }

         void SetSeed(unsigned int seed) { 
            fKey = (RNG::ukey_type){{seed, 0xdecafbad}};
            fCtr = (RNG::ctr_type){{0, 0xf00dcafe}};
            fLastElem = 0;
         }

         // method to generate an array of numbers
         template<class Iterator>  
         void GenerateArray(Iterator begin, Iterator end) { 
            unsigned int n = end-begin; 
            unsigned int mlast = n%fBuffer.size(); // remaining numbers to generate
            Iterator lastPos = end-mlast; 
            for (Iterator it = begin; it < lastPos; it += fBuffer.size() ) {
               FillBuffer(); 
               for (unsigned int i = 0; i < fBuffer.size() ; ++i)
                  *(it+i) = (fBuffer.v[i] >> 11) * R123_0x1p_53;
            }
            if (mlast>0) { 
               FillBuffer(); 
               for (unsigned int i = 0; i < mlast; ++i)
                  *(lastPos+i) = (fBuffer.v[i] >> 11) * R123_0x1p_53;
            }
            
         } 
         

      private: 

         ctr_type fCtr;
         ctr_type fBuffer;
         ukey_type fKey;
         R123_ULONG_LONG fNcount;
         size_t fLastElem;
         RNG  generate;

         // threefry2x32_key_t fKey;
         // threefry2x32_ctr_t fCtr;

         RNG::ctr_type::value_type fMax;
      };

    }  // end namespace Math
}  // end namespace ROOT

#endif
