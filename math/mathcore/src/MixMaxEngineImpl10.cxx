

// define number for used to Mixmax
#define ROOT_MM_N 10
#include "MixMaxEngineImpl.h"

// define the template instance we want to have in the librsary
//( need to be declared as extern template in the .h file)
namespace ROOT {
   namespace Math {
      template class MixMaxEngine<ROOT_MM_N,0>;
      template class MixMaxEngine<ROOT_MM_N,4>;
      template class MixMaxEngine<ROOT_MM_N,7>;
      template class MixMaxEngine<ROOT_MM_N,13>;
      template class MixMaxEngine<ROOT_MM_N,14>;
   }
}

