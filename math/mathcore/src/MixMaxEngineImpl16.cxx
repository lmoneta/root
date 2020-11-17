

// define number for used to Mixmax
#define ROOT_MM_N 16
#include "MixMaxEngineImpl.h"

// define the template instance we want to have in the librsary
//( need to be declared as extern template in the .h file)
namespace ROOT {
   namespace Math {
      template class MixMaxEngine<ROOT_MM_N,0>;
      template class MixMaxEngine<ROOT_MM_N,1>;
      template class MixMaxEngine<ROOT_MM_N,2>;
      template class MixMaxEngine<ROOT_MM_N,5>;
      template class MixMaxEngine<ROOT_MM_N,11>;
   }
}

