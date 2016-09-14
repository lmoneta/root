/*
 * mixmax.cxx version 0.1
 *
 *
 */

//#include "mixmax.h"
typedef unsigned long long myuint;

namespace mixmax {

   void get_skip_matrix_17(const myuint** skipMat) {

      const	myuint skipMat17[128][17] =
#include "mixmax_skip_N17.icc"
         ;

      for (int i=0; i<128; i++) { skipMat[i] = skipMat17[i];}
   }

   void get_skip_matrix_240(const myuint** skipMat) {
      const	myuint skipMat240[128][240] =
#include "mixmax_skip_N240.icc"
         ;

      for (int i=0; i<128; i++) { skipMat[i] = skipMat240[i];}
   }

   void get_skip_matrix_256_oldS(const myuint**  skipMat) {
      const	myuint skipMat256old[128][256] =
#include "mixmax_skip_N256.oldS.icc"
         ;

      for (int i=0; i<128; i++) { skipMat[i] = skipMat256old[i];}

   }
}
