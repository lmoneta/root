//Code generated automatically by TMVA for Inference of Model file [Flattentest.h5] at [Thu Aug 24 08:55:48 202] 

#ifndef ROOT_TMVA_SOFIE_FLATTENTEST
#define ROOT_TMVA_SOFIE_FLATTENTEST

#include<vector>
#include "TMVA/SOFIE_common.hxx"

namespace TMVA_SOFIE_Flattentest{
struct Session {
std::vector<float> fTensor_flatten1Reshape0 = std::vector<float>(4);
float * tensor_flatten1Reshape0 = fTensor_flatten1Reshape0.data();
std::vector<float> fTensor_reshape4Reshape0 = std::vector<float>(4);
float * tensor_reshape4Reshape0 = fTensor_reshape4Reshape0.data();


Session(std::string = "") {
}

std::vector<float> infer(float* tensor_reshape4input){
   ///--------Reshape operator

   std::copy( tensor_reshape4input, tensor_reshape4input + 4, tensor_reshape4Reshape0);
   ///--------Flatten operator

   std::copy( tensor_reshape4Reshape0, tensor_reshape4Reshape0 + 4, tensor_flatten1Reshape0);
   std::vector<float> ret (tensor_flatten1Reshape0, tensor_flatten1Reshape0 + 4);
   return ret;
}
};
} //TMVA_SOFIE_Flattentest

#endif  // ROOT_TMVA_SOFIE_FLATTENTEST
