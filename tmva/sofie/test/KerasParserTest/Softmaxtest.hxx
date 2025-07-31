//Code generated automatically by TMVA for Inference of Model file [Softmaxtest.h5] at [Thu Aug 24 08:55:43 202] 

#ifndef ROOT_TMVA_SOFIE_SOFTMAXTEST
#define ROOT_TMVA_SOFIE_SOFTMAXTEST

#include<algorithm>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_Softmaxtest{
namespace BLAS{
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
}//BLAS
struct Session {
std::vector<float> fTensor_dense13bias0 = std::vector<float>(64);
float * tensor_dense13bias0 = fTensor_dense13bias0.data();
std::vector<float> fTensor_dense13kernel0 = std::vector<float>(448);
float * tensor_dense13kernel0 = fTensor_dense13kernel0.data();
std::vector<float> fTensor_dense13Softmax0 = std::vector<float>(64);
float * tensor_dense13Softmax0 = fTensor_dense13Softmax0.data();
std::vector<float> fTensor_dense13Dense = std::vector<float>(64);
float * tensor_dense13Dense = fTensor_dense13Dense.data();
std::vector<float> fTensor_dense13bias0bcast = std::vector<float>(64);
float * tensor_dense13bias0bcast = fTensor_dense13bias0bcast.data();


Session(std::string filename ="") {
   if (filename.empty()) filename = "Softmaxtest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense13bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense13bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense13bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense13kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense13kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 448) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 448 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense13kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense13bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense13bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_dense13input){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 64;
   int op_0_k = 7;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 7;
   int op_0_ldb = 64;
   std::copy(tensor_dense13bias0bcast, tensor_dense13bias0bcast + 64, tensor_dense13Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense13kernel0, &op_0_ldb, tensor_dense13input, &op_0_lda, &op_0_beta, tensor_dense13Dense, &op_0_n);

   //------ SOFTMAX
   for (size_t n = 0; n < 1 ; n++){
         float sum = 0.;
         size_t index = 0+ n * 64;
         for (size_t i = 0; i < 64; i++) {
            tensor_dense13Softmax0[index + i*1] = std::exp(tensor_dense13Dense[index + i*1]);
            sum += tensor_dense13Softmax0[index + i*1];
         }
         for (size_t i = 0; i < 64; i++) {
            tensor_dense13Softmax0[index + i*1] /= sum;
         }
   }
   std::vector<float> ret (tensor_dense13Softmax0, tensor_dense13Softmax0 + 64);
   return ret;
}
};
} //TMVA_SOFIE_Softmaxtest

#endif  // ROOT_TMVA_SOFIE_SOFTMAXTEST
