//Code generated automatically by TMVA for Inference of Model file [MLPtest.h5] at [Thu Aug 24 08:55:44 202] 

#ifndef ROOT_TMVA_SOFIE_MLPTEST
#define ROOT_TMVA_SOFIE_MLPTEST

#include<algorithm>
#include<cmath>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_MLPtest{
namespace BLAS{
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
}//BLAS
struct Session {
std::vector<float> fTensor_dense4bias0 = std::vector<float>(2);
float * tensor_dense4bias0 = fTensor_dense4bias0.data();
std::vector<float> fTensor_dense3bias0 = std::vector<float>(64);
float * tensor_dense3bias0 = fTensor_dense3bias0.data();
std::vector<float> fTensor_dense4kernel0 = std::vector<float>(128);
float * tensor_dense4kernel0 = fTensor_dense4kernel0.data();
std::vector<float> fTensor_dense2bias0 = std::vector<float>(64);
float * tensor_dense2bias0 = fTensor_dense2bias0.data();
std::vector<float> fTensor_dense1bias0 = std::vector<float>(64);
float * tensor_dense1bias0 = fTensor_dense1bias0.data();
std::vector<float> fTensor_dense1kernel0 = std::vector<float>(4096);
float * tensor_dense1kernel0 = fTensor_dense1kernel0.data();
std::vector<float> fTensor_densebias0 = std::vector<float>(64);
float * tensor_densebias0 = fTensor_densebias0.data();
std::vector<float> fTensor_dense3kernel0 = std::vector<float>(4096);
float * tensor_dense3kernel0 = fTensor_dense3kernel0.data();
std::vector<float> fTensor_dense2kernel0 = std::vector<float>(4096);
float * tensor_dense2kernel0 = fTensor_dense2kernel0.data();
std::vector<float> fTensor_densekernel0 = std::vector<float>(448);
float * tensor_densekernel0 = fTensor_densekernel0.data();
std::vector<float> fTensor_dense4Sigmoid0 = std::vector<float>(2);
float * tensor_dense4Sigmoid0 = fTensor_dense4Sigmoid0.data();
std::vector<float> fTensor_dense4bias0bcast = std::vector<float>(2);
float * tensor_dense4bias0bcast = fTensor_dense4bias0bcast.data();
std::vector<float> fTensor_dense3Relu0 = std::vector<float>(64);
float * tensor_dense3Relu0 = fTensor_dense3Relu0.data();
std::vector<float> fTensor_dense2Dense = std::vector<float>(64);
float * tensor_dense2Dense = fTensor_dense2Dense.data();
std::vector<float> fTensor_dense4Dense = std::vector<float>(2);
float * tensor_dense4Dense = fTensor_dense4Dense.data();
std::vector<float> fTensor_denseSelu0 = std::vector<float>(64);
float * tensor_denseSelu0 = fTensor_denseSelu0.data();
std::vector<float> fTensor_dense1Dense = std::vector<float>(64);
float * tensor_dense1Dense = fTensor_dense1Dense.data();
std::vector<float> fTensor_dense3Dense = std::vector<float>(64);
float * tensor_dense3Dense = fTensor_dense3Dense.data();
std::vector<float> fTensor_denseDense = std::vector<float>(64);
float * tensor_denseDense = fTensor_denseDense.data();
std::vector<float> fTensor_dense3bias0bcast = std::vector<float>(64);
float * tensor_dense3bias0bcast = fTensor_dense3bias0bcast.data();
std::vector<float> fTensor_dense2Sigmoid0 = std::vector<float>(64);
float * tensor_dense2Sigmoid0 = fTensor_dense2Sigmoid0.data();
std::vector<float> fTensor_dense1bias0bcast = std::vector<float>(64);
float * tensor_dense1bias0bcast = fTensor_dense1bias0bcast.data();
std::vector<float> fTensor_densebias0bcast = std::vector<float>(64);
float * tensor_densebias0bcast = fTensor_densebias0bcast.data();
std::vector<float> fTensor_dense2bias0bcast = std::vector<float>(64);
float * tensor_dense2bias0bcast = fTensor_dense2bias0bcast.data();
std::vector<float> fTensor_dense1Tanh0 = std::vector<float>(64);
float * tensor_dense1Tanh0 = fTensor_dense1Tanh0.data();


Session(std::string filename ="") {
   if (filename.empty()) filename = "MLPtest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense4bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense4bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense4bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense3bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense3bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense3bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense4kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense4kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 128) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 128 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense4kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense2bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense2bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense2bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense1bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense1bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense1bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense1kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense1kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense1kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_densebias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_densebias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_densebias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense3kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense3kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense3kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense2kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense2kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense2kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_densekernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_densekernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 448) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 448 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_densekernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_densebias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_densebias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense1bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense1bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense2bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense2bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense3bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense3bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense4bias0,{ 2 }, { 1 , 2 });
      std::copy(data, data + 2, tensor_dense4bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_denseinput){

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
   std::copy(tensor_densebias0bcast, tensor_densebias0bcast + 64, tensor_denseDense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_densekernel0, &op_0_ldb, tensor_denseinput, &op_0_lda, &op_0_beta, tensor_denseDense, &op_0_n);
	for (int id = 0; id < 64 ; id++){
		tensor_denseSelu0[id] = 1.0507009873554804934193349852946 * (std::max(float(0.0), tensor_denseDense[id]) + std::min(0.0, 1.6732632423543772848170429916717 * (std::exp(tensor_denseDense[id])-1)));
	}

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 1;
   int op_2_n = 64;
   int op_2_k = 64;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 64;
   int op_2_ldb = 64;
   std::copy(tensor_dense1bias0bcast, tensor_dense1bias0bcast + 64, tensor_dense1Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense1kernel0, &op_2_ldb, tensor_denseSelu0, &op_2_lda, &op_2_beta, tensor_dense1Dense, &op_2_n);

//------ TANH
   for (int id = 0; id < 64 ; id++){
      tensor_dense1Tanh0[id] = std::tanh(tensor_dense1Dense[id]);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 1;
   int op_4_n = 64;
   int op_4_k = 64;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 64;
   int op_4_ldb = 64;
   std::copy(tensor_dense2bias0bcast, tensor_dense2bias0bcast + 64, tensor_dense2Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense2kernel0, &op_4_ldb, tensor_dense1Tanh0, &op_4_lda, &op_4_beta, tensor_dense2Dense, &op_4_n);
	for (int id = 0; id < 64 ; id++){
		tensor_dense2Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense2Dense[id]));
	}

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 1;
   int op_6_n = 64;
   int op_6_k = 64;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 64;
   int op_6_ldb = 64;
   std::copy(tensor_dense3bias0bcast, tensor_dense3bias0bcast + 64, tensor_dense3Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense3kernel0, &op_6_ldb, tensor_dense2Sigmoid0, &op_6_lda, &op_6_beta, tensor_dense3Dense, &op_6_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense3Relu0[id] = ((tensor_dense3Dense[id] > 0 )? tensor_dense3Dense[id] : 0);
   }

//--------- Gemm
   char op_8_transA = 'n';
   char op_8_transB = 'n';
   int op_8_m = 1;
   int op_8_n = 2;
   int op_8_k = 64;
   float op_8_alpha = 1;
   float op_8_beta = 1;
   int op_8_lda = 64;
   int op_8_ldb = 2;
   std::copy(tensor_dense4bias0bcast, tensor_dense4bias0bcast + 2, tensor_dense4Dense);
   BLAS::sgemm_(&op_8_transB, &op_8_transA, &op_8_n, &op_8_m, &op_8_k, &op_8_alpha, tensor_dense4kernel0, &op_8_ldb, tensor_dense3Relu0, &op_8_lda, &op_8_beta, tensor_dense4Dense, &op_8_n);
	for (int id = 0; id < 2 ; id++){
		tensor_dense4Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense4Dense[id]));
	}
   std::vector<float> ret (tensor_dense4Sigmoid0, tensor_dense4Sigmoid0 + 2);
   return ret;
}
};
} //TMVA_SOFIE_MLPtest

#endif  // ROOT_TMVA_SOFIE_MLPTEST
