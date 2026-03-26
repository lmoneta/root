//Code generated automatically by TMVA for Inference of Model file [BatchNormalizationtest.h5] at [Thu Aug 24 08:55:44 202] 

#ifndef ROOT_TMVA_SOFIE_BATCHNORMALIZATIONTEST
#define ROOT_TMVA_SOFIE_BATCHNORMALIZATIONTEST

#include<algorithm>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_BatchNormalizationtest{
namespace BLAS{
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void scopy_(const int *n, const float* x, const int *incx, float* y, const int* incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
}//BLAS
struct Session {
std::vector<float> fTensor_batchnormalization1movingvariance0 = std::vector<float>(64);
float * tensor_batchnormalization1movingvariance0 = fTensor_batchnormalization1movingvariance0.data();
std::vector<float> fTensor_batchnormalization1beta0 = std::vector<float>(64);
float * tensor_batchnormalization1beta0 = fTensor_batchnormalization1beta0.data();
std::vector<float> fTensor_batchnormalization1gamma0 = std::vector<float>(64);
float * tensor_batchnormalization1gamma0 = fTensor_batchnormalization1gamma0.data();
std::vector<float> fTensor_dense17kernel0 = std::vector<float>(448);
float * tensor_dense17kernel0 = fTensor_dense17kernel0.data();
std::vector<float> fTensor_dense17bias0 = std::vector<float>(64);
float * tensor_dense17bias0 = fTensor_dense17bias0.data();
std::vector<float> fTensor_batchnormalization1movingmean0 = std::vector<float>(64);
float * tensor_batchnormalization1movingmean0 = fTensor_batchnormalization1movingmean0.data();
std::vector<float> fTensor_batchnormalization1batchnormadd10 = std::vector<float>(64);
float * tensor_batchnormalization1batchnormadd10 = fTensor_batchnormalization1batchnormadd10.data();
std::vector<float> fTensor_dense17BiasAdd0 = std::vector<float>(64);
float * tensor_dense17BiasAdd0 = fTensor_dense17BiasAdd0.data();
std::vector<float> fTensor_dense17bias0bcast = std::vector<float>(64);
float * tensor_dense17bias0bcast = fTensor_dense17bias0bcast.data();


Session(std::string filename ="") {
   if (filename.empty()) filename = "BatchNormalizationtest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_batchnormalization1movingvariance0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_batchnormalization1movingvariance0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_batchnormalization1movingvariance0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_batchnormalization1beta0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_batchnormalization1beta0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_batchnormalization1beta0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_batchnormalization1gamma0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_batchnormalization1gamma0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_batchnormalization1gamma0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense17kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense17kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 448) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 448 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense17kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense17bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense17bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense17bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_batchnormalization1movingmean0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_batchnormalization1movingmean0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_batchnormalization1movingmean0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense17bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense17bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_dense17input){

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
   std::copy(tensor_dense17bias0bcast, tensor_dense17bias0bcast + 64, tensor_dense17BiasAdd0);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense17kernel0, &op_0_ldb, tensor_dense17input, &op_0_lda, &op_0_beta, tensor_dense17BiasAdd0, &op_0_n);
   constexpr int op_1_N =64;
   constexpr int op_1_incx = 1;
   constexpr int op_1_incy = 1;
   BLAS::scopy_(&op_1_N, tensor_dense17BiasAdd0, &op_1_incx,tensor_batchnormalization1batchnormadd10, &op_1_incy);

   float op_1_alpha = -1;
   BLAS::saxpy_(&op_1_N, &op_1_alpha, tensor_batchnormalization1movingmean0, &op_1_incx,tensor_batchnormalization1batchnormadd10, &op_1_incy);

    for (size_t i = 0; i < 64; i++) {
      tensor_batchnormalization1batchnormadd10[i] *= tensor_batchnormalization1gamma0[i] * tensor_batchnormalization1movingvariance0[i]; 
   }
   op_1_alpha = 1;
   BLAS::saxpy_(&op_1_N, &op_1_alpha, tensor_batchnormalization1beta0, &op_1_incx, tensor_batchnormalization1batchnormadd10, &op_1_incy);

   std::vector<float> ret (tensor_batchnormalization1batchnormadd10, tensor_batchnormalization1batchnormadd10 + 64);
   return ret;
}
};
} //TMVA_SOFIE_BatchNormalizationtest

#endif  // ROOT_TMVA_SOFIE_BATCHNORMALIZATIONTEST
