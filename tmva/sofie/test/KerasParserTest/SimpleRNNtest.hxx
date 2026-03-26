//Code generated automatically by TMVA for Inference of Model file [SimpleRNNtest.h5] at [Thu Aug 24 08:55:56 202] 

#ifndef ROOT_TMVA_SOFIE_SIMPLERNNTEST
#define ROOT_TMVA_SOFIE_SIMPLERNNTEST

#include<cmath>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_SimpleRNNtest{
namespace BLAS{
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_simplernnsimplernncell2recurrentkernel0 = std::vector<float>(4);
float * tensor_simplernnsimplernncell2recurrentkernel0 = fTensor_simplernnsimplernncell2recurrentkernel0.data();
std::vector<float> fTensor_simplernnsimplernncell2kernel0 = std::vector<float>(4);
float * tensor_simplernnsimplernncell2kernel0 = fTensor_simplernnsimplernncell2kernel0.data();
std::vector<float> fTensor_simplernnsimplernncell2bias0 = std::vector<float>(4);
float * tensor_simplernnsimplernncell2bias0 = fTensor_simplernnsimplernncell2bias0.data();
std::vector<float> fTensor_simplernntranspose10 = std::vector<float>(4);
float * tensor_simplernntranspose10 = fTensor_simplernntranspose10.data();
std::vector<float> fTensor_reshape2Reshape0 = std::vector<float>(4);
float * tensor_reshape2Reshape0 = fTensor_reshape2Reshape0.data();

std::vector<float> fVec_op_1_input = std::vector<float>(4);
std::vector<float> fVec_op_1_initial_hidden_state = std::vector<float>(2);
std::vector<float> fVec_op_1_feedforward = std::vector<float>(4);
std::vector<float> fVec_op_1_hidden_state = std::vector<float>(4);


Session(std::string filename ="") {
   if (filename.empty()) filename = "SimpleRNNtest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_simplernnsimplernncell2recurrentkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_simplernnsimplernncell2recurrentkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_simplernnsimplernncell2recurrentkernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_simplernnsimplernncell2kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_simplernnsimplernncell2kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_simplernnsimplernncell2kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_simplernnsimplernncell2bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_simplernnsimplernncell2bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_simplernnsimplernncell2bias0[i];
   f.close();
}

std::vector<float> infer(float* tensor_reshape2input){
   ///--------Reshape operator

   std::copy( tensor_reshape2input, tensor_reshape2input + 4, tensor_reshape2Reshape0);
   float * op_1_input = fVec_op_1_input.data();
   for(size_t seq = 0; seq < 2; seq++) {
      for(size_t batch = 0; batch < 1; batch++) {
         for(size_t i = 0; i < 2; i++) {
            op_1_input[seq * 2 + batch * 2 + i] = tensor_reshape2Reshape0[batch * 4 + seq * 2 + i];
         }
      }
   }
   float * op_1_feedforward = fVec_op_1_feedforward.data();
   float * op_1_hidden_state = fVec_op_1_hidden_state.data();
   char op_1_transA = 'N';
   char op_1_transB = 'T';
   int op_1_m = 2;
   int op_1_n = 2;
   int op_1_k = 2;
   float op_1_alpha = 1.;
   float op_1_beta = .0;
   int op_1_bias_size = 4;
   int op_1_incx = 1;
   int op_1_incy = 1;
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_simplernnsimplernncell2kernel0, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_feedforward, &op_1_n);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_simplernnsimplernncell2bias0, &op_1_incx, op_1_feedforward, &op_1_incy);
   for (size_t seq = 0; seq < 2; seq++) {
      size_t offset = seq * 2;
      size_t size = 2;
      size_t h_offset = seq * 2 + 0;
      std::copy(op_1_feedforward + offset, op_1_feedforward + offset + size, op_1_hidden_state + h_offset);
   }
   for (size_t seq = 0; seq < 2; seq++) {
      size_t index = seq;
      int m2 = 1;
      size_t offset = index * 2 + 0;
      size_t size = 2;
      if (seq == 0) {
      } else {
         size_t r_offset = 0;
         size_t previous_offset = (seq - 1) * 2 + 0;
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_simplernnsimplernncell2recurrentkernel0 + r_offset, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_hidden_state + offset, &op_1_n);
      }
      for (size_t i = offset; i < offset + size; i++) {
         float ex = std::exp(-2 * op_1_hidden_state[i]);
            op_1_hidden_state[i] = (1. - ex) / (1. + ex);
      }
   }
   for (size_t seq = 0; seq < 2; seq++) {
      for (size_t batch = 0; batch < 1; batch++) {
         size_t offset = seq * 2 + 0 + batch * 2;
         size_t y_offset = batch * 4 + seq * 2 + 0;
         std::copy(op_1_hidden_state + offset, op_1_hidden_state + offset + 2, tensor_simplernntranspose10 + y_offset);
      }
   }
   std::vector<float> ret (tensor_simplernntranspose10, tensor_simplernntranspose10 + 4);
   return ret;
}
};
} //TMVA_SOFIE_SimpleRNNtest

#endif  // ROOT_TMVA_SOFIE_SIMPLERNNTEST
