//Code generated automatically by TMVA for Inference of Model file [SimpleRNNtestWithBias.h5] at [Thu Aug 24 08:55:55 202] 

#ifndef ROOT_TMVA_SOFIE_SIMPLERNNTESTWITHBIAS
#define ROOT_TMVA_SOFIE_SIMPLERNNTESTWITHBIAS

#include<cmath>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_SimpleRNNtestWithBias{
namespace BLAS{
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_simplernn2simplernncellrecurrentkernel0 = std::vector<float>(4);
float * tensor_simplernn2simplernncellrecurrentkernel0 = fTensor_simplernn2simplernncellrecurrentkernel0.data();
std::vector<float> fTensor_simplernn2simplernncellkernel0 = std::vector<float>(4);
float * tensor_simplernn2simplernncellkernel0 = fTensor_simplernn2simplernncellkernel0.data();
std::vector<float> fTensor_simplernn2simplernncellbias0 = std::vector<float>(4);
float * tensor_simplernn2simplernncellbias0 = fTensor_simplernn2simplernncellbias0.data();
std::vector<float> fTensor_simplernn2transpose10 = std::vector<float>(4);
float * tensor_simplernn2transpose10 = fTensor_simplernn2transpose10.data();
std::vector<float> fTensor_reshape4Reshape0 = std::vector<float>(4);
float * tensor_reshape4Reshape0 = fTensor_reshape4Reshape0.data();

std::vector<float> fVec_op_1_input = std::vector<float>(4);
std::vector<float> fVec_op_1_initial_hidden_state = std::vector<float>(2);
std::vector<float> fVec_op_1_feedforward = std::vector<float>(4);
std::vector<float> fVec_op_1_hidden_state = std::vector<float>(4);


Session(std::string filename ="") {
   if (filename.empty()) filename = "SimpleRNNtestWithBias.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_simplernn2simplernncellrecurrentkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_simplernn2simplernncellrecurrentkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_simplernn2simplernncellrecurrentkernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_simplernn2simplernncellkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_simplernn2simplernncellkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_simplernn2simplernncellkernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_simplernn2simplernncellbias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_simplernn2simplernncellbias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_simplernn2simplernncellbias0[i];
   f.close();
}

std::vector<float> infer(float* tensor_reshape4input){
   ///--------Reshape operator

   std::copy( tensor_reshape4input, tensor_reshape4input + 4, tensor_reshape4Reshape0);
   float * op_1_input = fVec_op_1_input.data();
   for(size_t seq = 0; seq < 2; seq++) {
      for(size_t batch = 0; batch < 1; batch++) {
         for(size_t i = 0; i < 2; i++) {
            op_1_input[seq * 2 + batch * 2 + i] = tensor_reshape4Reshape0[batch * 4 + seq * 2 + i];
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
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_simplernn2simplernncellkernel0, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_feedforward, &op_1_n);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_simplernn2simplernncellbias0, &op_1_incx, op_1_feedforward, &op_1_incy);
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
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_simplernn2simplernncellrecurrentkernel0 + r_offset, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_hidden_state + offset, &op_1_n);
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
         std::copy(op_1_hidden_state + offset, op_1_hidden_state + offset + 2, tensor_simplernn2transpose10 + y_offset);
      }
   }
   std::vector<float> ret (tensor_simplernn2transpose10, tensor_simplernn2transpose10 + 4);
   return ret;
}
};
} //TMVA_SOFIE_SimpleRNNtestWithBias

#endif  // ROOT_TMVA_SOFIE_SIMPLERNNTESTWITHBIAS
