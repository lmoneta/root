//Code generated automatically by TMVA for Inference of Model file [GRUtestWithBias.h5] at [Thu Aug 24 08:55:50 202] 

#ifndef ROOT_TMVA_SOFIE_GRUTESTWITHBIAS
#define ROOT_TMVA_SOFIE_GRUTESTWITHBIAS

#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_GRUtestWithBias{
namespace BLAS{
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_gru2grucell2bias0 = std::vector<float>(36);
float * tensor_gru2grucell2bias0 = fTensor_gru2grucell2bias0.data();
std::vector<float> fTensor_gru2grucell2recurrentkernel0 = std::vector<float>(27);
float * tensor_gru2grucell2recurrentkernel0 = fTensor_gru2grucell2recurrentkernel0.data();
std::vector<float> fTensor_gru2grucell2kernel0 = std::vector<float>(18);
float * tensor_gru2grucell2kernel0 = fTensor_gru2grucell2kernel0.data();
std::vector<float> fTensor_gru2PartitionedCall1 = std::vector<float>(6);
float * tensor_gru2PartitionedCall1 = fTensor_gru2PartitionedCall1.data();
std::vector<float> fTensor_reshape8Reshape0 = std::vector<float>(4);
float * tensor_reshape8Reshape0 = fTensor_reshape8Reshape0.data();

std::vector<float> fVec_op_1_input = std::vector<float>(4);
std::vector<float> fVec_op_1_initial_hidden_state = std::vector<float>(3);
std::vector<float> fVec_op_1_initial_cell_state = std::vector<float>(3);
std::vector<float> fVec_op_1_f_update_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_f_reset_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_f_hidden_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_update_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_reset_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_hidden_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_feedback = std::vector<float>(3);
std::vector<float> fVec_op_1_hidden_state = std::vector<float>(6);


Session(std::string filename ="") {
   if (filename.empty()) filename = "GRUtestWithBias.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_gru2grucell2bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_gru2grucell2bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 36) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 36 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_gru2grucell2bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_gru2grucell2recurrentkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_gru2grucell2recurrentkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 27) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 27 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_gru2grucell2recurrentkernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_gru2grucell2kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_gru2grucell2kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 18) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 18 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_gru2grucell2kernel0[i];
   f.close();
}

std::vector<float> infer(float* tensor_reshape8input){
   ///--------Reshape operator

   std::copy( tensor_reshape8input, tensor_reshape8input + 4, tensor_reshape8Reshape0);
   float * op_1_input = fVec_op_1_input.data();
   for(size_t seq = 0; seq < 2; seq++) {
      for(size_t batch = 0; batch < 1; batch++) {
         for(size_t i = 0; i < 2; i++) {
            op_1_input[seq * 2 + batch * 2 + i] = tensor_reshape8Reshape0[batch * 4 + seq * 2 + i];
         }
      }
   }
   float * op_1_f_update_gate = fVec_op_1_f_update_gate.data();
   float * op_1_f_reset_gate = fVec_op_1_f_reset_gate.data();
   float * op_1_f_hidden_gate = fVec_op_1_f_hidden_gate.data();
   float * op_1_update_gate = fVec_op_1_update_gate.data();
   float * op_1_reset_gate = fVec_op_1_reset_gate.data();
   float * op_1_hidden_gate = fVec_op_1_hidden_gate.data();
   float * op_1_hidden_state = fVec_op_1_hidden_state.data();
   float * op_1_feedback = fVec_op_1_feedback.data();
   char op_1_transA = 'N';
   char op_1_transB = 'T';
   int op_1_m = 2;
   int op_1_m2 = 1;
   int op_1_n = 3;
   int op_1_k = 2;
   float op_1_alpha = 1.;
   float op_1_beta = 0.;
   int op_1_bias_size = 6;
   int op_1_incx = 1;
   int op_1_incy = 1;
   int op_1_feedback_size = 3;
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_gru2grucell2kernel0, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_f_update_gate, &op_1_n);
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_gru2grucell2kernel0 + 6, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_f_reset_gate, &op_1_n);
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_gru2grucell2kernel0 + 12, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_f_hidden_gate, &op_1_n);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_gru2grucell2bias0, &op_1_incx, op_1_f_update_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_gru2grucell2bias0 + 18, &op_1_incx, op_1_f_update_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_gru2grucell2bias0 + 6, &op_1_incx, op_1_f_reset_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_gru2grucell2bias0 + 24, &op_1_incx, op_1_f_reset_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_gru2grucell2bias0 + 12, &op_1_incx, op_1_f_hidden_gate, &op_1_incy);
   for (size_t seq = 0; seq < 2; seq++) {
      size_t offset = seq * 3;
      size_t gate_offset = seq * 3;
      std::copy(op_1_f_update_gate + offset, op_1_f_update_gate + offset + 3, op_1_update_gate + gate_offset);
      std::copy(op_1_f_reset_gate + offset, op_1_f_reset_gate + offset + 3, op_1_reset_gate + gate_offset);
      std::copy(op_1_f_hidden_gate + offset, op_1_f_hidden_gate + offset + 3, op_1_hidden_gate + gate_offset);
   }
   for (size_t seq = 0; seq < 2; seq++) {
      size_t index = seq;
      int m2 = 1;
      size_t offset = index * 3;
      if (seq == 0) {
      } else {
         size_t previous_offset = (seq - 1) * 3;
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_gru2grucell2recurrentkernel0, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_update_gate + offset, &op_1_n);
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_gru2grucell2recurrentkernel0 + 9, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_reset_gate + offset, &op_1_n);
      }
      for (size_t i = offset; i < offset + 3; i++) {
            op_1_update_gate[i] = 1. / (1. + exp(-op_1_update_gate[i]));
            op_1_reset_gate[i] = 1. / (1. + exp(-op_1_reset_gate[i]));
      }
      if (seq == 0) {
      } else {
         size_t previous_offset = (seq - 1) * 3;
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m2, &op_1_n, &op_1_alpha, tensor_gru2grucell2recurrentkernel0 + 18, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_beta, op_1_feedback, &op_1_n);
      }
      BLAS::saxpy_(&op_1_feedback_size, &op_1_alpha, tensor_gru2grucell2bias0 + 30, &op_1_incx, op_1_feedback, &op_1_incy);
      for (size_t i = 0; i < 3; i++) {
         op_1_feedback[i] *= op_1_reset_gate[i + offset];
      }
      BLAS::saxpy_(&op_1_feedback_size, &op_1_alpha, op_1_feedback, &op_1_incx, op_1_hidden_gate + offset, &op_1_incy);
      for (size_t i = offset; i < offset + 3; i++) {
         float ex = exp(-2 * op_1_hidden_gate[i]);
            op_1_hidden_gate[i] = (1. - ex) / (1. + ex);
      }
      for (size_t i = offset; i < offset + 3; i++) {
         op_1_hidden_state[i] = ( 1. - op_1_update_gate[i]) * op_1_hidden_gate[i];
      }
      if (seq == 0) {
      } else {
         size_t previous_offset = (seq - 1) * 3;
         for (size_t i = 0; i < 3; i++) {
            op_1_hidden_state[i + offset] += op_1_update_gate[i + offset] * op_1_hidden_state[i + previous_offset];
         }
      }
   }
   for (size_t seq = 0; seq < 2; seq++) {
      for (size_t batch = 0; batch < 1; batch++) {
         size_t offset = seq * 3 + 0 + batch * 3;
         size_t y_offset = batch * 6 + seq * 3 + 0;
         std::copy(op_1_hidden_state + offset, op_1_hidden_state + offset + 3, tensor_gru2PartitionedCall1 + y_offset);
      }
   }
   std::vector<float> ret (tensor_gru2PartitionedCall1, tensor_gru2PartitionedCall1 + 6);
   return ret;
}
};
} //TMVA_SOFIE_GRUtestWithBias

#endif  // ROOT_TMVA_SOFIE_GRUTESTWITHBIAS
