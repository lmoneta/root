//Code generated automatically by TMVA for Inference of Model file [LSTMtestWithBias.h5] at [Thu Aug 24 08:55:51 202] 

#ifndef ROOT_TMVA_SOFIE_LSTMTESTWITHBIAS
#define ROOT_TMVA_SOFIE_LSTMTESTWITHBIAS

#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_LSTMtestWithBias{
namespace BLAS{
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_lstm2lstmcellbias0 = std::vector<float>(24);
float * tensor_lstm2lstmcellbias0 = fTensor_lstm2lstmcellbias0.data();
std::vector<float> fTensor_lstm2lstmcellrecurrentkernel0 = std::vector<float>(36);
float * tensor_lstm2lstmcellrecurrentkernel0 = fTensor_lstm2lstmcellrecurrentkernel0.data();
std::vector<float> fTensor_lstm2lstmcellkernel0 = std::vector<float>(24);
float * tensor_lstm2lstmcellkernel0 = fTensor_lstm2lstmcellkernel0.data();
std::vector<float> fTensor_lstm2PartitionedCall1 = std::vector<float>(6);
float * tensor_lstm2PartitionedCall1 = fTensor_lstm2PartitionedCall1.data();
std::vector<float> fTensor_reshape7Reshape0 = std::vector<float>(4);
float * tensor_reshape7Reshape0 = fTensor_reshape7Reshape0.data();

std::vector<float> fVec_op_1_input = std::vector<float>(4);
std::vector<float> fVec_op_1_initial_hidden_state = std::vector<float>(3);
std::vector<float> fVec_op_1_initial_cell_state = std::vector<float>(3);
std::vector<float> fVec_op_1_ff_input_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_ff_output_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_ff_cell_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_ff_forget_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_input_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_output_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_cell_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_forget_gate = std::vector<float>(6);
std::vector<float> fVec_op_1_cell_state = std::vector<float>(6);
std::vector<float> fVec_op_1_new_cell_state = std::vector<float>(6);
std::vector<float> fVec_op_1_hidden_state = std::vector<float>(6);


Session(std::string filename ="") {
   if (filename.empty()) filename = "LSTMtestWithBias.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_lstm2lstmcellbias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_lstm2lstmcellbias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_lstm2lstmcellbias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_lstm2lstmcellrecurrentkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_lstm2lstmcellrecurrentkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 36) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 36 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_lstm2lstmcellrecurrentkernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_lstm2lstmcellkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_lstm2lstmcellkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_lstm2lstmcellkernel0[i];
   f.close();
}

std::vector<float> infer(float* tensor_reshape7input){
   ///--------Reshape operator

   std::copy( tensor_reshape7input, tensor_reshape7input + 4, tensor_reshape7Reshape0);
   float * op_1_input = fVec_op_1_input.data();
   for(size_t seq = 0; seq < 2; seq++) {
      for(size_t batch = 0; batch < 1; batch++) {
         for(size_t i = 0; i < 2; i++) {
            op_1_input[seq * 2 + batch * 2 + i] = tensor_reshape7Reshape0[batch * 4 + seq * 2 + i];
         }
      }
   }
   float * op_1_ff_input_gate = fVec_op_1_ff_input_gate.data();
   float * op_1_ff_output_gate = fVec_op_1_ff_output_gate.data();
   float * op_1_ff_cell_gate = fVec_op_1_ff_cell_gate.data();
   float * op_1_ff_forget_gate = fVec_op_1_ff_forget_gate.data();
   float * op_1_input_gate = fVec_op_1_input_gate.data();
   float * op_1_output_gate = fVec_op_1_output_gate.data();
   float * op_1_cell_gate = fVec_op_1_cell_gate.data();
   float * op_1_forget_gate = fVec_op_1_forget_gate.data();
   float * op_1_cell_state = fVec_op_1_cell_state.data();
   float * op_1_new_cell_state = fVec_op_1_new_cell_state.data();
   float * op_1_hidden_state = fVec_op_1_hidden_state.data();
   char op_1_transA = 'N';
   char op_1_transB = 'T';
   int op_1_m = 2;
   int op_1_n = 3;
   int op_1_k = 2;
   float  op_1_alpha = 1.;
   float  op_1_beta = 0.;
   int op_1_bias_size = 6;
   int op_1_incx = 1;
   int op_1_incy = 1;
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_lstm2lstmcellkernel0, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_ff_input_gate, &op_1_n);
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_lstm2lstmcellkernel0 + 6, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_ff_output_gate, &op_1_n);
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_lstm2lstmcellkernel0 + 18, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_ff_cell_gate, &op_1_n);
   BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &op_1_m, &op_1_k, &op_1_alpha, tensor_lstm2lstmcellkernel0 + 12, &op_1_k, op_1_input, &op_1_k, &op_1_beta, op_1_ff_forget_gate, &op_1_n);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_lstm2lstmcellbias0, &op_1_incx, op_1_ff_input_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_lstm2lstmcellbias0 + 6, &op_1_incx, op_1_ff_output_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_lstm2lstmcellbias0 + 18, &op_1_incx, op_1_ff_cell_gate, &op_1_incy);
   BLAS::saxpy_(&op_1_bias_size, &op_1_alpha, tensor_lstm2lstmcellbias0 + 12, &op_1_incx, op_1_ff_forget_gate, &op_1_incy);
   for (size_t seq = 0; seq < 2; seq++) {
      size_t ff_offset = seq * 3;
      size_t gate_offset = seq * 3;
      std::copy(op_1_ff_input_gate + ff_offset, op_1_ff_input_gate + ff_offset + 3, op_1_input_gate + gate_offset);
      std::copy(op_1_ff_output_gate + ff_offset, op_1_ff_output_gate + ff_offset + 3, op_1_output_gate + gate_offset);
      std::copy(op_1_ff_cell_gate + ff_offset, op_1_ff_cell_gate + ff_offset + 3, op_1_cell_gate + gate_offset);
      std::copy(op_1_ff_forget_gate + ff_offset, op_1_ff_forget_gate + ff_offset + 3, op_1_forget_gate + gate_offset);
   }
   for (size_t seq = 0; seq < 2; seq++) {
      size_t index = seq;
      int m2 = 1;
      size_t offset = index * 3;
      if (seq == 0) {
      } else {
         size_t previous_offset = (seq - 1) * 3;
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_lstm2lstmcellrecurrentkernel0, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_input_gate + offset, &op_1_n);
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_lstm2lstmcellrecurrentkernel0 + 9, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_output_gate + offset, &op_1_n);
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_lstm2lstmcellrecurrentkernel0 + 27, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_cell_gate + offset, &op_1_n);
         BLAS::sgemm_(&op_1_transB, &op_1_transA, &op_1_n, &m2, &op_1_n, &op_1_alpha, tensor_lstm2lstmcellrecurrentkernel0 + 18, &op_1_n, op_1_hidden_state + previous_offset, &op_1_n, &op_1_alpha, op_1_forget_gate + offset, &op_1_n);
      }
      for (size_t i = offset; i < offset + 3; i++) {
         float ex = exp(-2 * op_1_cell_gate[i]);
            op_1_cell_gate[i] = (1. - ex) / (1. + ex);
      }
      for (size_t i = offset; i < offset + 3; i++) {
            op_1_input_gate[i] = 1. / (1. + exp(-op_1_input_gate[i]));
      }
      for (size_t i = offset; i < offset + 3; i++) {
            op_1_forget_gate[i] = 1. / (1. + exp(-op_1_forget_gate[i]));
      }
      for (size_t i = offset; i < offset + 3; i++) {
         op_1_cell_state[i] = op_1_input_gate[i] * op_1_cell_gate[i];
      }
      if (seq == 0) {
      } else {
         size_t previous_offset = (seq - 1) * 3;
         for (size_t i = 0; i < 3; i++) {
            op_1_cell_state[i + offset] += op_1_forget_gate[i + offset] * op_1_cell_state[i + previous_offset];
         }
      }
      for (size_t i = offset; i < offset + 3; i++) {
            op_1_output_gate[i] = 1. / (1. + exp(-op_1_output_gate[i]));
      }
      std::copy(op_1_cell_state + offset, op_1_cell_state + offset + 3, op_1_new_cell_state + offset);
      for (size_t i = offset; i < offset + 3; i++) {
         float ex = exp(-2 * op_1_new_cell_state[i]);
            op_1_new_cell_state[i] = (1. - ex) / (1. + ex);
      }
      for (size_t i = offset; i < offset + 3; i++) {
         op_1_hidden_state[i] = op_1_output_gate[i] * op_1_new_cell_state[i];
      }
   }
   for (size_t seq = 0; seq < 2; seq++) {
      for (size_t batch = 0; batch < 1; batch++) {
         size_t offset = seq * 3 + 0 + batch * 3;
         size_t y_offset = batch * 6 + seq * 3 + 0;
         std::copy(op_1_hidden_state + offset, op_1_hidden_state + offset + 3, tensor_lstm2PartitionedCall1 + y_offset);
      }
   }
   std::vector<float> ret (tensor_lstm2PartitionedCall1, tensor_lstm2PartitionedCall1 + 6);
   return ret;
}
};
} //TMVA_SOFIE_LSTMtestWithBias

#endif  // ROOT_TMVA_SOFIE_LSTMTESTWITHBIAS
