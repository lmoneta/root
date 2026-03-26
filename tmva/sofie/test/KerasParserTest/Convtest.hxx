//Code generated automatically by TMVA for Inference of Model file [Convtest.h5] at [Thu Aug 24 08:55:47 202] 

#ifndef ROOT_TMVA_SOFIE_CONVTEST
#define ROOT_TMVA_SOFIE_CONVTEST

#include<cmath>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_Convtest{
namespace BLAS{
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_conv2d2bias0 = std::vector<float>(1);
float * tensor_conv2d2bias0 = fTensor_conv2d2bias0.data();
std::vector<float> fTensor_conv2d2kernel0 = std::vector<float>(4);
float * tensor_conv2d2kernel0 = fTensor_conv2d2kernel0.data();
std::vector<float> fTensor_conv2d2PostTrans = std::vector<float>(3);
float * tensor_conv2d2PostTrans = fTensor_conv2d2PostTrans.data();
std::vector<float> fTensor_conv2d2bias0bcast = std::vector<float>(3);
float * tensor_conv2d2bias0bcast = fTensor_conv2d2bias0bcast.data();
std::vector<float> fTensor_conv2d2Conv2D = std::vector<float>(3);
float * tensor_conv2d2Conv2D = fTensor_conv2d2Conv2D.data();
std::vector<float> fTensor_conv2d2Sigmoid0 = std::vector<float>(3);
float * tensor_conv2d2Sigmoid0 = fTensor_conv2d2Sigmoid0.data();
std::vector<float> fTensor_conv2d2PreTrans = std::vector<float>(4);
float * tensor_conv2d2PreTrans = fTensor_conv2d2PreTrans.data();
std::vector<float> fTensor_reshape2Reshape0 = std::vector<float>(4);
float * tensor_reshape2Reshape0 = fTensor_reshape2Reshape0.data();

std::vector<float> fVec_op_2_f = std::vector<float>(4);
std::vector<float> fVec_op_2_xcol = std::vector<float>(12);


Session(std::string filename ="") {
   if (filename.empty()) filename = "Convtest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_conv2d2bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_conv2d2bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_conv2d2bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_conv2d2kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_conv2d2kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_conv2d2kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_conv2d2bias0, { 1 , 1 , 1 }, { 1 , 1 , 1 , 3 });
      std::copy(data, data + 3, tensor_conv2d2bias0bcast);
      delete[] data;
   }
}

std::vector<float> infer(float* tensor_reshape2input){
   ///--------Reshape operator

   std::copy( tensor_reshape2input, tensor_reshape2input + 4, tensor_reshape2Reshape0);
   ///------- Transpose operator

   for (size_t id = 0; id < 4 ; id++){
      tensor_conv2d2PreTrans[id] = tensor_reshape2Reshape0[ ( id / 4 ) * 4 + ( (id % 4) / 2 ) * 2 + ( (id % 2) ) * 1 + ( (id % 4) / 4 )];
   }

//----  operator Conv op_2
   float * op_2_f = fVec_op_2_f.data();
   for (std::size_t oc = 0; oc < 1; oc++) {
      for (std::size_t ic = 0; ic < 1; ic++) {
         for (std::size_t kh = 0; kh < 2; kh++) {
            for (std::size_t kw = 0; kw < 2; kw++) {
               op_2_f[oc * 4 + ic * 4 + kh * 2 + kw * 1  ] = tensor_conv2d2kernel0[oc * 4 + ic * 4 + kh * 2 + kw ];
            }
         }
      }
   }
   char op_2_transA = 'N';
   char op_2_transB = 'N';
   int op_2_m = 3;
   int op_2_n = 1;
   int op_2_k = 4;
   float op_2_alpha = 1.0;
   float op_2_beta = 0.0;
   float * op_2_xcol = fVec_op_2_xcol.data();
   for (size_t n = 0; n < 1; n++) {
      size_t out_offset = n * 3;
      size_t x_offset = n * 4;
      TMVA::Experimental::SOFIE::UTILITY::Im2col<float>(tensor_conv2d2PreTrans + x_offset,1,2,2,2,2,0,1,1,1,1,1,op_2_xcol);

       BLAS::sgemm_(&op_2_transA, &op_2_transB, &op_2_m, &op_2_n, &op_2_k, &op_2_alpha, op_2_xcol, &op_2_m,
         op_2_f, &op_2_k, &op_2_beta, tensor_conv2d2Conv2D + out_offset, &op_2_m);
   int op_2_size = 3;
   float op_2_gamma = 1.0;
   int op_2_incx = 1;
   int op_2_incy = 1;
   BLAS::saxpy_(&op_2_size, &op_2_gamma, tensor_conv2d2bias0bcast, &op_2_incx, tensor_conv2d2Conv2D + out_offset, &op_2_incy);
   }
   ///------- Transpose operator

   for (size_t id = 0; id < 3 ; id++){
      tensor_conv2d2PostTrans[id] = tensor_conv2d2Conv2D[ ( id / 3 ) * 3 + ( (id % 1) ) * 3 + ( (id % 3) / 3 ) * 3 + ( (id % 3) / 1 )];
   }
	for (int id = 0; id < 3 ; id++){
		tensor_conv2d2Sigmoid0[id] = 1 / (1 + std::exp( - tensor_conv2d2PostTrans[id]));
	}
   std::vector<float> ret (tensor_conv2d2Sigmoid0, tensor_conv2d2Sigmoid0 + 3);
   return ret;
}
};
} //TMVA_SOFIE_Convtest

#endif  // ROOT_TMVA_SOFIE_CONVTEST
