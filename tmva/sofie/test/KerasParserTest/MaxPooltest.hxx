//Code generated automatically by TMVA for Inference of Model file [MaxPooltest.h5] at [Thu Aug 24 08:55:47 202] 

#ifndef ROOT_TMVA_SOFIE_MAXPOOLTEST
#define ROOT_TMVA_SOFIE_MAXPOOLTEST

#include<algorithm>
#include<cmath>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_MaxPooltest{
namespace BLAS{
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
}//BLAS
struct Session {
std::vector<float> fTensor_dense2bias0 = std::vector<float>(64);
float * tensor_dense2bias0 = fTensor_dense2bias0.data();
std::vector<float> fTensor_dense2kernel0 = std::vector<float>(64);
float * tensor_dense2kernel0 = fTensor_dense2kernel0.data();
std::vector<float> fTensor_dense2Tanh0 = std::vector<float>(64);
float * tensor_dense2Tanh0 = fTensor_dense2Tanh0.data();
std::vector<float> fTensor_dense2bias0bcast = std::vector<float>(64);
float * tensor_dense2bias0bcast = fTensor_dense2bias0bcast.data();
std::vector<float> fTensor_flatten2Reshape0 = std::vector<float>(1);
float * tensor_flatten2Reshape0 = fTensor_flatten2Reshape0.data();
std::vector<float> fTensor_reshape5Reshape0 = std::vector<float>(4);
float * tensor_reshape5Reshape0 = fTensor_reshape5Reshape0.data();
std::vector<float> fTensor_maxpooling2d2PostTrans = std::vector<float>(1);
float * tensor_maxpooling2d2PostTrans = fTensor_maxpooling2d2PostTrans.data();
std::vector<float> fTensor_dense2Dense = std::vector<float>(64);
float * tensor_dense2Dense = fTensor_dense2Dense.data();
std::vector<float> fTensor_maxpooling2d2MaxPooling2D = std::vector<float>(1);
float * tensor_maxpooling2d2MaxPooling2D = fTensor_maxpooling2d2MaxPooling2D.data();
std::vector<float> fTensor_maxpooling2d2PreTrans = std::vector<float>(4);
float * tensor_maxpooling2d2PreTrans = fTensor_maxpooling2d2PreTrans.data();

std::vector<float> fVec_op_2_xpad = std::vector<float>(4);

Session(std::string filename ="") {
   if (filename.empty()) filename = "MaxPooltest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
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
   if (tensor_name != "tensor_dense2kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense2kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense2kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense2bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense2bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_reshape5input){
   ///--------Reshape operator

   std::copy( tensor_reshape5input, tensor_reshape5input + 4, tensor_reshape5Reshape0);
   ///------- Transpose operator

   for (size_t id = 0; id < 4 ; id++){
      tensor_maxpooling2d2PreTrans[id] = tensor_reshape5Reshape0[ ( id / 4 ) * 4 + ( (id % 4) / 2 ) * 2 + ( (id % 2) ) * 1 + ( (id % 4) / 4 )];
   }

//----  operator MaxPool  op_2
{
   constexpr int hsize = 2;
   constexpr int hmin = 0;
   constexpr int hmax = 1;
   constexpr int kh = 2;
   constexpr int wsize = 2;
   constexpr int wmin = 0;
   constexpr int wmax = 1;
   constexpr int kw = 2;
   size_t outIndex = 0;
   for (size_t n = 0; n < 1; n++) {
      size_t inputOffset = n*4;
      for (int i = hmin; i < hmax; i+=2) {
         for (int j = wmin; j < wmax; j+=2) {
            float value = -INFINITY;
            for (int l = i;  l < i + kh; l++) {
               if (l < 0 || l >= hsize) continue;
               for (int m = j; m < j + kw; m++) {
                  if (m < 0 || m >= wsize) continue;
                     int index = inputOffset + l*wsize + m;
                     auto xval = tensor_maxpooling2d2PreTrans[index];
                     if (xval > value) value = xval;
                  }
               }
            tensor_maxpooling2d2MaxPooling2D[outIndex++] = value;
         }
      }
   }
   }
   ///------- Transpose operator

   for (size_t id = 0; id < 1 ; id++){
      tensor_maxpooling2d2PostTrans[id] = tensor_maxpooling2d2MaxPooling2D[ ( id / 1 ) * 1 + ( (id % 1) ) * 1 + ( (id % 1) / 1 ) * 1 + ( (id % 1) / 1 )];
   }
   ///--------Flatten operator

   std::copy( tensor_maxpooling2d2PostTrans, tensor_maxpooling2d2PostTrans + 1, tensor_flatten2Reshape0);

//--------- Gemm
   char op_5_transA = 'n';
   char op_5_transB = 'n';
   int op_5_m = 1;
   int op_5_n = 64;
   int op_5_k = 1;
   float op_5_alpha = 1;
   float op_5_beta = 1;
   int op_5_lda = 1;
   int op_5_ldb = 64;
   std::copy(tensor_dense2bias0bcast, tensor_dense2bias0bcast + 64, tensor_dense2Dense);
   BLAS::sgemm_(&op_5_transB, &op_5_transA, &op_5_n, &op_5_m, &op_5_k, &op_5_alpha, tensor_dense2kernel0, &op_5_ldb, tensor_flatten2Reshape0, &op_5_lda, &op_5_beta, tensor_dense2Dense, &op_5_n);

//------ TANH
   for (int id = 0; id < 64 ; id++){
      tensor_dense2Tanh0[id] = std::tanh(tensor_dense2Dense[id]);
   }
   std::vector<float> ret (tensor_dense2Tanh0, tensor_dense2Tanh0 + 64);
   return ret;
}
};
} //TMVA_SOFIE_MaxPooltest

#endif  // ROOT_TMVA_SOFIE_MAXPOOLTEST
