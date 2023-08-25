//Code generated automatically by TMVA for Inference of Model file [CNNtest.h5] at [Thu Aug 24 08:55:45 202] 

#ifndef ROOT_TMVA_SOFIE_CNNTEST
#define ROOT_TMVA_SOFIE_CNNTEST

#include<algorithm>
#include<cmath>
#include<vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_CNNtest{
namespace BLAS{
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
	extern "C" void saxpy_(const int * n, const float * alpha, const float * x,
	                         const int * incx, float * y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_dense1bias0 = std::vector<float>(2);
float * tensor_dense1bias0 = fTensor_dense1bias0.data();
std::vector<float> fTensor_dense1kernel0 = std::vector<float>(128);
float * tensor_dense1kernel0 = fTensor_dense1kernel0.data();
std::vector<float> fTensor_densebias0 = std::vector<float>(64);
float * tensor_densebias0 = fTensor_densebias0.data();
std::vector<float> fTensor_conv2dkernel0 = std::vector<float>(90);
float * tensor_conv2dkernel0 = fTensor_conv2dkernel0.data();
std::vector<float> fTensor_densekernel0 = std::vector<float>(40960);
float * tensor_densekernel0 = fTensor_densekernel0.data();
std::vector<float> fTensor_conv2dbias0 = std::vector<float>(10);
float * tensor_conv2dbias0 = fTensor_conv2dbias0.data();
std::vector<float> fTensor_dense1Sigmoid0 = std::vector<float>(2);
float * tensor_dense1Sigmoid0 = fTensor_dense1Sigmoid0.data();
std::vector<float> fTensor_dense1Dense = std::vector<float>(2);
float * tensor_dense1Dense = fTensor_dense1Dense.data();
std::vector<float> fTensor_denseTanh0 = std::vector<float>(64);
float * tensor_denseTanh0 = fTensor_denseTanh0.data();
std::vector<float> fTensor_denseDense = std::vector<float>(64);
float * tensor_denseDense = fTensor_denseDense.data();
std::vector<float> fTensor_conv2dRelu0 = std::vector<float>(2560);
float * tensor_conv2dRelu0 = fTensor_conv2dRelu0.data();
std::vector<float> fTensor_densebias0bcast = std::vector<float>(64);
float * tensor_densebias0bcast = fTensor_densebias0bcast.data();
std::vector<float> fTensor_flattenReshape0 = std::vector<float>(640);
float * tensor_flattenReshape0 = fTensor_flattenReshape0.data();
std::vector<float> fTensor_reshapeReshape0 = std::vector<float>(256);
float * tensor_reshapeReshape0 = fTensor_reshapeReshape0.data();
std::vector<float> fTensor_maxpooling2dPostTrans = std::vector<float>(640);
float * tensor_maxpooling2dPostTrans = fTensor_maxpooling2dPostTrans.data();
std::vector<float> fTensor_maxpooling2dPreTrans = std::vector<float>(2560);
float * tensor_maxpooling2dPreTrans = fTensor_maxpooling2dPreTrans.data();
std::vector<float> fTensor_conv2dPostTrans = std::vector<float>(2560);
float * tensor_conv2dPostTrans = fTensor_conv2dPostTrans.data();
std::vector<float> fTensor_conv2dConv2D = std::vector<float>(2560);
float * tensor_conv2dConv2D = fTensor_conv2dConv2D.data();
std::vector<float> fTensor_conv2dPreTrans = std::vector<float>(256);
float * tensor_conv2dPreTrans = fTensor_conv2dPreTrans.data();
std::vector<float> fTensor_maxpooling2dMaxPooling2D = std::vector<float>(640);
float * tensor_maxpooling2dMaxPooling2D = fTensor_maxpooling2dMaxPooling2D.data();
std::vector<float> fTensor_dense1bias0bcast = std::vector<float>(2);
float * tensor_dense1bias0bcast = fTensor_dense1bias0bcast.data();
std::vector<float> fTensor_conv2dbias0bcast = std::vector<float>(2560);
float * tensor_conv2dbias0bcast = fTensor_conv2dbias0bcast.data();

std::vector<float> fVec_op_2_f = std::vector<float>(90);
std::vector<float> fVec_op_2_xcol = std::vector<float>(2304);

std::vector<float> fVec_op_6_xpad = std::vector<float>(2560);

Session(std::string filename ="") {
   if (filename.empty()) filename = "CNNtest.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()){
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   int length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense1bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense1bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_dense1bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense1kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense1kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 128) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 128 , read " + std::to_string(length) ;
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
   if (tensor_name != "tensor_conv2dkernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_conv2dkernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 90) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 90 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_conv2dkernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_densekernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_densekernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 40960) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 40960 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_densekernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_conv2dbias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_conv2dbias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 10) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 10 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
    for (int i =0; i < length; ++i) 
       f >> tensor_conv2dbias0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_conv2dbias0, { 10 , 1 , 1 }, { 1 , 10 , 16 , 16 });
      std::copy(data, data + 2560, tensor_conv2dbias0bcast);
      delete[] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_densebias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_densebias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense1bias0,{ 2 }, { 1 , 2 });
      std::copy(data, data + 2, tensor_dense1bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_reshapeinput){
   ///--------Reshape operator

   std::copy( tensor_reshapeinput, tensor_reshapeinput + 256, tensor_reshapeReshape0);
   ///------- Transpose operator

   for (size_t id = 0; id < 256 ; id++){
      tensor_conv2dPreTrans[id] = tensor_reshapeReshape0[ ( id / 256 ) * 256 + ( (id % 256) / 16 ) * 16 + ( (id % 16) ) * 1 + ( (id % 256) / 256 )];
   }

//----  operator Conv op_2
   float * op_2_f = fVec_op_2_f.data();
   for (std::size_t oc = 0; oc < 10; oc++) {
      for (std::size_t ic = 0; ic < 1; ic++) {
         for (std::size_t kh = 0; kh < 3; kh++) {
            for (std::size_t kw = 0; kw < 3; kw++) {
               op_2_f[oc * 9 + ic * 9 + kh * 3 + kw * 1  ] = tensor_conv2dkernel0[oc * 9 + ic * 9 + kh * 3 + kw ];
            }
         }
      }
   }
   char op_2_transA = 'N';
   char op_2_transB = 'N';
   int op_2_m = 256;
   int op_2_n = 10;
   int op_2_k = 9;
   float op_2_alpha = 1.0;
   float op_2_beta = 0.0;
   float * op_2_xcol = fVec_op_2_xcol.data();
   for (size_t n = 0; n < 1; n++) {
      size_t out_offset = n * 2560;
      size_t x_offset = n * 256;
      TMVA::Experimental::SOFIE::UTILITY::Im2col<float>(tensor_conv2dPreTrans + x_offset,1,16,16,3,3,1,1,1,1,1,1,op_2_xcol);

       BLAS::sgemm_(&op_2_transA, &op_2_transB, &op_2_m, &op_2_n, &op_2_k, &op_2_alpha, op_2_xcol, &op_2_m,
         op_2_f, &op_2_k, &op_2_beta, tensor_conv2dConv2D + out_offset, &op_2_m);
   int op_2_size = 2560;
   float op_2_gamma = 1.0;
   int op_2_incx = 1;
   int op_2_incy = 1;
   BLAS::saxpy_(&op_2_size, &op_2_gamma, tensor_conv2dbias0bcast, &op_2_incx, tensor_conv2dConv2D + out_offset, &op_2_incy);
   }
   ///------- Transpose operator

   for (size_t id = 0; id < 2560 ; id++){
      tensor_conv2dPostTrans[id] = tensor_conv2dConv2D[ ( id / 2560 ) * 2560 + ( (id % 10) ) * 256 + ( (id % 2560) / 160 ) * 16 + ( (id % 160) / 10 )];
   }

//------ RELU
   for (int id = 0; id < 2560 ; id++){
      tensor_conv2dRelu0[id] = ((tensor_conv2dPostTrans[id] > 0 )? tensor_conv2dPostTrans[id] : 0);
   }
   ///------- Transpose operator

   for (size_t id = 0; id < 2560 ; id++){
      tensor_maxpooling2dPreTrans[id] = tensor_conv2dRelu0[ ( id / 2560 ) * 2560 + ( (id % 256) / 16 ) * 160 + ( (id % 16) ) * 10 + ( (id % 2560) / 256 )];
   }

//----  operator MaxPool  op_6
{
   constexpr int hsize = 16;
   constexpr int hmin = 0;
   constexpr int hmax = 15;
   constexpr int kh = 2;
   constexpr int wsize = 16;
   constexpr int wmin = 0;
   constexpr int wmax = 15;
   constexpr int kw = 2;
   size_t outIndex = 0;
   for (size_t n = 0; n < 10; n++) {
      size_t inputOffset = n*256;
      for (int i = hmin; i < hmax; i+=2) {
         for (int j = wmin; j < wmax; j+=2) {
            float value = -INFINITY;
            for (int l = i;  l < i + kh; l++) {
               if (l < 0 || l >= hsize) continue;
               for (int m = j; m < j + kw; m++) {
                  if (m < 0 || m >= wsize) continue;
                     int index = inputOffset + l*wsize + m;
                     auto xval = tensor_maxpooling2dPreTrans[index];
                     if (xval > value) value = xval;
                  }
               }
            tensor_maxpooling2dMaxPooling2D[outIndex++] = value;
         }
      }
   }
   }
   ///------- Transpose operator

   for (size_t id = 0; id < 640 ; id++){
      tensor_maxpooling2dPostTrans[id] = tensor_maxpooling2dMaxPooling2D[ ( id / 640 ) * 640 + ( (id % 10) ) * 64 + ( (id % 640) / 80 ) * 8 + ( (id % 80) / 10 )];
   }
   ///--------Flatten operator

   std::copy( tensor_maxpooling2dPostTrans, tensor_maxpooling2dPostTrans + 640, tensor_flattenReshape0);

//--------- Gemm
   char op_9_transA = 'n';
   char op_9_transB = 'n';
   int op_9_m = 1;
   int op_9_n = 64;
   int op_9_k = 640;
   float op_9_alpha = 1;
   float op_9_beta = 1;
   int op_9_lda = 640;
   int op_9_ldb = 64;
   std::copy(tensor_densebias0bcast, tensor_densebias0bcast + 64, tensor_denseDense);
   BLAS::sgemm_(&op_9_transB, &op_9_transA, &op_9_n, &op_9_m, &op_9_k, &op_9_alpha, tensor_densekernel0, &op_9_ldb, tensor_flattenReshape0, &op_9_lda, &op_9_beta, tensor_denseDense, &op_9_n);

//------ TANH
   for (int id = 0; id < 64 ; id++){
      tensor_denseTanh0[id] = std::tanh(tensor_denseDense[id]);
   }

//--------- Gemm
   char op_11_transA = 'n';
   char op_11_transB = 'n';
   int op_11_m = 1;
   int op_11_n = 2;
   int op_11_k = 64;
   float op_11_alpha = 1;
   float op_11_beta = 1;
   int op_11_lda = 64;
   int op_11_ldb = 2;
   std::copy(tensor_dense1bias0bcast, tensor_dense1bias0bcast + 2, tensor_dense1Dense);
   BLAS::sgemm_(&op_11_transB, &op_11_transA, &op_11_n, &op_11_m, &op_11_k, &op_11_alpha, tensor_dense1kernel0, &op_11_ldb, tensor_denseTanh0, &op_11_lda, &op_11_beta, tensor_dense1Dense, &op_11_n);
	for (int id = 0; id < 2 ; id++){
		tensor_dense1Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense1Dense[id]));
	}
   std::vector<float> ret (tensor_dense1Sigmoid0, tensor_dense1Sigmoid0 + 2);
   return ret;
}
};
} //TMVA_SOFIE_CNNtest

#endif  // ROOT_TMVA_SOFIE_CNNTEST
