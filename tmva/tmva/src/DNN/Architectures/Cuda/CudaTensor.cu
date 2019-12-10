// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 13/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/////////////////////////////////////////////
// Implementation of the TCudaTensor class. //
/////////////////////////////////////////////

#include "TMVA/DNN/Architectures/Cuda/CudaTensor.h"
#include "TMVA/DNN/Architectures/Cuda/Device.h"

#include <algorithm>
#include <cassert>
#include <iostream>

namespace TMVA {
namespace DNN  {


// Static members.
//____________________________________________________________________________
/*template<typename AFloat>
size_t                   TCudaTensor<AFloat>::fInstances        = 0;*/
/*template<typename AFloat>
cublasHandle_t           TCudaTensor<AFloat>::fCublasHandle     = nullptr;*/
/*template<typename AFloat>
cudnnHandle_t            TCudaTensor<AFloat>::fCudnnHandle      = nullptr;*/
template<typename AFloat>
std::vector<cudnnHandle_t> TCudaTensor<AFloat>::fCudnnHandle(1);
/*template<typename AFloat>
cudnnTensorDescriptor_t  TCudaTensor<AFloat>::fTensorDescriptor = nullptr;*/
template<typename AFloat>
cudnnDataType_t          TCudaTensor<AFloat>::fDataType         = CUDNN_DATA_FLOAT;
/*template<typename AFloat>
AFloat                   * TCudaTensor<AFloat>::fDeviceReturn   = nullptr;*/
/*template<typename AFloat>
AFloat                   * TCudaTensor<AFloat>::fOnes           = nullptr;*/
/*template<typename AFloat>
curandState_t            * TCudaTensor<AFloat>::fCurandStates   = nullptr;*/
/*template<typename AFloat>
size_t                   TCudaTensor<AFloat>::fNCurandStates    = 0;*/
/*template<typename AFloat>
size_t                   TCudaTensor<AFloat>::fNOnes            = 0;*/
/*template<typename AFloat>
std::vector<std::vector<int> >         TCudaTensor<AFloat>::fStreamIndxs(std:vector<int>(), std::vector<int>());*/
template<typename AFloat>
std::vector<int>         TCudaTensor<AFloat>::fInstances(1,0);

/// This information is needed for the multi-dimensional indexing. See here:
/// https://en.wikipedia.org/wiki/Row-_and_column-major_order
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.strides.html
template<typename AFloat>
std::vector<std::size_t> TCudaTensor<AFloat>::ComputeStridesFromShape(const std::vector<std::size_t> &shape, 
   bool rowmajorLayout)
{
   const auto size = shape.size();
   std::vector<std::size_t> strides(size);
   if (rowmajorLayout)  {
      for (std::size_t i = 0; i < size; i++) {
         if (i == 0) {
            strides[size - 1 - i] = 1;
         } else {
            strides[size - 1 - i] = strides[size - 1 - i + 1] * shape[size - 1 - i + 1];
         }
      }
   } else  {
      for (std::size_t i = 0; i < size; i++) {
         if (i == 0) {
            strides[i] = 1;
         } else {
            strides[i] = strides[i - 1] * shape[i - 1];
         }
      }
   }
   return strides;
}

// Constructors.
//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor()
    : fShape(), fStrides(), fNDim(0), fSize(0), fElementBuffer(), fStreamIndx(0), fTensorDescriptor(nullptr)
{
   InitializeCuda();
}


//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const std::vector<size_t> & shape,
                                 TCudaTensor::MemoryLayout layout,
                                 int device, int streamIndx)
    : fShape(shape), fStrides(shape.size()), fNDim(shape.size()), fDevice(device), fStreamIndx(streamIndx),
      fTensorDescriptor(nullptr), fMemoryLayout(layout)
{
   fStrides = ComputeStridesFromShape(fShape, layout==MemoryLayout::RowMajor);
   
   fSize = (layout==MemoryLayout::RowMajor) ? fStrides.front()*fShape.front() : 
                                              fStrides.back()*fShape.back(); 

   fElementBuffer = TCudaDeviceBuffer<AFloat>(fSize, 0);
   
   InitializeCuda();
}

//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const AFloat * host_data, const std::vector<size_t> & shape,
                                 TCudaTensor::MemoryLayout layout,
                                 int device, int streamIndx)
   : TCudaTensor(shape, layout, device, streamIndx)
{
   // do I need to allocate this buffer ???? 
   // is not a mem leak
   // AFloat * buffer = new AFloat[fSize];
   // size_t index = 0;
   // for (size_t j = 0; j < fSize; ++j) {
   //       buffer[j] = static_cast<AFloat>(host_data[j]);
   //    }
   // }

   cudaMemcpy(fElementBuffer, host_data, fSize * sizeof(AFloat),
              cudaMemcpyHostToDevice);
}

//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(TCudaDeviceBuffer<AFloat> buffer,  
                                 const std::vector<size_t> & shape,
                                 TMVA::Experimental::MemoryLayout layout,
                                 int device, int streamIndx)
   : fNDim(shape.size()), fElementBuffer(buffer), fShape(shape), fStrides( shape.size()), fDevice(device), 
     fStreamIndx(streamIndx), fTensorDescriptor(nullptr), fMemoryLayout(layout)
{
   fStrides = ComputeStridesFromShape(fShape, layout==MemoryLayout::RowMajor);
   
   fSize = (layout==MemoryLayout::RowMajor) ? fStrides.front()*fShape.front() : 
                                              fStrides.back()*fShape.back();  
   InitializeCuda();  
}

//____________________________________________________________________________
//FIXME: Go to shared_ptr implementation of instance tracking
template <typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const TCudaTensor<AFloat>& oldTensor, size_t /*dim*/) :
   TCudaTensor(oldTensor.fShape, oldTensor.fMemoryLayout, oldTensor.fDevice, oldTensor.fStreamIndx)
{
   // No deep copy
   fStrides       = oldTensor.fStrides;
   fElementBuffer = oldTensor.fElementBuffer;   
        
   InitializeCuda();
}

//____________________________________________________________________________
template <typename AFloat>
TCudaTensor<AFloat>::TCudaTensor(const TCudaMatrix<AFloat>& matrix, size_t dim) :
   TCudaTensor( matrix.GetDeviceBuffer(), {matrix.GetNrows(), matrix.GetNcols()}, MemoryLayout::ColumnMajor)
{
   // No deep copy
   if (dim > 2) {
      // change shape from (nrows,ncols) to (nrows,ncols,1,1)
      // this works onlt for coolum major layout since this is same of TCudaMatrix
      fShape.insert(fShape.end(), dim-2, 1);
      fStrides.insert(fStrides.end(),dim-2,fSize); 
      fNDim = dim; 
   }

   InitializeCuda();
}

//____________________________________________________________________________
template <typename AFloat>
TCudaTensor<AFloat>::~TCudaTensor() 
{
   if (fTensorDescriptor && fTensorDescriptor.use_count() == 1 ) { 
      CUDNNCHECK(cudnnDestroyTensorDescriptor(fTensorDescriptor->fCudnnDesc));
   }
   fInstances[fStreamIndx]--;

   // When all tensors in a streamIndx are destroyed, release cudnn resources 
   //if (--fInstances[fStreamIndx] <= 0) CUDNNCHECK(cudnnDestroy(fCudnnHandle[fStreamIndx]));
//#endif
}

template<typename AFloat>
TCudaTensor<AFloat>::operator TMatrixT<AFloat>() const
{
   // this should work only for size 2 or 4 tensors
   if (fNDim < 4) {
      TCudaMatrix<AFloat> temp = GetMatrix();
      return temp;
   }
   // we can convert directy to TMatrix 
   assert(fNDim == 4); 
   size_t nRows = fShape[0]*fShape[1];
   size_t nCols = fShape[2]*fShape[3]; 
   TMatrixT<AFloat> hostMatrix( nRows, nCols ); 

   
   cudaMemcpy(hostMatrix.GetMatrixArray(), fElementBuffer, fSize * sizeof(AFloat),
           cudaMemcpyDeviceToHost);

   return hostMatrix;
}
//____________________________________________________________________________
template <typename AFloat>
inline void TCudaTensor<AFloat>::InitializeCuda()
{
   // If fNDim >= 4, a cudnn tensor is required and we initialize cudnn
   if (fNDim >= 2 && fSize > 0) {
      // Also check whether a new streamIndx has been opened
      if (fInstances.size() - 1 < fStreamIndx) {
         // If need to resize once, need probably to resize more often
         fInstances.resize(2*fStreamIndx + 1, 0);
         fCudnnHandle.resize(2*fStreamIndx + 1, nullptr);
      }
      if (fInstances[fStreamIndx] == 0) {
        CUDNNCHECK(cudnnCreate(&fCudnnHandle[fStreamIndx]));
        //CUDNNCHECK(cudnnSetStream(fCudnnHandle[fStreamIndx], fElementBuffer.GetComputeStream()));
        
        //cublasCreate(&fCublasHandle);
        //CUDACHECK(cudaMalloc(& fDeviceReturn, sizeof(AFloat)));
        //CUDACHECK(cudaMalloc(& fCurandStates, TDevice::NThreads(*this)));
      }
      // if (TDevice::NThreads(*this) > (int) fNCurandStates) {
      //     fNCurandStates = TDevice::NThreads(*this);
      //     if (fCurandStates) {
      //         cudaFree(fCurandStates);
      //     }
      //     cudaMalloc(&fCurandStates, TDevice::NThreads(*this) * sizeof(curandState_t));
      //     InitializeCurandStates();
      // }

      if (!fTensorDescriptor) { 
         fTensorDescriptor = std::make_shared<TensorDescriptor>();
         CUDNNCHECK(cudnnCreateTensorDescriptor(&(fTensorDescriptor->fCudnnDesc)));
      }
        
      fInstances[fStreamIndx]++;
      
      // Prevent template specialization of entire class
      if      (std::is_same<AFloat, double>::value) {fDataType = CUDNN_DATA_DOUBLE;}
      else if (std::is_same<AFloat, float>::value)  {fDataType = CUDNN_DATA_FLOAT;}

      SetTensorDescriptor();
   }

}
template<typename AFloat>
void TCudaTensor<AFloat>::SetTensorDescriptor() {
      if (!fTensorDescriptor) return; 
      if (fSize == 0) return;

      // cuDNN NdTensor format has a minsize of 4 tensor dimensions
      // 4D tensor is more performant on lower dimensions and supports all folowing operations
      //if (fNDim == 4) {
      Shape_t shape = fShape; 
      if (fNDim < 4 && fNDim > 1 ) { 
         // add 1 to tensor 
         if (fMemoryLayout == MemoryLayout::RowMajor)  
            shape.insert(shape.end(),4-fNDim, 1);
         else 
            shape.insert(shape.begin(),4-fNDim,1);
      } else if (fNDim > 4) { 
         std::cout << "Error : Dim = "<< fNDim 
         <<". Currently only 4D tensors are supported for the TMVA cuDNN backend."
         << std::endl;
      }
      if (fMemoryLayout == MemoryLayout::RowMajor)  {
            CUDNNCHECK(cudnnSetTensor4dDescriptor(fTensorDescriptor->fCudnnDesc,
                                               CUDNN_TENSOR_NCHW,// Layout of the tensor in memory
                                               fDataType,
                                               (int)shape[0],   // batch size
                                               (int)shape[1],   // no. channels
                                               (int)shape[2],   // image height
                                               (int)shape[3])); // image width
      }
      else {
            CUDNNCHECK(cudnnSetTensor4dDescriptor(fTensorDescriptor->fCudnnDesc,
                       CUDNN_TENSOR_NCHW,// Layout of the tensor in memory
                       fDataType,
                       (int)shape[3],   // batch size
                       (int)shape[2],   // no. channels
                       (int)shape[1],   // image height
                       (int)shape[0])); // image width
      }
      
      // Some operations in cudnn may not work with this tensor description
      //else if 
 
        /*CUDNNCHECK(cudnnSetTensorNdDescriptor(fTensorDescriptor,
                                              fDataType,
                                              (int)fNDim,
                                              (int *)fShape.data(),
                                              (int *)fStrides.data()));*/
      //}
   
#ifdef NDEBUG
      size_t tensorSize;
      CUDNNCHECK(cudnnGetTensorSizeInBytes(fTensorDescriptor->fCudnnDesc, &tensorSize));
      assert(fSize == tensorSize/sizeof(AFloat));

        //    int n,c,h,w = 0; 
   // int s1,s2,s3,s4 = 0; 
   // cudnnDataType_t  dataType; 
   // cudnnGetTensor4dDescriptor( fTensorDescriptor, &dataType,&n,&c,&h,&w,&s1,&s2,&s3,&s4 );
   // std::vector<size_t>  shape_input = {n,c,h,w}; 
   // assert (shape_input == GetShape());

#endif

 
   }

//____________________________________________________________________________
template<typename AFloat>
void TCudaTensor<AFloat>::InitializeCurandStates()
{
   // dim3 blockDims = TDevice::BlockDims2D();
   // dim3 gridDims  = TDevice::GridDims2D(*this);
   // CurandInitializationKernel<<<gridDims, blockDims>>>(time(nullptr), fCurandStates);
}

template<typename AFloat>
void TCudaTensor<AFloat>::Print(const char * name, bool truncate) const
{
      //TCudaBuffer<AFloat> hostBuffer (fSize);
      //fElementBuffer.CopyTo(hostBuffer);
    #if 0  
      AFloat hostBuffer[fSize]; 

      cudaMemcpy(hostBuffer, fElementBuffer, fSize * sizeof(AFloat),
                 cudaMemcpyDeviceToHost);
      
      for (size_t i = 0; i < fSize; i++) std::cout << hostBuffer[i] << "  ";
   #endif
   PrintShape(name);
   size_t n = fSize; 
   if (n > 10 && truncate) n = 10; 
   std::cout << "Data : { ";
   for (size_t i = 0; i < n; ++i ) {
      AFloat * elementPointer = fElementBuffer + i; 
      std::cout << AFloat( TCudaDeviceReference<AFloat>(elementPointer) );
      if (i < n-1) std::cout << " , "; 
   }
   if (n < fSize) std::cout << "............   } "; 
   std::cout << " } " << std::endl;
}
template<typename AFloat>
void TCudaTensor<AFloat>::PrintShape(const char * name) const
{
      std::string memlayout = (GetLayout() == MemoryLayout::RowMajor) ? "RowMajor" : "ColMajor"; 
      std::cout << name << " shape : { ";
      for (size_t i = 0; i < fNDim-1; ++i ) 
         std::cout << fShape[i] << " , ";
      std::cout << fShape.back() << " } " << " Layout : " << memlayout << std::endl;
}
#if 0
// Conversion to RTensor
//____________________________________________________________________________
template<typename AFloat>
TCudaTensor<AFloat>::operator Experimental::RTensor<AFloat>() const
{
   std::vector<size_t> shape(fNDims, fNDims + fDim)
   
   Experimental::RTensor<AFloat> hostTensor( shape)

   AFloat * buffer = new AFloat[fSize];
   cudaMemcpy(buffer, fElementBuffer, fSize * sizeof(AFloat),
              cudaMemcpyDeviceToHost);

   int index = 0;
   for (int j = 0; j < fSize; j++) {
         hostTensor.GetData()[j] = static_cast<AFloat>(buffer[j]);
      }
   }

   delete[] buffer;
   return hostTensor;
}
#endif
// Explicit Instantiations.

template class TCudaTensor<float>;
template class TCudaTensor<double>;

} // namespace DNN
} // namespace TMVA
