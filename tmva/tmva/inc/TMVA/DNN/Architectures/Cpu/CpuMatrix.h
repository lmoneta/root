// @(#)root/tmva/tmva/dnn:$Id$
// Author: Simon Pfreundschuh 20/07/16

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh                                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

//////////////////////////////////////////////////////////
// Definition of the CpuMatrix class used to represent  //
// weight and bias matrices in neural nets.             //
//////////////////////////////////////////////////////////

#ifndef TMVA_DNN_ARCHITECTURES_CPU_CPUMATRIX
#define TMVA_DNN_ARCHITECTURES_CPU_CPUMATRIX

#ifdef R__USE_IMT
#define DL_USE_MTE  // use MT with tbb
#endif

#include <cstddef>
#include <vector>

#include "TMatrix.h"
#include "TMVA/Config.h"
#include "CpuBuffer.h"
#include <TMVA/Config.h>
#include <TMVA/RTensor.hxx>

// #define DEBUG_TMVA_TCPUMATRIX
#if defined(DEBUG_TMVA_TCPUMATRIX)
/*
 * Debug(!) function for printing matrices.
 * 
 * Prints the input expression `mat` using preprocessor directives (with
 * `#mat`). E.g. `PrintMatrix(matA, "Test")` _could_ generate
 * "matA is null pointer".
 *
 * Note: This is a preprocessor macro. It does _not_ respect namespaces.
 * 
 * @param mat  Matrix to print
 * @param text Name of matrix
 */
#define TMVA_DNN_PrintTCpuMatrix(mat, text)                                                                \
{                                                                                                      \
   auto _dpointer = mat.GetRawDataPointer();                                                           \
   if (_dpointer == NULL) {                                                                            \
      std::cout << #mat << " is null pointer" << std::endl;                                            \
      exit(1);                                                                                         \
   }                                                                                                   \
   auto _nrows = mat.GetNrows();                                                                       \
   auto _ncols = mat.GetNcols();                                                                       \
   std::cout << "---------------------" << text << " " << #mat << "(" << _nrows << "," << _ncols << ")"\
             << "--------------------" << std::endl;                                                   \
   for (size_t _i = 0; _i < _nrows; _i++) {                                                            \
      for (size_t _j = 0; _j < _ncols; _j++) {                                                         \
         std::cout << mat(_i, _j);                                                                     \
         if (_j < _ncols - 1) std::cout << ",";                                                        \
      }                                                                                                \
      std::cout << std::endl;                                                                          \
   }                                                                                                   \
}
#else
#define TMVA_DNN_PrintTCpuMatrix(mat, text)
#endif

namespace TMVA
{
namespace DNN
{

/** The TCpuMatrix class.
 *
 * Matrix class for multi-threaded CPU architectures. Uses the TCpuBuffer
 * class to store the matrices in column-major format for compatibility with
 * BLAS. Provides Map and MapFrom member functions to simplify the application of
 * activation functions and derivatives to matrices.
 *
 * Copying and assignment of TCpuMatrix objects only performs shallow copies, i.e.
 * copying is fast and the resulting objects share the element data.
 *
 * \tparam AFloat The floating point type used to represent the matrix elements.
 */
//______________________________________________________________________________
template<typename AFloat>
class TCpuMatrix
{
private:
   static std::vector<AFloat> fOnes;  ///< Vector filled with ones used for BLAS calls.

   TCpuBuffer<AFloat> fBuffer; ///< The buffer holding the matrix elements
                               ///< in column-major format.
   size_t     fNCols;
   size_t     fNRows;

public:

   //friend class TCpuTensor<AFloat>; 

   /** Returns pointer to a vector holding only ones with a guaranteed length
    *  of the number of columns of every instantiated CpuMatrix object. */
   static const AFloat * GetOnePointer() {return fOnes.data();}

   static size_t GetOnePointerSize() { return fOnes.size(); }

   static void InitializeOneVector( size_t n); 

   /** Construct matrix and allocate space for its elements. */
   TCpuMatrix(size_t nRows, size_t nCols);
   /** Construct a TCpuMatrix object by (deeply) copying from a
    *  TMatrixT<Double_t> matrix. */
   TCpuMatrix(const TMatrixT<AFloat> &);
   /** Construct a m-times-n matrix from the given buffer. The size must of
    *  course match. */
   TCpuMatrix(const TCpuBuffer<AFloat> &buffer, size_t m, size_t n);

   //N.B the default copy constructor does a shallow copy (NOT a deep one) !
   TCpuMatrix(const TCpuMatrix  &)             = default;
   TCpuMatrix(      TCpuMatrix &&)             = default;
   TCpuMatrix & operator=(const TCpuMatrix &)  = default;
   TCpuMatrix & operator=(TCpuMatrix &&)       = default;
   ~TCpuMatrix()                               = default;

   /** Clear content of the matrix and initialize to zero elements
    */
   void Zero();

   /** Convert to a TMatrixT<AFloat_t> object. Performs a deep copy of the matrix
    *  elements. */
   operator TMatrixT<AFloat>() const;

   /** Map the given function over the matrix elements. Executed in parallel
    *  using TThreadExecutor. */
   template <typename Function_t>
   void Map(Function_t &f);

   /** Same as maps but takes the input values from the matrix \p A and writes
    *  the results in this matrix. */
   template <typename Function_t>
   void MapFrom(Function_t &f, const TCpuMatrix & A);

   size_t GetNrows() const {return fNRows;}
   size_t GetNcols() const {return fNCols;}
   size_t GetNoElements() const {return fNRows * fNCols;}

   /** Return matrix element in row \p i and column \p j. */
   AFloat   operator()(size_t i, size_t j) const {return fBuffer[j * fNRows + i];}
   AFloat & operator()(size_t i, size_t j)       {return fBuffer[j * fNRows + i];}

   /** Return raw pointer to the elements stored contiguously in column-major
    *  order. */
   AFloat *       GetRawDataPointer()        {return fBuffer;}
   const AFloat * GetRawDataPointer()  const {return fBuffer;}

   static Executor &GetThreadExecutor() { return TMVA::Config::Instance().GetThreadExecutor(); }

    // static function to get the number of elements for task
   static size_t GetNWorkItems(size_t nelements);

   // print matrix
   void Print() const {
      TCpuMatrix cpuMatrix = *this;
      TMVA_DNN_PrintTCpuMatrix(cpuMatrix,"CpuMatrix");
   }
   

private:

   void Initialize();

};

template<typename AFloat>
std::vector<AFloat> TCpuMatrix<AFloat>::fOnes {};


// Inline Functions.
//______________________________________________________________________________
template<typename AFloat>
size_t TCpuMatrix<AFloat>::GetNWorkItems(size_t nElements) 
{
   // nElements should have at least 100 
   // const size_t nWorkers = TMVA::Config::Instance().GetNCpu();
   // return  (nElements > nWorkers) ?  (int) nElements/nWorkers : 1;
   const size_t minElements = 1000;
   const size_t nCpu = TMVA::Config::Instance().GetNCpu();
   if (nElements <= minElements) return nElements;
   if (nElements < nCpu*minElements) {
      size_t nt = nElements/minElements;
      return nElements/nt; 
   }
   return nElements/nCpu; 
   // if (nElements < nCpu*20) return nElements/nCpu;
   // return nElements/(nCpu*10); 
}

   
//______________________________________________________________________________
template<typename AFloat>
template<typename Function_t>
inline void TCpuMatrix<AFloat>::Map(Function_t &f)
{
   AFloat  *data = GetRawDataPointer();
   size_t nelements =  GetNoElements();
   size_t nsteps = TCpuMatrix<AFloat>::GetNWorkItems(nelements);

   auto ff = [data, &nsteps, &nelements, &f](UInt_t workerID)
   {
      size_t jMax = std::min(workerID+nsteps,nelements); 
      for (size_t j = workerID; j < jMax; ++j) {
         data[j] = f(data[j]);
      }
      return 0;
   };

   if (nsteps < nelements) {
      TMVA::Config::Instance().GetThreadExecutor().Foreach(ff, ROOT::TSeqI(0,nelements,nsteps));

      // for (size_t i = 0;  i < nelements; i+=nsteps)
      //    ff(i);

   }
   else {
      R__ASSERT(nelements == nsteps);
      ff(0);
   }
}

//______________________________________________________________________________
template<typename AFloat>
template<typename Function_t>
inline void TCpuMatrix<AFloat>::MapFrom(Function_t &f, const TCpuMatrix &A)
{
         AFloat  *dataB = GetRawDataPointer();
   const AFloat  *dataA = A.GetRawDataPointer();

   size_t nelements =  GetNoElements();
   R__ASSERT(nelements == A.GetNoElements() );
   size_t nsteps = TCpuMatrix<AFloat>::GetNWorkItems(nelements);

   auto ff = [&dataB, &dataA,  &nsteps, &nelements, &f](UInt_t workerID)
   {
      size_t jMax = std::min(workerID+nsteps,nelements); 
      for (size_t j = workerID; j < jMax; ++j) {
         dataB[j] = f(dataA[j]);
      }
      return 0;
   };
   if (nsteps < nelements) { 
      TMVA::Config::Instance().GetThreadExecutor().Foreach(ff, ROOT::TSeqI(0,nelements,nsteps));
      // for (size_t i = 0;  i < nelements; i+=nsteps)
      //    ff(i);

   }
   else {
      R__ASSERT(nelements == nsteps);
      ff(0);
   }
}
//______________________________________________________________________________
template<typename AFloat>
void TCpuMatrix<AFloat>::Zero()  
{
   for (size_t j = 0; j < fNCols; j++) {
      for (size_t i = 0; i < fNRows; i++) {
         (*this)(i, j) = 0;
      }
   }
}


//CPU Tensor Class 
// It is a simple wrapper for TMVA RTensor based on 
// memory owned by CPU Buffer
// We need to keep a pointer for CPUBuffer for fast conversion
// without copying to TCpuMatrix
// also provides compatibility with old interface

template<typename AFloat>
class TCpuTensor { 

private: 

    TCpuBuffer<AFloat> fBuffer; ///< The buffer holding the matrix elements
    TMVA::Experimental::RTensor<AFloat>    fTensor; /// The tensor view 


public:

   friend class TCpuMatrix<AFloat>;

   using Shape_t =  typename TMVA::Experimental::RTensor<AFloat>::Shape_t;
   using MemoryLayout = TMVA::Experimental::MemoryLayout; 

   // default constructor
   //  TCpuTensor() : 
   //    fBuffer(0),
   //    fTensor( fBuffer,  {0} )
   // { }

   /** constructors from batch size, depth, height*width */
   TCpuTensor(size_t bsize, size_t depth, size_t hw) : 
      fBuffer( bsize * depth * hw ),
      fTensor( fBuffer,  { depth, hw, bsize}, MemoryLayout::ColumnMajor )
   { }

   /** constructors from batch size, depth, height, width */
   TCpuTensor(size_t bsize, size_t depth, size_t height, size_t width) : 
      fBuffer( bsize * depth * height * width ),
      fTensor( fBuffer,  { depth, height, width, bsize}, MemoryLayout::ColumnMajor )
   { }
  
   /** constructors from a shape. Allocate the buffer in this case and initialize its elements to zero*/
   TCpuTensor(const Shape_t & shape, MemoryLayout memlayout = MemoryLayout::ColumnMajor ) : 
      fBuffer( TMVA::Experimental::Internal::GetSizeFromShape(shape)),
      fTensor( fBuffer,  shape, memlayout )
   { }

  

   /** constructors from a TCpuBuffer and a rtensor */
   TCpuTensor(const TCpuBuffer<AFloat> &buffer, const  TMVA::Experimental::RTensor<AFloat>  & tensor ) : 
      fBuffer(buffer),
      fTensor( tensor)
   { }

   /** constructors from a TCpuBuffer and a shape */
   TCpuTensor(const TCpuBuffer<AFloat> &buffer, const Shape_t & shape, MemoryLayout memlayout = MemoryLayout::ColumnMajor ) : 
      fBuffer(buffer),
      fTensor( fBuffer, shape, memlayout)
   { }

    /** constructors from a TCpuMatrix */
   TCpuTensor(const TCpuMatrix<AFloat> &matrix) : 
      fBuffer(matrix.fBuffer),
      fTensor( fBuffer, { matrix.GetNrows(), matrix.GetNcols() }, MemoryLayout::ColumnMajor )
   { }

   /** Return raw pointer to the elements stored contiguously in column-major
    *  order. */
   AFloat *       GetRawDataPointer()        {return fBuffer;}
   const AFloat * GetRawDataPointer()  const {return fBuffer;}

   AFloat *       GetData() { return fTensor.GetData();}
   const AFloat * GetData() const { return fTensor.GetData();}

   Shape_t GetShape() const { return fTensor.GetShape(); }

   size_t GetNoElements() const { return fTensor.GetSize(); }
   size_t GetSize() const { return fTensor.GetSize(); }

   // return the size of the first dimension (if in row order) or last dimension if in column order
   size_t GetFirstSize() const { 
      return (fTensor.GetMemoryLayout() == MemoryLayout::RowMajor) ? fTensor.GetShape().first() : fTensor.GetShape().back(); 
   }

   // Matrix conversion for tensors of shape 2 
   TCpuMatrix<AFloat> GetMatrix() const { 
      assert( fTensor.GetShape().size() == 2);
      return TCpuMatrix<AFloat>(fBuffer, fTensor.GetShape()[0], fTensor.GetShape()[1]);
   }
   // copy tensor content in a given matrix. No copy is done, just swapping the pointers
   static void TensorToMatrix(const TMVA::Experimental::RTensor<AFloat> & tensor, TCpuMatrix<AFloat> & matrix) { 
      assert (tensor.GetShape().size() == 2);
      assert( tensor.GetShape()[0] == matrix.GetNRows() );
      assert( tensor.GetShape()[1] == matrix.GetNCols() ); 
      matrix.fBuffer = tensor.fBuffer; 
   }
   // copy matrix content in a given tensor. No copy is done, just swapping the pointers
   static void MatrixToTensor(const TCpuMatrix<AFloat> & matrix, TMVA::Experimental::RTensor<AFloat> & tensor) { 
      if (tensor.GetShape().size() != 2) tensor.fTensor.Squeeze();
      assert (tensor.GetShape().size() == 2);
      assert( tensor.GetShape()[0] == matrix.GetNRows() );
      assert( tensor.GetShape()[1] == matrix.GetNCols() ); 
      tensor.fBuffer = matrix.fBuffer;  
   }

   // return slices in the first dimension (if row wise) or last dimension if colun wise
   // so single event slides 
   TCpuTensor<AFloat> operator()( size_t i ) { 
      typename TMVA::Experimental::RTensor<AFloat>::Slice_t slice( fTensor.GetShape().size()); 
      size_t buffsize = 1;
      if (fTensor.GetMemoryLayout() == MemoryLayout::RowMajor) { 
         slice[0][0] = i; 
         slice[0][1] = i+1; 
         assert(i < fTensor.GetShape()[0] ); 
      
         for (size_t j = 1; j < slice.size(); ++j) { 
            slice[j][0] = 0; 
            slice[j][1] = fTensor.GetShape()[j]; 
            buffsize *= fTensor.GetShape()[j]; 
         }
      } else { 
         slice.back()[0] = i; 
         slice.back()[1] = i+1; 
         assert(i < fTensor.GetShape().back() ); 
      
         for (size_t j = 0; j < slice.size()-1; ++j) { 
            slice[j][0] = 0; 
            slice[j][1] = fTensor.GetShape()[j]; 
            buffsize *= fTensor.GetShape()[j]; 
         }
      }  
      size_t offset = i * buffsize; 
      return TCpuTensor<AFloat>( fBuffer.GetSubBuffer(offset, buffsize)   , fTensor.Slice(slice));
   }
   TCpuTensor<AFloat>  At(size_t i ) { return (*this)(i); }

   // set all the tensor contents to zero
   void Zero() {
      AFloat * data = fBuffer; 
      for (size_t i = 0; i < fBuffer.GetSize(); ++i)
         data[i] = 0; 
   }
   
   // access single element - using first 
   // AFloat   operator()(size_t k, size_t i, size_t j) const {return fBuffer[j * fNRows + i];}
   // AFloat &   operator()(size_t k, size_t i, size_t j) {return fBuffer[j * fNRows + i];}


/** Map the given function over the matrix elements. Executed in parallel
    *  using TThreadExecutor. */
   template <typename Function_t>
   void Map(Function_t &f);

   /** Same as maps but takes the input values from the tensor \p A and writes
    *  the results in this tensor. */
   template <typename Function_t>
   void MapFrom(Function_t &f, const TCpuTensor<AFloat> & A);

};

//______________________________________________________________________________
template<typename AFloat>
template<typename Function_t>
inline void TCpuTensor<AFloat>::Map(Function_t &f)
{
   AFloat  *data = GetRawDataPointer();
   size_t nelements =  GetNoElements();
   size_t nsteps = TCpuMatrix<AFloat>::GetNWorkItems(nelements);

   auto ff = [data, &nsteps, &nelements, &f](UInt_t workerID)
   {
      size_t jMax = std::min(workerID+nsteps,nelements); 
      for (size_t j = workerID; j < jMax; ++j) {
         data[j] = f(data[j]);
      }
      return 0;
   };

   if (nsteps < nelements) {
      TMVA::Config::Instance().GetThreadExecutor().Foreach(ff, ROOT::TSeqI(0,nelements,nsteps));

      // for (size_t i = 0;  i < nelements; i+=nsteps)
      //    ff(i);

   }
   else {
      R__ASSERT(nelements == nsteps);
      ff(0);
   }
}

//______________________________________________________________________________
template<typename AFloat>
template<typename Function_t>
inline void TCpuTensor<AFloat>::MapFrom(Function_t &f, const TCpuTensor<AFloat> &A)
{
         AFloat  *dataB = GetRawDataPointer();
   const AFloat  *dataA = A.GetRawDataPointer();

   size_t nelements =  GetNoElements();
   R__ASSERT(nelements == A.GetNoElements() );
   size_t nsteps = TCpuMatrix<AFloat>::GetNWorkItems(nelements);

   auto ff = [&dataB, &dataA,  &nsteps, &nelements, &f](UInt_t workerID)
   {
      size_t jMax = std::min(workerID+nsteps,nelements); 
      for (size_t j = workerID; j < jMax; ++j) {
         dataB[j] = f(dataA[j]);
      }
      return 0;
   };
   if (nsteps < nelements) { 
      TMVA::Config::Instance().GetThreadExecutor().Foreach(ff, ROOT::TSeqI(0,nelements,nsteps));
      // for (size_t i = 0;  i < nelements; i+=nsteps)
      //    ff(i);

   }
   else {
      R__ASSERT(nelements == nsteps);
      ff(0);
   }
}


} // namespace DNN
} // namespace TMVA

#endif
