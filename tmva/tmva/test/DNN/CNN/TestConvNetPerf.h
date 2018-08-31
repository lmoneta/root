//
// Created by L. Moneta 
//
/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Im2Col method                                                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Vladimir Ilievsky     <ilievski.vladimir@live.com>  - CERN, Switzerland   *
 *      Manos Stergiadis      <em.stergiadis@gmail.com>  - CERN, Switzerland      *
 *                                                                                *
 * Copyright (c) 2005-2015:                                                       *
 *      CERN, Switzerland                                                         *
 *      U. of Victoria, Canada                                                    *
 *      MPI-K Heidelberg, Germany                                                 *
 *      U. of Bonn, Germany                                                       *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in LICENSE           *
 * (http://tmva.sourceforge.net/LICENSE)                                          *
 **********************************************************************************/

////////////////////////////////////////////////////////////////////
// Functions to test the performances of the CNN                  //
////////////////////////////////////////////////////////////////////
#ifndef ROOT_TESTCONVNETPERF_H
#define ROOT_TESTCONVNETPERF_H

#include "../Utility.h"
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"


using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;




#include <chrono>


template<typename Architecture>
struct ConvNetPerfTest {
   
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;

   size_t batchSize = 32; 
   int nrep = 1000; 
   
 // conv layer output dimension matrix 
   size_t n0 = 12;
   size_t n1 = 32;
   size_t n2 = 32;


   TDeepNet<Architecture> createConvNet(size_t imgDepth, size_t imgHeight, size_t imgWidth,
                                        size_t batchDepth, size_t batchHeight, size_t batchWidth) 
   {


      Net_t convNet(batchSize, imgDepth, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                    ELossFunction::kMeanSquaredError, EInitialization::kGauss);


      // using of size 3 plus padding 1 gives output size = input size 
      size_t depth1 = n0;
      size_t filterHeightConv1 = 3;
      size_t filterWidthConv1 = 3;
      size_t strideRowsConv1 = 1;
      size_t strideColsConv1 = 1;
      size_t zeroPaddingHeight1 = 1;
      size_t zeroPaddingWidth1 = 1;

      EActivationFunction fConv1 = EActivationFunction::kRelu;

      //EActivationFunction fConv1 = ActivationFunctions[rand() % ActivationFunctions.size()];

      // add 2 layers to test also activation backward
      convNet.AddConvLayer(depth1, filterHeightConv1, filterWidthConv1, strideRowsConv1, strideColsConv1, zeroPaddingHeight1,
                       zeroPaddingWidth1, fConv1);
      convNet.AddConvLayer(depth1, filterHeightConv1, filterWidthConv1, strideRowsConv1, strideColsConv1, zeroPaddingHeight1,
                       zeroPaddingWidth1, fConv1);

      //constructConvNet(convNet);
      convNet.Initialize();
      return convNet;
   }


      
   static void createTensor(std::vector<Matrix_t> & X, size_t nb, size_t n, size_t m)
   {
      for (size_t i = 0; i < nb; i++) {
         X.emplace_back(n, m);
         randomMatrix(X[i]);
      }
   }

   // test functions
   
   double testConvForwardPass()
   {

      std::vector<Matrix_t> X;
      createTensor(X, batchSize, n0, n1*n2);

      auto convNet = createConvNet(n0, n1, n2, batchSize, n0, n1*n2); 

      std::chrono::time_point<std::chrono::system_clock> tstart, tend;
      tstart = std::chrono::system_clock::now();
      for (int i= 0; i < nrep; ++i)
         convNet.Forward(X); 
      
      tend = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = tend - tstart;
      return elapsed_seconds.count(); 
   }

   double testConvBackwardPass()
   {

      std::vector<Matrix_t> X;
      createTensor(X, batchSize, n0, n1*n2);

      auto convNet = createConvNet(n0, n1, n2, batchSize, n0, n1*n2);
      // do a forward pass 
      convNet.Forward(X);

      Matrix_t Y(batchSize, convNet.GetOutputWidth());
      Matrix_t W(batchSize, 1);
      randomMatrix(Y);
      randomMatrix(W);


      std::chrono::time_point<std::chrono::system_clock> tstart, tend;
      tstart = std::chrono::system_clock::now();
      for (int i= 0; i < nrep; ++i)
         convNet.Backward(X, Y, W); 
      
      tend = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = tend - tstart;
      return elapsed_seconds.count(); 
   }



};

   // bool testEvalFunction()
   // {

   // }





#endif
