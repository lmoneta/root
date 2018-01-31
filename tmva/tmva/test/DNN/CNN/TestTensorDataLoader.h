// @(#)root/tmva/tmva/cnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  :                                                                       *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Testing Tensor Data Loader Features                                       *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Vladimir Ilievski      <ilievski.vladimir@live.com>  - CERN, Switzerland  *
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

#ifndef TMVA_TEST_DNN_TEST_CNN_TEST_TENSOR_DATA_LOADER_H
#define TMVA_TEST_DNN_TEST_CNN_TEST_TENSOR_DATA_LOADER_H

#include "../Utility.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TBranch.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "TMVA/Event.h"

#include "TMVA/DNN/TensorDataLoader.h"
#include "TMVA/DNN/DeepNet.h"

#include <vector>
#include <iostream>

using namespace TMVA::DNN;
using namespace TMVA::DNN::CNN;

/** Test that the tensor data loader loads all data in the data set by summing
 *  up all elements batch wise and comparing to the result obtained by summing
 *  over the complete dataset.
 */
//______________________________________________________________________________
template <typename Architecture_t>
auto testSum() -> typename Architecture_t::Scalar_t
{
   using Scalar_t = typename Architecture_t::Scalar_t;
   using Matrix_t = typename Architecture_t::Matrix_t;
   using DataLoader_t = TTensorDataLoader<TensorInput, Architecture_t>;

   size_t nSamples = 10000;
   size_t nRows = 1000;
   size_t nCols = 1;

   // Create the tensor input
   std::vector<TMatrixT<Double_t>> inputTensor;
   std::vector<Matrix_t> inputTensorArch;
   inputTensor.reserve(nSamples);
   inputTensorArch.reserve(nSamples);

   for (size_t i = 0; i < nSamples; i++) {
      inputTensor.emplace_back(nRows, nCols);
      inputTensorArch.emplace_back(nRows, nCols);
   }

   for (size_t i = 0; i < nSamples; i++) {
      for (size_t j = 0; j < nRows; j++) {
         for (size_t k = 0; k < nCols; k++) {
            inputTensor[i](j, k) = i;
            inputTensorArch[i](j, k) = i;
         }
      }
   }

   // Create the output
   TMatrixT<Double_t> output(nSamples, 1);
   Matrix_t outputArch(nSamples, 1);
   for (size_t i = 0; i < nSamples; i++) {
      output(i, 0) = i;
      outputArch(i, 0) = i;
   }

   // Create the weights
   TMatrixT<Double_t> weights(nSamples, 1);
   Matrix_t weightsArch(nSamples, 1);
   for (size_t i = 0; i < nSamples; i++) {
      weights(i, 0) = i;
      weightsArch(i, 0) = i;
   }

   TensorInput input(inputTensor, output, weights);
   DataLoader_t loader(input, nSamples, 5, 5, 1000, 1, 1);

   Matrix_t Sum(1, 1), SumTotal(1, 1);
   Scalar_t sum = 0.0, sumTotal = 0.0;

   for (auto b : loader) {
      for (size_t i = 0; i < b.GetInput().size(); i++) {
         Architecture_t::SumColumns(Sum, b.GetInput()[i]);
         sum += Sum(0, 0);
      }

      Architecture_t::SumColumns(Sum, b.GetOutput());
      sum += Sum(0, 0);

      Architecture_t::SumColumns(Sum, b.GetWeights());
      sum += Sum(0, 0);
   }

   for (size_t i = 0; i < inputTensorArch.size(); i++) {
      Architecture_t::SumColumns(SumTotal, inputTensorArch[i]);
      sumTotal += SumTotal(0, 0);
   }

   Architecture_t::SumColumns(SumTotal, outputArch);
   sumTotal += SumTotal(0, 0);

   Architecture_t::SumColumns(SumTotal, weightsArch);
   sumTotal += SumTotal(0, 0);

   return fabs(sumTotal - sum) / sumTotal;
}

//______________________________________________________________________________
template <typename Architecture_t>
auto testIdentity() -> typename Architecture_t::Scalar_t
{
   using Scalar_t = typename Architecture_t::Scalar_t;
   using Net_t = TDeepNet<Architecture_t>;
   using DataLoader_t = TTensorDataLoader<TensorInput, Architecture_t>;

   size_t nSamples = 2000;
   size_t nRows = 3;
   size_t nCols = 1024;

   std::vector<TMatrixT<Double_t>> inputTensor;
   inputTensor.reserve(nSamples);

   for (size_t i = 0; i < nSamples; i++) {
      inputTensor.emplace_back(nRows, nCols);
   }

   for (size_t i = 0; i < nSamples; i++) {
      randomMatrix(inputTensor[i]);
   }

   // Create the output
   TMatrixT<Double_t> output(nSamples, 1);
   for (size_t i = 0; i < nSamples; i++) {
      output(i, 0) = 0;
   }

   // Create the weights
   TMatrixT<Double_t> weights(nSamples, 1);
   for (size_t i = 0; i < nSamples; i++) {
      weights(i, 0) = 1;
   }

   TensorInput input(inputTensor, output, weights);
   DataLoader_t loader(input, nSamples, 5, 5, 3, 1024, 1);

   Net_t convNet(5, 3, 32, 32, 5, 3, 1024, ELossFunction::kMeanSquaredError, EInitialization::kIdentity);
   constructLinearConvNet(convNet);
   convNet.Initialize();

   Scalar_t maximumError = 0.0;
   for (auto b : loader) {
      auto linputTensor = b.GetInput();
      auto outputMatrix = b.GetOutput();
      auto weightMatrix = b.GetWeights();
      Scalar_t error = convNet.Loss(linputTensor, outputMatrix, weightMatrix);
      maximumError = std::max(error, maximumError);
   }

   return maximumError;
}

//______________________________________________________________________________
template <typename Architecture_t>
void testDataLoaderDataSet()
{
   using DeepNet_t = TDeepNet<Architecture_t>;
   using DataLoader_t = TTensorDataLoader<TensorInput, Architecture_t>;

   // Load the photon input
   TFile *photonInput(0);
   TString photonFileName =
      "/Users/vladimirilievski/Desktop/Vladimir/GSoC/ROOT-CI/common-version/root/tmva/tmva/test/DNN/CNN/"
      "dataset/SinglePhotonPt50_FEVTDEBUG_n250k_IMG_CROPS32.root";
   photonInput = TFile::Open(photonFileName);
   TTreeReader photonReader("RHTree", photonInput);

   // Load the electron input
   TFile *electronInput(0);
   TString electronFileName = "/Users/vladimirilievski/Desktop/Vladimir/GSoC/ROOT-CI/common-version/root/tmva/tmva/"
                              "test/DNN/CNN/dataset/SingleElectronPt50_FEVTDEBUG_n250k_IMG_CROPS32.root";
   electronInput = TFile::Open(electronFileName);
   TTreeReader electronReader("RHTree", electronInput);

   // Calculate the number of samples
   size_t nSamples = (size_t)(photonReader.GetEntries(true) + electronReader.GetEntries(true));
   size_t nChannels = 1;
   size_t imgWidth = 32;
   size_t imgHeight = 32;
   size_t batchSize = 256;

   // Initialize the input
   std::vector<TMatrixT<Double_t>> inputTensor;
   inputTensor.reserve(nSamples);
   for (size_t i = 0; i < nSamples; i++) {
      inputTensor.emplace_back(nChannels, imgWidth * imgHeight);
   }

   size_t nOutput = 2;
   // Initialize the output and the weights
   TMatrixT<Double_t> output(nSamples, nOutput);
   TMatrixT<Double_t> weights(nSamples, 1);

   // Get the Energy Red Branches
   TTreeReaderArray<Real_t> photonEnergyRedBr(photonReader, "EBenergyRed");
   TTreeReaderArray<Real_t> electronEnergyRedBr(electronReader, "EBenergyRed");

   size_t counter = 0;

   // Read the photons first
   while (photonReader.Next()) {
      // Set the input
      for (size_t i = 0; i < photonEnergyRedBr.GetSize(); i++) {
         inputTensor[counter](0, i) = static_cast<Double_t>(photonEnergyRedBr[i]);
      }

      // Set the output
      output(counter, 0) = 1.0;
      output(counter, 1) = 0.0;

      // Set the weights
      weights(counter, 0) = 1.0;

      counter++;
   }

   // Read the electrons
   while (electronReader.Next()) {
      // Set the input
      for (size_t i = 0; i < electronEnergyRedBr.GetSize(); i++) {
         inputTensor[counter](0, i) = static_cast<Double_t>(electronEnergyRedBr[i]);
      }

      // Set the output
      output(counter, 0) = 0.0;
      output(counter, 1) = 1.0;

      // Set the weights
      weights(counter, 0) = 1.0;

      counter++;
   }

   // Construct the Linear Net
   size_t batchDepth = batchSize;
   size_t batchHeight = nChannels;
   size_t batchWidth = imgHeight * imgWidth;
   size_t nOutputs = 2;

   TensorInput input(inputTensor, output, weights);
   DataLoader_t loader(input, nSamples, batchSize, batchDepth, batchHeight, batchWidth, nOutputs);

   DeepNet_t convNet(batchSize, nChannels, imgHeight, imgWidth, batchDepth, batchHeight, batchWidth,
                     ELossFunction::kMeanSquaredError, EInitialization::kIdentity);
   constructConvNet(convNet);
   convNet.Initialize();

   for (auto b : loader) {
      auto linputTensor = b.GetInput();
      auto outputMatrix = b.GetOutput();
      auto weightMatrix = b.GetWeights();

      convNet.Forward(linputTensor);
      convNet.Backward(linputTensor, outputMatrix, weightMatrix);
   }
}

#endif
