// @(#)root/tmva $Id$
// Author: Simon Pfreundschuh

/*************************************************************************
 * Copyright (C) 2016, Simon Pfreundschuh, Saurav Shekhar                *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

////////////////////////////////////////////////////////////////////
// Generic tests of the backpropagation algorithm.                //
//                                                                //
// All tests randomly generate a net with identity activation     //
// functions, i.e.  which is completely linear and then tests the //
// computed gradients for each layer using numerical              //
// derivation. The restriction to linear nets is to avoid the     //
// required division by the finite difference interval used to    //
// approximate the numerical derivatives, which would otherwise   //
// cause precision loss.                                          //
////////////////////////////////////////////////////////////////////

#include <iostream>
#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"
#include "TMVA/DNN/BatchNormLayer.h"
#include "Utility.h"

using namespace TMVA::DNN;

/*! Compute the loss of the net as a function of the weight at index (i,j) in
 *  layer l. dx is added as an offset to the current value of the weight. */
//______________________________________________________________________________
template <typename Architecture>
auto evaluate_net_weight(TDeepNet<Architecture> &net, std::vector<typename Architecture::Matrix_t> & X,
                         const typename Architecture::Matrix_t &Y, const typename Architecture::Matrix_t &W, size_t l,
                         size_t k, size_t i, size_t j, typename Architecture::Scalar_t dx) ->
   typename Architecture::Scalar_t
{
    using Scalar_t = typename Architecture::Scalar_t;

    net.GetLayerAt(l)->GetWeightsAt(k).operator()(i,j) += dx;
    Scalar_t res = net.Loss(X, Y, W, false, false);
    net.GetLayerAt(l)->GetWeightsAt(k).operator()(i,j) -= dx;
    return res;
}


   
// TODO pass as function params
size_t tbatchSize = 2, timeSteps = 1, inputSize = 2, outputSize = 2;

/*! Generate a random net, perform forward and backward propagation and check
 *  the weight gradients using numerical differentiation. Returns the maximum
 *  relative gradient error and also prints it to stdout. */
//______________________________________________________________________________
template <typename Architecture>
auto testBackpropagationWeights(typename Architecture::Scalar_t dx)
-> typename Architecture::Scalar_t
{
   using Scalar_t = typename Architecture::Scalar_t;
   using Matrix_t = typename Architecture::Matrix_t;
   using Net_t = TDeepNet<Architecture>;
   // using FCLayer_t  = TDenseLayer<Architecture>;

   // Random net.
   Net_t net(tbatchSize, timeSteps, tbatchSize, inputSize, 0, 0, 0, ELossFunction::kMeanSquaredError,
             EInitialization::kGauss);
   // FCLayer_t* l1 = net.AddDenseLayer(outputSize, EActivationFunction::kIdentity);
   net.AddDenseLayer(outputSize, EActivationFunction::kIdentity);
    
    auto & layers = net.GetLayers();
    auto bnlayer = new TBatchNormLayer<Architecture>(tbatchSize, outputSize);
    layers.push_back( bnlayer);
    //net.AddBatchNormLayer()


   // Random training data.
   std::vector<Matrix_t> X(timeSteps, Matrix_t(tbatchSize, inputSize)); // T x B x D
   Matrix_t Y(tbatchSize, outputSize), weights(tbatchSize, 1);
   net.Initialize();
   randomBatch(X[0]);
   randomMatrix(Y);
   fillMatrix(weights, 1.0);

   net.Forward(X);
   net.Backward(X, Y, weights);

   Scalar_t maximum_error = 0.0;

   // Compute derivatives for all weights using finite differences and
   // compare to result obtained from backpropagation.
   for (size_t l = 0; l < net.GetDepth(); l++) {
      std::cout << "\rTesting weight gradients:      layer: " << l << " / " << net.GetDepth();
      std::cout << std::flush;
      auto layer = net.GetLayerAt(l);
      auto &W = layer->GetWeightGradientsAt(0);

      int i = 0; 
      for (size_t j = 0; j < layer->GetInputWidth(); j++) {
         auto f = [&net, &X, &Y, &weights, l, i, j](Scalar_t x) {
            return evaluate_net_weight(net, X, Y, weights, l, 0, i, j, x);
         };
         Scalar_t dy = finiteDifference(f, dx) / (2.0 * dx);
         Scalar_t dy_ref = W(0, j);

         // Compute the relative error if dy != 0.
         Scalar_t error;
         if (std::fabs(dy_ref) > 1e-15) {
            error = std::fabs((dy - dy_ref) / dy_ref);
         } else {
            error = std::fabs(dy - dy_ref);
         }
         
         maximum_error = std::max(error, maximum_error);
      }
   }

   std::cout << "\rTesting weight gradients:      ";
   std::cout << "maximum relative error: " << print_error(maximum_error) << std::endl;
   return maximum_error;
}

