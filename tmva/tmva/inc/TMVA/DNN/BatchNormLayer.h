
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TBatchNormLayer                                                           *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Dense Layer Class                                                         *
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

#ifndef TMVA_DNN_BatchNormLayer
#define TMVA_DNN_BatchNormLayer

#include "TMatrix.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/Functions.h"

#include <iostream>
#include <iomanip>

namespace TMVA {
namespace DNN {
/** \class TBatchNormLayer

Generic layer class.

This generic layer class represents a dense layer of a neural network with
a given width n and activation function f. The activation function of each
layer is given by \f$\mathbf{u} = \mathbf{W}\mathbf{x} + \boldsymbol{\theta}\f$.

In addition to the weight and bias matrices, each layer allocates memory
for its activations and the corresponding first partial fDerivatives of
the activation function as well as the gradients of the weights and biases.

The layer provides member functions for the forward propagation of
activations through the given layer.
*/
template <typename Architecture_t>
class TBatchNormLayer : public VGeneralLayer<Architecture_t> {
public:
   using Scalar_t = typename Architecture_t::Scalar_t;
   //using Matrix_t = typename Architecture_t::Matrix_t;

    using Matrix_t = TMatrixD;
private:
   std::vector<Matrix_t> fDerivatives; ///< First fDerivatives of the activations of this layer.

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

   EActivationFunction fF; ///< Activation function of the layer.
   ERegularization fReg;   ///< The regularization method.
   Scalar_t fWeightDecay;  ///< The weight decay.

    Matrix_t xmu;
    Matrix_t xhat;
    Matrix_t gammax;
    Matrix_t out;
    Matrix_t dgamma;
    Matrix_t dgammax;
    Matrix_t dx;
    Matrix_t dbeta;
    Matrix_t var;
    Matrix_t sqrtvar;
    Matrix_t ivar;
    double epsilon;
    double gamma;
    double beta;
    double mu_Training;
    double var_Training;
    
public:
   /*! Constructor */
   TBatchNormLayer(size_t BatchSize, size_t InputWidth, size_t Width, EInitialization init, Scalar_t DropoutProbability,
               EActivationFunction f, ERegularization reg, Scalar_t weightDecay);

   /*! Copy the dense layer provided as a pointer */
   TBatchNormLayer(TBatchNormLayer<Architecture_t> *layer);

   /*! Copy Constructor */
   TBatchNormLayer(const TBatchNormLayer &);

   /*! Destructor */
   ~TBatchNormLayer();

   /*! Compute activation of the layer for the given input. The input
    * must be in 3D tensor form with the different matrices corresponding to
    * different events in the batch. Computes activations as well as
    * the first partial derivative of the activation function at those
    * activations. */
   void Forward(std::vector<Matrix_t> &input, bool inTraining = true);

   /*! Compute weight, bias and activation gradients. Uses the precomputed
    *  first partial derviatives of the activation function computed during
    *  forward propagation and modifies them. Must only be called directly
    *  a the corresponding call to Forward(...). */
   void Backward(std::vector<Matrix_t> &gradients_backward, const std::vector<Matrix_t> &activations_backward);

   /*! Printing the layer info. */
   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void *parent);

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void *parent);

   const std::vector<Matrix_t> &GetDerivatives() const { return fDerivatives; }
   std::vector<Matrix_t> &GetDerivatives() { return fDerivatives; }

   Matrix_t &GetDerivativesAt(size_t i) { return fDerivatives[i]; }
   const Matrix_t &GetDerivativesAt(size_t i) const { return fDerivatives[i]; }

   EActivationFunction GetActivationFunction() const { return fF; }
   ERegularization GetRegularization() const { return fReg; }
   Scalar_t GetWeightDecay() const { return fWeightDecay; }
};

//
//
//  The Dense Layer Class - Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::TBatchNormLayer(size_t batchSize, size_t inputWidth, size_t width, EInitialization init,
                                         Scalar_t dropoutProbability, EActivationFunction f, ERegularization reg,
                                         Scalar_t weightDecay)
   : VGeneralLayer<Architecture_t>(batchSize, 1, 1, inputWidth, 1, 1, width, 1, width, inputWidth, 1, width, 1, 1,
                                   batchSize, width, init),
     fDerivatives(), fDropoutProbability(dropoutProbability), fF(f), fReg(reg), fWeightDecay(weightDecay)
{
   fDerivatives.emplace_back(batchSize, width);
}

//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::TBatchNormLayer(TBatchNormLayer<Architecture_t> *layer)
   : VGeneralLayer<Architecture_t>(layer), fDerivatives(), fDropoutProbability(layer->GetDropoutProbability()),
     fF(layer->GetActivationFunction()), fReg(layer->GetRegularization()), fWeightDecay(layer->GetWeightDecay())
{
   fDerivatives.emplace_back(layer->GetBatchSize(), layer->GetWidth());
}

//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::TBatchNormLayer(const TBatchNormLayer &layer)
   : VGeneralLayer<Architecture_t>(layer), fDerivatives(), fDropoutProbability(layer.fDropoutProbability), fF(layer.fF),
     fReg(layer.fReg), fWeightDecay(layer.fWeightDecay)
{
   fDerivatives.emplace_back(layer.fBatchSize, layer.fWidth);
}

//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::~TBatchNormLayer()
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TBatchNormLayer<Architecture_t>::Forward(std::vector<Matrix_t> &x, bool inTraining) -> void
{
    Matrix_t & input = x[0];
    
    Matrix_t & gamma = this->GetWeightsAt(0);
    Matrix_t & beta = this->GetWeightsAt(1);
    
    
    int n = input.GetNrows();
    int d = input.GetNcols();
    
    for (int k = 0; k < d; ++k) {
        double mean = 0;
        for ( int i = 0; i < n; i++ ) {
            mean = mean + input(i,k);
        }
        mean = mean / n;
        
        for( int i = 0; i < n; i++ ) {
            
            var(i,k)  = 0;
            xmu(i,k) = input(i,k) - mean;
        }
        double sq = 0;
        
        
        for ( int i = 0; i < n; i++ ) {
            
            sq = sq + (xmu(i,k) * xmu(i,k));
        }
        
        for( int i = 0; i < n; i++ ) {
            
            
            var(i,k) = sq / n;
            
            var(i,k) = var(i,k) + epsilon;
            
            sqrtvar(i,k) = std::sqrt(var(i,k));
            ivar(i,k) = 1./sqrtvar(i,k);
            xhat(i,k) = xmu(i,k) * ivar(i,k);
            double gammax = gamma(0,k) * xhat(i,k);
            out(i,k) = gammax + beta(0,k);
        }
    }
}
    
//______________________________________________________________________________
template <typename Architecture_t>
auto TBatchNormLayer<Architecture_t>::Backward(std::vector<Matrix_t> &gradients_backward, const std::vector<Matrix_t> &activations_backward) -> void
{
    const Matrix_t & dout = activations_backward[0];
    
    int d = dout.GetNcols();
    int n = dout.GetNrows();
    
    dgamma.ResizeTo(1,d);
    dx.ResizeTo(n,d);
    dbeta.ResizeTo(1,d);
    
    TMatrixD dxhat(n,d);
    
    for( int k = 0; k < d; k++) {
        for ( int i = 0; i < n; i++ ) {
            dgamma(0,k) += dout(i,k) * xhat(i,k);
            dxhat(i,k) = dout(i,k) * gamma;
        }
    }
    std::vector<double> divar(d);
    std::vector<double> dmu(d);
    std::vector<double> dx1(n);
    
    for( int k = 0; k < d; k++) {
        divar[k] = 0;
        for ( int i = 0; i < n; i++ ) {
            divar[k] += dxhat(i,k) * xmu(i,k);
        }
        for( int i = 0; i < n; i++ ) {
            double dxmu1 = dxhat(i,k) * ivar(i,k);
            double dsqrtvar = -1. /(pow(sqrtvar(i,k),2)) * divar[k];
            double dvar = 0.5 * 1. /sqrt(var(i,k) + epsilon) * dsqrtvar;
            double dsq = dvar / n;
            double dxmu2 = 2 * xmu(i,k) * dsq;
            dmu[k] += -1 * dxmu1 + dxmu2;
            dx1[i] = (dxmu1 + dxmu2);
        }
        double dx2 = dmu[k] / n;
        
        for( int i = 0; i < n; i++ ) {
            dx(i,k) = dx1[i] + dx2;
            
        }
        
    }
    
}

//______________________________________________________________________________
template <typename Architecture_t>
void TBatchNormLayer<Architecture_t>::Print() const
{
   std::cout << " DENSE Layer: \t";
   std::cout << " ( Input =" << std::setw(6) << this->GetWeightsAt(0).GetNcols();  // input size 
   std::cout << " , Width =" << std::setw(6) << this->GetWeightsAt(0).GetNrows() << " ) ";  // layer width
   if (this->GetOutput().size() > 0) {
      std::cout << "\tOutput = ( " << std::setw(2) << this->GetOutput().size() << " ," << std::setw(6) << this->GetOutput()[0].GetNrows() << " ," << std::setw(6) << this->GetOutput()[0].GetNcols() << " ) ";
   }
   std::vector<std::string> activationNames = { "Identity","Relu","Sigmoid","Tanh","SymmRelu","SoftSign","Gauss" };
   std::cout << "\t Activation Function = ";
   std::cout << activationNames[ static_cast<int>(fF) ];
   if (fDropoutProbability != 1.) std::cout << "\t Dropout prob. = " << fDropoutProbability;
   std::cout << std::endl;
}

//______________________________________________________________________________

template <typename Architecture_t>
void TBatchNormLayer<Architecture_t>::AddWeightsXMLTo(void *parent)
{
  // write layer width activation function + weigbht and bias matrices

   auto layerxml = gTools().xmlengine().NewChild(parent, 0, "BatchNormLayer");

   gTools().xmlengine().NewAttr(layerxml, 0, "Width", gTools().StringFromInt(this->GetWidth()));

   int activationFunction = static_cast<int>(this -> GetActivationFunction());
   gTools().xmlengine().NewAttr(layerxml, 0, "ActivationFunction",
                                TString::Itoa(activationFunction, 10));
   // write weights and bias matrix 
   this->WriteMatrixToXML(layerxml, "Weights", this -> GetWeightsAt(0));
   this->WriteMatrixToXML(layerxml, "Biases",  this -> GetBiasesAt(0));
}

//______________________________________________________________________________
template <typename Architecture_t>
void TBatchNormLayer<Architecture_t>::ReadWeightsFromXML(void *parent)
{
   // Read layer weights and biases from XML
   this->ReadMatrixXML(parent,"Weights", this -> GetWeightsAt(0));
   this->ReadMatrixXML(parent,"Biases", this -> GetBiasesAt(0));
   
}

} // namespace DNN
} // namespace TMVA

#endif
