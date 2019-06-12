
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

   
   
   
    Scalar_t fMomentum;  ///< The weight decay.
    Scalar_t fEpsilon;

    Matrix_t xmu;
    Matrix_t xhat;
    
    
    Matrix_t dgammax;
    Matrix_t var;
    Matrix_t sqrtvar;
    Matrix_t ivar;
    
    
    //std::vector<Scalar_t> fGamma;
    //std::vector<Scalar_t> fBeta;
    std::vector<Scalar_t> fMu_Training;
    std::vector<Scalar_t> fVar_Training;
    
public:
   /*! Constructor */
    TBatchNormLayer(size_t BatchSize, size_t InputWidth,  Scalar_t epsilon = 0.0001);

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
   void Backward(std::vector<Matrix_t> &gradients_backward, const std::vector<Matrix_t> &activations_backward, std::vector<Matrix_t> &inp1, std::vector<Matrix_t> &inp2);
    

   /*! Printing the layer info. */
   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
    virtual void AddWeightsXMLTo(void *parent) {}

   /*! Read the information and the weights about the layer from XML node. */
    virtual void ReadWeightsFromXML(void *parent) {}

   /* initialize weights */
   virtual void Initialize() {}
   
   //Scalar_t GetWeightDecay() const { return fWeightDecay; }
};

//
//
//  The Dense Layer Class - Implementation
//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::TBatchNormLayer(size_t batchSize, size_t inputWidth, Scalar_t epsilon)
   : VGeneralLayer<Architecture_t>(batchSize, 1, 1, inputWidth, 1, 1, inputWidth,
                                   2, 1, inputWidth,   // weight tensor dim.
                                   1, 1, 1,   // bias
                                   1 , batchSize, inputWidth,  // output tensor
                                   EInitialization::kZero),
   fEpsilon(epsilon),
    xmu(batchSize, inputWidth),
    xhat(batchSize, inputWidth),
    dgammax(batchSize, inputWidth),
    var(batchSize, inputWidth),
    sqrtvar(batchSize, inputWidth),
    ivar(batchSize, inputWidth),
    fMu_Training(inputWidth),
    fVar_Training(inputWidth)
    
    //
    
{
    fMu_Training.assign(inputWidth, 0.);
    fVar_Training.assign(inputWidth, 1.);
    Matrix_t & gamma = this->GetWeightsAt(0);
    for (int i = 0; i < gamma.GetNcols(); ++i)
        gamma(0,i) = 1;
    Matrix_t & beta = this->GetWeightsAt(1);
    for (int i = 0; i < gamma.GetNcols(); ++i)
        beta(0,i) = 0;

    printf("constructed gamma \n");
    this->GetWeightsAt(0).Print();
}

//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::TBatchNormLayer(TBatchNormLayer<Architecture_t> *layer)
   : VGeneralLayer<Architecture_t>(layer)
{
   // to be implemented
    printf("Error - copy ctor not implmented\n");
}

//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::TBatchNormLayer(const TBatchNormLayer &layer)
   : VGeneralLayer<Architecture_t>(layer)
{
   // to be implmeented
    printf("Error - copy ctor not implmented\n");
}

//______________________________________________________________________________
template <typename Architecture_t>
TBatchNormLayer<Architecture_t>::~TBatchNormLayer()
{
   // Nothing to do here.
}

//______________________________________________________________________________
template <typename Architecture_t>
auto TBatchNormLayer<Architecture_t>::Forward(std::vector<Matrix_t> &x, bool /* inTraining */) -> void
{
    Matrix_t & input = x[0];
    
    Matrix_t & gamma = this->GetWeightsAt(0);
    Matrix_t & beta = this->GetWeightsAt(1);

    //printf("gamma and beta \n");
    //gamma.Print();
    //beta.Print(); 
    
    Matrix_t & out = this->GetOutputAt(0);
    double epsilon = fEpsilon;
    
    int n = input.GetNrows();
    int d = input.GetNcols();
    
    for (int k = 0; k < d; ++k) {
        double mean = 0;
        for ( int i = 0; i < n; i++ ) {
            mean = mean + input(i,k);
        }
        mean = mean / n;
        
        for( int i = 0; i < n; i++ ) {
            
            var(0,k)  = 0;
            xmu(i,k) = input(i,k) - mean;
        }
        double sq = 0;
        
        
        for ( int i = 0; i < n; i++ ) {
            
            sq = sq + (xmu(i,k) * xmu(i,k));
        }
        var(0,k) = sq / n;
        var(0,k) = var(0,k) + epsilon;
        sqrtvar(0,k) = std::sqrt(var(0,k));
        ivar(0,k) = 1./sqrtvar(0,k);
        for( int i = 0; i < n; i++ ) {
            xhat(i,k) = xmu(i,k) * ivar(0,k);
            double gammax = gamma(0,k) * xhat(i,k);
            out(i,k) = gammax + beta(0,k);
        }
    }
    //std::cout << "xhat \n";
    //xhat.Print();
    // std::cout << "out \n";
    // out.Print();
    // std::cout << "beta \n";
    // beta.Print(); 
}
    
//______________________________________________________________________________
template <typename Architecture_t>
auto TBatchNormLayer<Architecture_t>::Backward(std::vector<Matrix_t> &gradients_backward, const std::vector<Matrix_t> &activations_backward,std::vector<Matrix_t> &, std::vector<Matrix_t> &) -> void
{
    const Matrix_t & dout = this->GetActivationGradients()[0];

       
    double epsilon = fEpsilon;
    
    int d = dout.GetNcols();
    int n = dout.GetNrows();
    
    const Matrix_t & gamma = this->GetWeightsAt(0);
    
    Matrix_t & dgamma = this->GetWeightGradientsAt(0);
    Matrix_t & dbeta = this->GetWeightGradientsAt(1);
    Matrix_t & dx = gradients_backward[0];
    Matrix_t & dh = gradients_backward[0];
    const Matrix_t & x = activations_backward[0];
    
    TMatrixD dxhat(n,d);
    
    printf("doing backpropagation \n");
    printf("input dout\n");
    dout.Print();  
     
    for( int k = 0; k < d; k++) {
        dgamma(0,k) = 0;
        dbeta(0,k) = 0;
        for ( int i = 0; i < n; i++ ) {
            dbeta(0,k) += dout(i,k);
            dgamma(0,k) += dout(i,k) * xhat(i,k);
            dxhat(i,k) = dout(i,k) * gamma(0,k);
        }
    }
    std::vector<double> divar(d);
    std::vector<double> dmu(d);

    printf("computed dbeta \n");
    dbeta.Print();
    printf("computed dgamma \n");
    dgamma.Print(); 
    
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
    printf("output gradient method1 \n");
    dx.Print();
    
    double npSumDy = 0;
    double npSumDyHMu = 0;
    
    for ( int k = 0; k < d; k++) {
        for ( int i = 0; i < n; i++) {
            npSumDy += dout(i,k);
            npSumDyHMu += dout(i,k) * xmu(i,k);
        }
        for ( int i = 0; i < n; i++) {
            dx(i,k) = (1./double(n) * gamma(0,k) * ivar(0,k)) * (n * dout(i,k) - npSumDy - xmu(i,k) * pow(var(0,k),-1.0) * npSumDyHMu);
        }
        
    }
    
    
    printf("output gradient method2 \n");
    dx.Print(); 
}

//______________________________________________________________________________
template <typename Architecture_t>
void TBatchNormLayer<Architecture_t>::Print() const
{
   std::cout << " BATCH NORM Layer: \t";
   std::cout << " ( Input =" << std::setw(6) << this->GetWeightsAt(0).GetNcols() << " ) ";
   std::cout << std::endl;
}

//______________________________________________________________________________
/*
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
 */

} // namespace DNN
} // namespace TMVA

#endif
