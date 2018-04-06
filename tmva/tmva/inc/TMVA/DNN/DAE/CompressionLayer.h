// @(#)root/tmva/tmva/dnn:$Id$
// Author: Akshay Vashistha

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : TCompressionLayer                                                      *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Compressed Layer for DeepAutoEncoders                                     *
 *                                                                                *
 * Authors (alphabetical):                                                        *
 *      Akshay Vashistha <akshayvashistha1995@gmail.com>  - CERN, Switzerland     *
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


#ifndef TMVA_DAE_COMPRESSION_LAYER
#define TMVA_DAE_COMPRESSION_LAYER

#include "TMatrix.h"

#include "TMVA/DNN/GeneralLayer.h"
#include "TMVA/DNN/Functions.h"

#include <iostream>
#include <vector>

namespace TMVA {
namespace DNN {
namespace DAE {

  /** \class TCompressionLayer
    *  Used to Compress the input values.
    *  Here input is stored into compressed number of hidden units.
    *  These are basically the features after every layer.
  */

template <typename Architecture_t>
class TCompressionLayer : public VGeneralLayer<Architecture_t> {
public:
   using Matrix_t = typename Architecture_t::Matrix_t;
   using Scalar_t = typename Architecture_t::Scalar_t;


   size_t fVisibleUnits; ///< total number of visible units

   size_t fHiddenUnits; ///< tital number of hidden units.

   Scalar_t fDropoutProbability; ///< Probability that an input is active.

   size_t fType; ///< Type of layer

   EActivationFunction fF; ///< Activation function of the layer.

   /*! Constructor. */
   TCompressionLayer(size_t BatchSize, size_t VisibleUnits, size_t HiddenUnits,
                     Scalar_t DropoutProbability, EActivationFunction f,
                     std::vector<Matrix_t> Weights, std::vector<Matrix_t> Biases);

   /*! Copy the compression layer provided as a pointer */
   TCompressionLayer(TCompressionLayer<Architecture_t> *layer);

   /*! Copy constructor. */
   TCompressionLayer(const TCompressionLayer &);

   /* Destructor */
   ~TCompressionLayer();

   /* This forward pass compresses the input. Updated weights and biases are used from previous layers. */
   void Forward(std::vector<Matrix_t> &input, bool applyDropout = false);

   /* Not required in this layer */
   void Backward(std::vector<Matrix_t> &gradients_backward,
                 const std::vector<Matrix_t> &activations_backward,
                 std::vector<Matrix_t> &inp1,
                 std::vector<Matrix_t> &inp2);

   void Print() const;

   /*! Writes the information and the weights about the layer in an XML node. */
   virtual void AddWeightsXMLTo(void * /* parent */ ) {}

   /*! Read the information and the weights about the layer from XML node. */
   virtual void ReadWeightsFromXML(void * /* parent */ ) {}


   size_t GetVisibleUnits() const { return fVisibleUnits; }
   size_t GetHiddenUnits() const {return fHiddenUnits;}
   size_t GetType() const {return fType;}
   Scalar_t GetDropoutProbability() const { return fDropoutProbability; }
   EActivationFunction GetActivationFunction() const { return fF; }

};

//______________________________________________________________________________
template <typename Architecture_t>
TCompressionLayer<Architecture_t>::TCompressionLayer(size_t batchSize, size_t visibleUnits, size_t hiddenUnits,
                                                     Scalar_t dropoutProbability, EActivationFunction f,
                                                     std::vector<Matrix_t> weights, std::vector<Matrix_t> biases)
   : VGeneralLayer<Architecture_t>(batchSize, 1, 1, 0, 0, 0, 0, 1, {hiddenUnits}, {visibleUnits}, 2,
                                   {hiddenUnits, visibleUnits}, {1, 1}, batchSize, hiddenUnits, 1,
                                   EInitialization::kZero),
     fVisibleUnits(visibleUnits), fHiddenUnits(hiddenUnits), fDropoutProbability(dropoutProbability), fType(2), fF(f)

{
   Architecture_t::Copy(this->GetWeightsAt(0),weights[0]);
   Architecture_t::Copy(this->GetBiasesAt(0),biases[0]);
 }
//______________________________________________________________________________
 template <typename Architecture_t>
 TCompressionLayer<Architecture_t>::TCompressionLayer(TCompressionLayer<Architecture_t> *layer)
    : VGeneralLayer<Architecture_t>(layer), fVisibleUnits(layer->GetVisibleUnits()),
      fHiddenUnits(layer->GetHiddenUnits()), fDropoutProbability(layer->GetDropoutProbability()), fType(2),
      fF(layer->GetActivationFunction())
 {
    Architecture_t::Copy(this->GetWeightsAt(0), layer->weights[0]);
    Architecture_t::Copy(this->GetBiasesAt(0), layer->biases[0]);
    // Output Tensor will be created in General Layer
}
//______________________________________________________________________________
template <typename Architecture_t>
TCompressionLayer<Architecture_t>::TCompressionLayer(const TCompressionLayer &compress)
   : VGeneralLayer<Architecture_t>(compress), fVisibleUnits(compress.fVisibleUnits),
     fHiddenUnits(compress.fHiddenUnits), fDropoutProbability(compress.fDropoutProbability), fType(2), fF(compress.fF)

{
   Architecture_t::Copy(this->GetWeightsAt(0), compress.weights[0]);
   Architecture_t::Copy(this->GetBiasesAt(0), compress.biases[0]);
   // Output Tensor will be created in General Layer
}
//______________________________________________________________________________

template <typename Architecture_t> TCompressionLayer<Architecture_t>::~TCompressionLayer() {}

//______________________________________________________________________________
template <typename Architecture_t>
auto TCompressionLayer<Architecture_t>::Forward(std::vector<Matrix_t> &input, bool /*applyDropout*/) -> void
{

   for (size_t i = 0; i < this->GetBatchSize(); i++) {

      Architecture_t::EncodeInput(input[i], this->GetOutputAt(i), this->GetWeightsAt(0));
      Architecture_t::AddBiases(this->GetOutputAt(i), this->GetBiasesAt(0));
      evaluate<Architecture_t>(this->GetOutputAt(i), fF);
   }

}
//______________________________________________________________________________
template <typename Architecture_t>
auto inline TCompressionLayer<Architecture_t>::Backward(std::vector<Matrix_t> & /*gradients_backward*/,
                                                        const std::vector<Matrix_t> & /*activations_backward*/,
                                                        std::vector<Matrix_t> & /*inp1*/, std::vector<Matrix_t> &
                                                        /*inp2*/) -> void

{
}
//______________________________________________________________________________
template<typename Architecture_t>
auto TCompressionLayer<Architecture_t>::Print() const
-> void
{
   std::cout << "Batch Size: " << this->GetBatchSize() << "\n"
               << "Input Units: " << this->GetVisibleUnits() << "\n"
               << "Hidden units " << this->GetHiddenUnits() << "\n";

      std::cout<<"Compressed Input: "<<std::endl;
      for (Int_t j = 0; j < this->GetWeightsAt(0).GetNrows(); j++) {
         for (Int_t k = 0; k < this->GetWeightsAt(0).GetNcols(); k++) {
            std::cout<<this->GetWeightsAt(0)(j,k)<<"\t";
         }
         std::cout<<std::endl;
      }

}
//______________________________________________________________________________

}// namespace DAE
}//namespace DNN
}//namespace TMVA
#endif /* TMVA_DAE_Compression_LAYER*/
