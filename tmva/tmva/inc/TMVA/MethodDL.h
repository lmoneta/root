// @(#)root/tmva/tmva/dnn:$Id$
// Author: Vladimir Ilievski

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : MethodDL                                                              *
 * Web    : http://tmva.sourceforge.net                                           *
 *                                                                                *
 * Description:                                                                   *
 *      Deep Neural Network Method                                                *
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

#ifndef ROOT_TMVA_MethodDL
#define ROOT_TMVA_MethodDL

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// MethodDL                                                             //
//                                                                      //
// Method class for all Deep Learning Networks                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TString.h"

#include "TMVA/MethodBase.h"
#include "TMVA/Types.h"

#include "TMVA/DNN/Functions.h"
#include "TMVA/DNN/DeepNet.h"

#include <vector>

namespace TMVA {

/*! All of the options that can be specified in the training string */
struct TTrainingSettings {
   size_t batchSize;
   size_t testInterval;
   size_t convergenceSteps;
   DNN::ERegularization regularization;
   Double_t learningRate;
   Double_t momentum;
   Double_t weightDecay;
   std::vector<Double_t> dropoutProbabilities;
   bool multithreading;
};

class MethodDL : public MethodBase {

private:
   // Key-Value vector type, contining the values for the training options
   using KeyValueVector_t = std::vector<std::map<TString, TString>>;

   /*! The option handling methods */
   void DeclareOptions();
   void ProcessOptions();

   void Init();

   /*! After calling the ProcesOptions(), all of the options are parsed,
    *  so using the parsed options, and given the architecture and the
    *  type of the layers, we build the Deep Network passed as
    *  a reference in the function. */
   template <typename Architecture_t, typename Layer_t>
   void CreateDeepNet(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                      std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets);

   template <typename Architecture_t, typename Layer_t>
   void ParseDenseLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                        std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseConvLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                       std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseMaxPoolLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                          std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString,
                          TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseReshapeLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                          std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString,
                          TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseRnnLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                      std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);

   template <typename Architecture_t, typename Layer_t>
   void ParseLstmLayer(DNN::TDeepNet<Architecture_t, Layer_t> &deepNet,
                       std::vector<DNN::TDeepNet<Architecture_t, Layer_t>> &nets, TString layerString, TString delim);

   DNN::EInitialization fWeightInitialization; ///< The initialization method
   DNN::EOutputFunction fOutputFunction;       ///< The output function for making the predictions
   DNN::ELossFunction fLossFunction;           ///< The loss function

   TString fLayoutString;               ///< The string defining the layout of the CNN
   TString fErrorStrategy;              ///< The string defining the error strategy for training
   TString fTrainingStrategyString;     ///< The string defining the training strategy
   TString fWeightInitializationString; ///< The string defining the weight initialization method
   TString fArchitectureString;         ///< The string defining the architecure: CPU or GPU
   bool fResume;

   KeyValueVector_t fSettings;                       ///< Map for the training strategy
   std::vector<TTrainingSettings> fTrainingSettings; ///< The vector defining each training strategy

   ClassDef(MethodDL, 0);

protected:
   // provide a help message
   void GetHelpMessage() const;

public:
   /*! Constructor */
   MethodDL(const TString &jobName, const TString &methodTitle, DataSetInfo &theData, const TString &theOption);

   /*! Constructor */
   MethodDL(DataSetInfo &theData, const TString &theWeightFile);

   /*! Virtual Destructor */
   virtual ~MethodDL();

   /*! Function for parsing the training settings, provided as a string
    *  in a key-value form.  */
   KeyValueVector_t ParseKeyValueString(TString parseString, TString blockDelim, TString tokenDelim);

   /*! Check the type of analysis the deep learning network can do */
   Bool_t HasAnalysisType(Types::EAnalysisType type, UInt_t numberClasses, UInt_t numberTargets);

   /*! Methods for training the deep learning network */
   void Train();
   void TrainGpu();
   void TrainCpu();

   Double_t GetMvaValue(Double_t *err = 0, Double_t *errUpper = 0);

   /*! Methods for writing and reading weights */
   using MethodBase::ReadWeightsFromStream;
   void AddWeightsXMLTo(void *parent) const;
   void ReadWeightsFromXML(void *wghtnode);
   void ReadWeightsFromStream(std::istream &);

   /* Create ranking */
   const Ranking *CreateRanking();

   /* Getters */

   DNN::EInitialization GetWeightInitialization() { return fWeightInitialization; }
   DNN::EOutputFunction GetOutputFunction() { return fOutputFunction; }
   DNN::ELossFunction GetLossFunction() { return fLossFunction; }

   TString GetLayoutString() { return fLayoutString; }
   TString GetErrorStrategyString() { return fErrorStrategy; }
   TString GetTrainingStrategyString() { return fTrainingStrategyString; }
   TString GetWeightInitializationString() { return fWeightInitializationString; }
   TString GetArchitectureString() { return fArchitectureString; }

   const std::vector<TTrainingSettings> &GetTrainingSettings() const { return fTrainingSettings; }
   std::vector<TTrainingSettings> &GetTrainingSettings() { return fTrainingSettings; }
   const KeyValueVector_t &GetKeyValueSettings() const { return fSettings; }
   KeyValueVector_t &GetKeyValueSettings() { return fSettings; }

   /** Setters */
   void SetWeightInitialization(DNN::EInitialization weightInitialization)
   {
      fWeightInitialization = weightInitialization;
   }
   void SetOutputFunction(DNN::EOutputFunction outputFunction) { fOutputFunction = outputFunction; }
   void SetErrorStrategyString(TString errorStrategy) { fErrorStrategy = errorStrategy; }
   void SetTrainingStrategyString(TString trainingStrategyString) { fTrainingStrategyString = trainingStrategyString; }
   void SetWeightInitializationString(TString weightInitializationString)
   {
      fWeightInitializationString = weightInitializationString;
   }
   void SetArchitectureString(TString architectureString) { fArchitectureString = architectureString; }
   void SetLayoutString(TString layoutString) { fLayoutString = layoutString; }
};

} // namespace TMVA

#endif