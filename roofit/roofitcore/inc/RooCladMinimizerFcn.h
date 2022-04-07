/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   AL, Alfio Lazzaro,   INFN Milan,        alfio.lazzaro@mi.infn.it        *
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl   *
 *                                                                           *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef roofit_roofitcore_RooCladMinimizerFcn_h
#define roofit_roofitcore_RooCladMinimizerFcn_h

#include "Math/IFunction.h"
#include "Fit/ParameterSettings.h"
#include "Fit/FitResult.h"

#include "RooAbsReal.h"
#include "RooArgList.h"

#include <fstream>
#include <vector>

#include <RooAbsMinimizerFcn.h>

template<typename T> class TMatrixTSym;
using TMatrixDSym = TMatrixTSym<double>;

// forward declaration
class RooMinimizer;

// class RooCladMinimizerFcn : public RooAbsMinimizerFcn, public ROOT::Math::IBaseFunctionMultiDim {
class RooCladMinimizerFcn : public ROOT::Math::IMultiGradFunction, public RooAbsMinimizerFcn {

public:
   RooCladMinimizerFcn(RooAbsReal *funct, RooMinimizer *context, bool verbose = false);
   RooCladMinimizerFcn(const RooCladMinimizerFcn &other);
   ~RooCladMinimizerFcn() override;

   ROOT::Math::IBaseFunctionMultiDim *Clone() const override;

   /// IMultiGradFunction overrides necessary for Minuit
   void Gradient(const double *x, double *grad) const override;
   void GradientWithPrevResult(const double *x, double *grad, double *previous_grad, double *previous_g2,
                               double *previous_gstep) const override;

   unsigned int NDim() const override { return getNDim(); }

   std::string getFunctionName() const override;
   std::string getFunctionTitle() const override;

   void setOptimizeConstOnFunction(RooAbsArg::ConstOpCode opcode, Bool_t doAlsoTrackingOpt) override;

   void setOffsetting(Bool_t flag) override;

private:
   /// IMultiGradFunction overrides necessary for Minuit
   double DoDerivative(const double *x, unsigned int icoord) const override;
   double DoDerivativeWithPrevResult(const double *x, unsigned int i_component, double *previous_grad,
                                     double *previous_g2, double *previous_gstep) const override;

   double DoEval(const double *x) const override;

   RooAbsReal *_funct;
};

#endif
