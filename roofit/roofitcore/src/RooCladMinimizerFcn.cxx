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

//////////////////////////////////////////////////////////////////////////////
/// \class RooCladMinimizerFcn
/// RooCladMinimizerFcn is an interface to the ROOT::Math::IBaseFunctionMultiDim,
/// a function that ROOT's minimisers use to carry out minimisations.
///

#include "RooCladMinimizerFcn.h"

#include "RooAbsArg.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooNaNPacker.h"

#include "TClass.h"
#include "TMatrixDSym.h"

#include <fstream>
#include <iomanip>

using namespace std;


namespace {

// Helper function that wraps RooAbsArg::getParameters and directly returns the
// output RooArgSet. To be used in the initializer list of the RooCladMinimizerFcn
// constructor.
RooArgSet getParameters(RooAbsReal const& funct) {
    RooArgSet out;
    funct.getParameters(nullptr, out);
    return out;
}

} // namespace


RooCladMinimizerFcn::RooCladMinimizerFcn(RooAbsReal *funct, RooMinimizer* context,
            bool verbose) :
  RooAbsMinimizerFcn(getParameters(*funct), context, verbose), _funct(funct)
{}



RooCladMinimizerFcn::RooCladMinimizerFcn(const RooCladMinimizerFcn& other) : ROOT::Math::IMultiGradFunction(other), RooAbsMinimizerFcn(other),
  _funct(other._funct)
{}


RooCladMinimizerFcn::~RooCladMinimizerFcn()
{}


ROOT::Math::IBaseFunctionMultiDim* RooCladMinimizerFcn::Clone() const
{
  return new RooCladMinimizerFcn(*this) ;
}

void RooCladMinimizerFcn::setOptimizeConstOnFunction(RooAbsArg::ConstOpCode opcode, Bool_t doAlsoTrackingOpt)
{
   _funct->constOptimizeTestStatistic(opcode, doAlsoTrackingOpt);
}

/// Evaluate function given the parameters in `x`.
double RooCladMinimizerFcn::DoEval(const double *x) const {

  // Set the parameter values for this iteration
  for (unsigned index = 0; index < _nDim; index++) {
    if (_logfile) (*_logfile) << x[index] << " " ;
    SetPdfParamVal(index,x[index]);
  }

  // Calculate the function for these parameters
  RooAbsReal::setHideOffset(kFALSE) ;
  double fvalue = _funct->getVal();
  RooAbsReal::setHideOffset(kTRUE) ;

  if (!std::isfinite(fvalue) || RooAbsReal::numEvalErrors() > 0 || fvalue > 1e30) {
    printEvalErrors();
    RooAbsReal::clearEvalErrorLog() ;
    _numBadNLL++ ;

    if (_doEvalErrorWall) {
      const double badness = RooNaNPacker::unpackNaN(fvalue);
      fvalue = (std::isfinite(_maxFCN) ? _maxFCN : 0.) + _recoverFromNaNStrength * badness;
    }
  } else {
    if (_evalCounter > 0 && _evalCounter == _numBadNLL) {
      // This is the first time we get a valid function value; while before, the
      // function was always invalid. For invalid  cases, we returned values > 0.
      // Now, we offset valid values such that they are < 0.
      _funcOffset = -fvalue;
    }
    fvalue += _funcOffset;
    _maxFCN = std::max(fvalue, _maxFCN);
  }

  // Optional logging
  if (_logfile)
    (*_logfile) << setprecision(15) << fvalue << setprecision(4) << endl;
  if (_verbose) {
    cout << "\nprevFCN" << (_funct->isOffsetting()?"-offset":"") << " = " << setprecision(10)
         << fvalue << setprecision(4) << "  " ;
    cout.flush() ;
  }

  _evalCounter++ ;

  return fvalue;
}

std::string RooCladMinimizerFcn::getFunctionName() const
{
   return _funct->GetName();
}

std::string RooCladMinimizerFcn::getFunctionTitle() const
{
   return _funct->GetTitle();
}

void RooCladMinimizerFcn::setOffsetting(Bool_t flag)
{
   _funct->enableOffsetting(flag);
}

void RooCladMinimizerFcn::Gradient(const double *x, double *grad) const {
   DoEval(x);
   _funct->evaluateGradient(grad);
}

void RooCladMinimizerFcn::GradientWithPrevResult(const double *x, double *grad, double * /*previous_grad*/, double * /*previous_g2*/,
                            double * /*previous_gstep*/) const {
   Gradient(x, grad);
}

double RooCladMinimizerFcn::DoDerivative(const double *x, unsigned int icoord) const {
   double grad[_nDim];
   Gradient(x, grad);
   return grad[icoord];
}

double RooCladMinimizerFcn::DoDerivativeWithPrevResult(const double *x, unsigned int i_component, double *previous_grad,
                                  double *previous_g2, double *previous_gstep) const {
   double grad[_nDim];
   Gradient(x, grad);
   return grad[i_component];
}
