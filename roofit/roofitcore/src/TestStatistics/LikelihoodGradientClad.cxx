/*
 * Project: RooFit
 * Authors:
 *   PB, Patrick Bos, Netherlands eScience Center, p.bos@esciencecenter.nl
 *
 * Copyright (c) 2021, CERN
 *
 * Redistribution and use in source and binary forms,
 * with or without modification, are permitted according to the terms
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)
 */

#include "LikelihoodGradientClad.h"

#include "RooFit/MultiProcess/JobManager.h"
#include "RooFit/MultiProcess/Messenger.h"
#include "RooFit/MultiProcess/Queue.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"

#include "Minuit2/MnStrategy.h"

namespace RooFit {
namespace TestStatistics {

LikelihoodGradientClad::LikelihoodGradientClad(std::shared_ptr<RooAbsL> likelihood,
                                             std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean,
                                             std::size_t N_dim, RooMinimizer *minimizer)
   : LikelihoodGradientWrapper(std::move(likelihood), std::move(calculation_is_clean), N_dim, minimizer)
{
   // Note to future maintainers: take care when storing the minimizer_fcn pointer. The
   // RooAbsMinimizerFcn subclasses may get cloned inside MINUIT, which means the pointer
   // should also somehow be updated in this class.
   minuit_internal_x_.reserve(N_dim);

   _grad.resize(N_dim);
}

LikelihoodGradientClad::LikelihoodGradientClad(const LikelihoodGradientClad &other)
   : LikelihoodGradientWrapper(other),
     minuit_internal_x_(other.minuit_internal_x_)
{
   _grad.resize(minimizer_->getNPar());
}

LikelihoodGradientClad *LikelihoodGradientClad::clone() const
{
   return new LikelihoodGradientClad(*this);
}

void LikelihoodGradientClad::synchronizeParameterSettings(
   const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings)
{
   LikelihoodGradientWrapper::synchronizeParameterSettings(parameter_settings);
}

void LikelihoodGradientClad::synchronizeParameterSettings(
   ROOT::Math::IMultiGenFunction *function, const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings)
{
}

void LikelihoodGradientClad::fillGradient(double *grad)
{
   if (!calculation_is_clean_->gradient) {
   likelihood_->evaluateGradient(_grad.data());
   }

   for(std::size_t i = 0; i < minimizer_->getNPar(); ++i) grad[i] = _grad[i];
}

void LikelihoodGradientClad::fillGradientWithPrevResult(double *grad, double *previous_grad, double *previous_g2,
                                                       double *previous_gstep)
{
   if (!calculation_is_clean_->gradient) {
   likelihood_->evaluateGradient(_grad.data());
   }

   for(std::size_t i = 0; i < minimizer_->getNPar(); ++i) grad[i] = _grad[i];
}

void LikelihoodGradientClad::updateMinuitInternalParameterValues(const std::vector<double> &minuit_internal_x)
{
   minuit_internal_x_ = minuit_internal_x;
}

} // namespace TestStatistics
} // namespace RooFit
