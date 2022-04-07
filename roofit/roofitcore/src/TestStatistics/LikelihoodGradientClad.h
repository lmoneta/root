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

#ifndef ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientClad
#define ROOT_ROOFIT_TESTSTATISTICS_LikelihoodGradientClad

#include "RooFit/MultiProcess/Job.h"
#include "RooFit/TestStatistics/LikelihoodGradientWrapper.h"

#include "Math/MinimizerOptions.h"
#include "Minuit2/NumericalDerivator.h"
#include "Minuit2/MnMatrix.h"

#include <vector>

namespace RooFit {
namespace TestStatistics {

class LikelihoodGradientClad : public LikelihoodGradientWrapper {
public:
   LikelihoodGradientClad(std::shared_ptr<RooAbsL> likelihood,
                         std::shared_ptr<WrapperCalculationCleanFlags> calculation_is_clean, std::size_t N_dim,
                         RooMinimizer *minimizer);
   LikelihoodGradientClad *clone() const override;
   LikelihoodGradientClad(const LikelihoodGradientClad &other);

   void fillGradient(double *grad) override;
   void fillGradientWithPrevResult(double *grad, double *previous_grad, double *previous_g2,
                                   double *previous_gstep) override;

   enum class GradientCalculatorMode { ExactlyMinuit2, AlmostMinuit2 };

private:
   void synchronizeParameterSettings(ROOT::Math::IMultiGenFunction *function,
                                     const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings) override;
   // this overload must also be overridden here so that the one above doesn't trigger a overloaded-virtual warning:
   void synchronizeParameterSettings(const std::vector<ROOT::Fit::ParameterSettings> &parameter_settings) override;

   void synchronizeWithMinimizer(const ROOT::Math::MinimizerOptions &options) override {}

   void updateMinuitInternalParameterValues(const std::vector<double> &minuit_internal_x) override;

   bool usesMinuitInternalValues() override { return false; }

   // members

   std::vector<double> minuit_internal_x_;
   std::vector<double> _grad;
};

} // namespace TestStatistics
} // namespace RooFit

#endif // ROOT_ROOFIT_LikelihoodGradientClad
