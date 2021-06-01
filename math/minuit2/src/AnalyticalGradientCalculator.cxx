// @(#)root/minuit2:$Id$
// Authors: M. Winkler, F. James, L. Moneta, A. Zsenei, E.G.P. Bos   2003-2017

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2005 LCG ROOT Math team,  CERN/PH-SFT                *
 * Copyright (c) 2017 Patrick Bos, Netherlands eScience Center        *
 *                                                                    *
 **********************************************************************/

#include "Minuit2/AnalyticalGradientCalculator.h"
#include "Minuit2/FCNGradientBase.h"
#include "Minuit2/MnUserTransformation.h"
#include "Minuit2/FunctionGradient.h"
#include "Minuit2/MinimumParameters.h"
#include "Minuit2/MnMatrix.h"
#include "Minuit2/MnPrint.h"

namespace ROOT {
namespace Minuit2 {


    FunctionGradient AnalyticalGradientCalculator::operator()(const MinimumParameters& par) const {
      // evaluate analytical gradient. take care of parameter transformations

      std::vector<double> grad = fGradCalc.Gradient(fTransformation(par.Vec()));
      assert(grad.size() == fTransformation.Parameters().size());

      MnAlgebraicVector v(par.Vec().size());
      for(unsigned int i = 0; i < par.Vec().size(); i++) {
        unsigned int ext = fTransformation.ExtOfInt(i);
        if(fTransformation.Parameter(ext).HasLimits()) {
          //double dd = (fTransformation.Parameter(ext).Upper() - fTransformation.Parameter(ext).Lower())*0.5*cos(par.Vec()(i));
          //       const ParameterTransformation * pt = fTransformation.transformation(ext);
          //       double dd = pt->dInt2ext(par.Vec()(i), fTransformation.Parameter(ext).Lower(), fTransformation.Parameter(ext).Upper() );
          double dd = fTransformation.DInt2Ext(i, par.Vec()(i));
          v(i) = dd*grad[ext];
        } else {
          v(i) = grad[ext];
        }
      }

      MnPrint print("AnalyticalGradientCalculator");
      print.Debug("User given gradient in Minuit2", v);

      // check for 2nd derivative and step-size from the external gradient
      // function and use them if present
      // N.B.: for the time being we only allow both at the same time, since
      //       FunctionGradient only has ctors for two cases: 1. gradient only,
      //       2. grad, g2 & gstep.
      if (fGradCalc.hasG2ndDerivative() && fGradCalc.hasGStepSize()) {
        std::vector<double> g2 = fGradCalc.G2ndDerivative(fTransformation(par.Vec()));
        std::vector<double> gstep = fGradCalc.GStepSize(fTransformation(par.Vec()));

        MnAlgebraicVector vg2(par.Vec().size());
        MnAlgebraicVector vgstep(par.Vec().size());
        for(unsigned int i = 0; i < par.Vec().size(); i++) {
          unsigned int ext = fTransformation.ExtOfInt(i);
          if(fTransformation.Parameter(ext).HasLimits()) {
            double int2ext_g2 = fTransformation.D2Int2Ext(i, par.Vec()(i));
            double int2ext_gstep = fTransformation.GStepInt2Ext(i, par.Vec()(i));
            vg2(i) = int2ext_g2 * g2[ext];
            vgstep(i) = int2ext_gstep * gstep[ext];
          } else {
            vg2(i) = g2[ext];
            vgstep(i) = gstep[ext];
          }
        }

        return FunctionGradient(v, vg2, vgstep);
      } else {
        return FunctionGradient(v);
      }
    }

    FunctionGradient AnalyticalGradientCalculator::operator()(const MinimumParameters& par, const FunctionGradient&) const {
      // needed from base class
      return (*this)(par);
    }

    bool AnalyticalGradientCalculator::CheckGradient() const {
      // check to be sure FCN implements analytical gradient
      return fGradCalc.CheckGradient();
    }

  }  // namespace Minuit2
}  // namespace ROOT
