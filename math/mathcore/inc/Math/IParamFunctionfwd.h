// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Nov 14 14:38:52 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Forward declarations for template class  IParamFunction class

#ifndef ROOT_Math_IParamFunctionfwd
#define ROOT_Math_IParamFunctionfwd

#ifndef ROOT_Math_IFunctionfwd
#include "Math/IFunctionfwd.h"
#endif

namespace ROOT {

   namespace Math {

      class IParametricFunctionOneDim;
      class IParametricGradFunctionOneDim;
      class IParametricFunctionMultiDim;
      class IParametricGradFunctionMultiDim;
      template<class T>
      class IParametricFunctionMultiDimTempl;
      template<class T>
      class IParametricGradFunctionMultiDimTempl;
      
      
      typedef IParametricFunctionOneDim        IParamFunction;
      typedef IParametricFunctionMultiDim      IParamMultiFunction;
      template<class T>
      using IParamMultiFunctionTempl = IParametricFunctionMultiDimTempl<T>;

      typedef IParametricGradFunctionOneDim        IParamGradFunction;
      typedef IParametricGradFunctionMultiDim      IParamMultiGradFunction;
      template<class T>
      using IParamMultiGradFunctionTempl =IParametricGradFunctionMultiDimTempl<T>;


   } // end namespace Math

} // end namespace ROOT


#endif /* ROOT_Math_IParamFunctionfwd */
