// @(#)root/mathcore:$Id$
// Author: L. Moneta Tue Sep  5 09:13:32 2006

/**********************************************************************
 *                                                                    *
 * Copyright (c) 2006  LCG ROOT Math Team, CERN/PH-SFT                *
 *                                                                    *
 *                                                                    *
 **********************************************************************/

// Header file for class Chi2FCN

#ifndef ROOT_Fit_Chi2FCN
#define ROOT_Fit_Chi2FCN

#ifndef ROOT_Fit_BasicFCN
#include "Fit/BasicFCN.h"
#endif

#ifndef ROOT_Math_IParamFunction
#include "Math/IParamFunction.h"
#endif

#include "Math/IFunction.h"
#include "Math/IFunctionfwd.h"

#ifndef ROOT_Fit_BinData
#include "Fit/BinData.h"
#endif


#ifndef ROOT_Fit_FitUtil
#include "Fit/FitUtil.h"
#endif
//#define ROOT_FIT_PARALLEL

#ifdef ROOT_FIT_PARALLEL
#ifndef ROOT_Fit_FitUtilParallel
#include "Fit/FitUtilParallel.h"
#endif
#endif

#include <memory>
#include <VecCore/VecCore>

/**
@defgroup FitMethodFunc Fit Method Classes

Classes describing Fit Method functions
@ingroup Fit
*/

using Double_v = typename vecCore::backend::VcVector::Double_v; 


namespace ROOT {


   namespace Fit {

template<class T>
struct SumOfT{
  static double Sum(T &x){
    return x.sum();
  }
};

template<>
struct SumOfT<double>{
  static double Sum( double &x){return x;}
};

//___________________________________________________________________________________
/**
   Chi2FCN class for binnned fits using the least square methods

   @ingroup  FitMethodFunc
*/
template<class FunType>
class Chi2FCN : public BasicFCN<FunType, BinData> {

public:

   typedef typename FunType::BackendType T;
   typedef  BasicFCN<FunType, BinData> BaseFCN; 

   typedef  ::ROOT::Math::BasicFitMethodFunction<FunType> BaseObjFunction;
   typedef typename  BaseObjFunction::BaseFunction BaseFunction;

   //typedef  typename ::ROOT::Math::ParamFunctionTrait<FunType>::PFType IModelFunction;
   typedef  ::ROOT::Math::IParamMultiFunctionTempl<T> IModelFunction;
   typedef typename BaseObjFunction::Type_t Type_t;

   /**
      Constructor from data set (binned ) and model function
   */
   Chi2FCN (const std::shared_ptr<BinData> & data, const std::shared_ptr<IModelFunction> & func) :
      BaseFCN( data, func),
      fNEffPoints(0),
      fGrad ( std::vector<double> ( func->NPar() ) )
   { }

   /**
      Same Constructor from data set (binned ) and model function but now managed by the user
      we clone the function but not the data
   */
   Chi2FCN ( const BinData & data, const IModelFunction & func) :
      BaseFCN(std::shared_ptr<BinData>(const_cast<BinData*>(&data), DummyDeleter<BinData>()), std::shared_ptr<IModelFunction>(dynamic_cast<IModelFunction*>(func.Clone() ) ) ),
      fNEffPoints(0),
      fGrad ( std::vector<double> ( func.NPar() ) )
   { }

   /**
      Destructor (no operations)
   */
   virtual ~Chi2FCN ()  {}
   /**
      Copy constructor
   */
   Chi2FCN(const Chi2FCN & f) :
      BaseFCN(f.DataPtr(), f.ModelFunctionPtr() ),
      fNEffPoints( f.fNEffPoints ),
      fGrad( f.fGrad)
   {  }

   /**
      Assignment operator
   */
   Chi2FCN & operator = (const Chi2FCN & rhs) {
      SetData(rhs.DataPtr() );
      SetModelFunction(rhs.ModelFunctionPtr() );
      fNEffPoints = rhs.fNEffPoints;
      fGrad = rhs.fGrad; 
   }

   /* 
      clone the function
    */
   virtual BaseFunction * Clone() const {
      // return new Chi2FCN(*this);
      return nullptr; 
   }



   using BaseObjFunction::operator();


   /// i-th chi-square residual
   virtual double DataElement(const double * x, unsigned int i, double * g) const {
      if (i==0) this->UpdateNCalls();
      return 0.0;//FitUtil::EvaluateChi2Residual(BaseFCN::ModelFunction(), BaseFCN::Data(), x, i, g);
   }

   // need to be virtual to be instantiated
   virtual void Gradient(const double *x, double *g) const {
      // evaluate the chi2 gradient
      // FitUtil::EvaluateChi2Gradient(BaseFCN::ModelFunction(), BaseFCN::Data(), x, g, fNEffPoints);
   }

   /// get type of fit method function
   virtual  typename BaseObjFunction::Type_t Type() const { return BaseObjFunction::kLeastSquare; }



protected:

   /// set number of fit points (need to be called in const methods, make it const)                                                                                                      
   virtual void SetNFitPoints(unsigned int n) const { fNEffPoints = n; }
   
private:

   /**
      Evaluation of the  function (required by interface)
    */
   virtual double DoEval (const double * p) const {
//       this->UpdateNCalls();
// #ifdef ROOT_FIT_PARALLEL
//       return FitUtilParallel::EvaluateChi2(BaseFCN::ModelFunction(), BaseFCN::Data(), x, fNEffPoints);
// #else
//       if (!BaseFCN::Data().HaveCoordErrors() )
//          return FitUtil::EvaluateChi2(BaseFCN::ModelFunction(), BaseFCN::Data(), x, fNEffPoints);
//       else
//          return FitUtil::EvaluateChi2Effective(BaseFCN::ModelFunction(), BaseFCN::Data(), x, fNEffPoints);
// #endif
    T tmp{};
    for(int i=0; i<BaseFCN::Data().Size(); i+=4){
        auto ux = vecCore::FromPtr<T>(BaseFCN::Data().GetCoordComponent(i,0));
        auto &f = BaseFCN::ModelFunction(); 
        tmp += f(&ux, p);
    }
    return SumOfT<T>::Sum(tmp);
   }

   
   // for derivatives
   virtual double  DoDerivative(const double * x, unsigned int icoord ) const {
      Gradient(x, fGrad.data());
      return fGrad[icoord];
   }


   mutable unsigned int fNEffPoints;  // number of effective points used in the fit

   mutable std::vector<double> fGrad; // for derivatives


};

      // define useful typedef's
      typedef Chi2FCN<ROOT::Math::IMultiGenFunction> Chi2Function;
      typedef Chi2FCN<ROOT::Math::IMultiGradFunction> Chi2GradFunction;


   } // end namespace Fit

} // end namespace ROOT


#endif /* ROOT_Fit_Chi2FCN */
