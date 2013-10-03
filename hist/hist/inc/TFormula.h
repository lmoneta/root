// @(#)root/hist:$Id$
// Author: Maciej Zimnoch   30/09/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/
// ---------------------------------- TFormula.h
#ifndef ROOT_TFormula
#define ROOT_TFormula


#ifndef ROOT_TNamed
#include "TNamed.h"
#endif
 #ifndef ROOT_TBits
#include "TBits.h"
#endif
#ifndef ROOT_TObjArray
#include "TObjArray.h"
#endif
#include "TMethodCall.h"
#include <vector>
#include <list>
#include <map>

using namespace std;

class TFormulaFunction
{
public:
   TString  fName;
   TString  fBody;
   Int_t    fNargs;
   Bool_t   fFound;
   Bool_t   fFuncCall;
   TString  GetName()    { return fName; }
   TString  GetBody()    { return fBody; }
   Int_t    GetNargs()   { return fNargs;}
   Bool_t   IsFuncCall() { return fFuncCall;}
   TFormulaFunction(){}
   TFormulaFunction(const TString &name, const TString &body, int numArgs)
      : fName(name),fBody(body),fNargs(numArgs),fFound(false),fFuncCall(true) {}
   TFormulaFunction(const TString& name)
      : fName(name),fNargs(0),fFound(false),fFuncCall(false){}
   Bool_t operator<(const TFormulaFunction &rhv) const
   {  
      return fName < rhv.fName && fBody < rhv.fBody;
   }
   Bool_t operator==(const TFormulaFunction &rhv) const
   {  
      return fName == rhv.fName && fBody == rhv.fBody && fNargs == rhv.fNargs;
   }
};
class TFormulaVariable
{
public:
   TString fName;
   Double_t fValue;
   Int_t fArrayPos;
   Bool_t fFound;
   TString  GetName()      { return fName; }
   Double_t GetValue()     { return fValue; }
   Int_t    GetArrayPos()  { return fArrayPos; }
   TFormulaVariable():fName(""),fValue(-1),fArrayPos(-1),fFound(false){}
   TFormulaVariable(const TString &name, Double_t value, Int_t pos)
   : fName(name), fValue(value), fArrayPos(pos),fFound(false) {}
   Bool_t operator<(const TFormulaVariable &rhv) const
   {
      return fName < rhv.fName;
   }
};

class TFormula : public TNamed
{
private:

   TString           fClingInput;
   vector<Double_t>  fClingVariables;
   vector<Double_t>  fClingParameters;
   Bool_t            fReadyToExecute;
   Bool_t            fClingInitialized;
   Bool_t            fAllParametersSetted;
   TMethodCall*      fMethod;
   TString           fNamePrefix;

   void     InputFormulaIntoCling();
   void     PrepareEvalMethod();
   void     FillDefaults();
   void     HandlePolN(TString &formula);
   void     HandleParametrizedFunctions(TString &formula);
   void     HandleExponentiation(TString &formula);
   void     HandleLinear(TString &formula);
   Bool_t   IsDefaultVariableName(const TString &name);
protected:
   
   list<TFormulaFunction>         fFuncs;
   map<TString,TFormulaVariable>  fVars;
   map<TString,TFormulaVariable>  fParams;
   map<TString,Double_t>          fConsts;
   map<TString,TString>           fFunctionsShortcuts;
   TString                        fFormula;
   Int_t                          fNdim;
   Int_t                          fNpar;
   Int_t                          fNumber;
   TObjArray                      fLinearParts;

   Bool_t IsOperator(const char c);
   Bool_t IsBracket(const char c);
   Bool_t IsFunctionNameChar(const char c);
   void   ExtractFunctors(TString &formula);
   void   PreProcessFormula(TString &formula);
   void   ProcessFormula(TString &formula);
   Bool_t PrepareFormula(TString &formula);


public:
   enum {
      kNotGlobal     = BIT(10),  // don't store in gROOT->GetListOfFunction
      kNormalized    = BIT(14),   // set to true if the TFormula (ex gausn) is normalized
      kLinear        = BIT(16)    //set to true if the TFormula is for linear fitting
   };
                  TFormula();
   virtual        ~TFormula();
   TFormula&      operator=(const TFormula &rhs);
                  TFormula(const TString &name, TString formula);
                  TFormula(const TFormula &formula);
                  TFormula(const char *name, Int_t nparams, Int_t ndims);

   void           AddVariable(const TString &name, Double_t value);
   void           AddVariables(const pair<TString,Double_t> *vars, const Int_t size);
   void           AddParameter(const TString &name, Double_t value);
   Double_t       GetVariable(const TString &name);
   void           Copy(TObject &f1) const;
   Double_t       GetParameter(const TString &name);
   Double_t       GetParameter(Int_t param);
   Double_t*      GetParameters() const;
   void           GetParameters(Double_t *params);
   void           SetParameter(const TString &name, Double_t value);
   void           SetParameter(Int_t param, Double_t value);
   void           SetParameters(const Double_t *params,Int_t size);
   void           SetParameters(const Double_t *params);
   void           SetParameters(const pair<TString,Double_t> *params, const Int_t size);
   void           SetParameters(Double_t p0,Double_t p1,Double_t p2=0,Double_t p3=0,Double_t p4=0,
                                     Double_t p5=0,Double_t p6=0,Double_t p7=0,Double_t p8=0,
                                     Double_t p9=0,Double_t p10=0); // *MENU*
   void           SetParName(Int_t ipar, const char *name);
   void           SetParNames(const char *name0="p0",const char *name1="p1",const char
                             *name2="p2",const char *name3="p3",const char
                             *name4="p4", const char *name5="p5",const char *name6="p6",const char *name7="p7",const char
                             *name8="p8",const char *name9="p9",const char *name10="p10"); // *MENU*
   void           SetVariable(const TString &name, Double_t value);
   void           SetVariables(const pair<TString,Double_t> *vars, const Int_t size);
   TString        GetExpFormula() const { return fFormula; }
   const char    *GetParName(Int_t ipar) const;
   const TObject *GetLinearPart(Int_t i);
   Int_t          GetNdim() const {return fNdim;}
   Int_t          GetNpar() const {return fNpar;}
   Bool_t         IsValid() const { return fReadyToExecute && fAllParametersSetted; }
   Bool_t         IsLinear() const { return TestBit(kLinear); } 
   Int_t          GetNumber() const { return fNumber; }
   Double_t       EvalPar(const Double_t *x, const Double_t *params=0); // Back compatibility
   Double_t       Eval();
   Double_t       Eval(Double_t x);
   Double_t       Eval(Double_t x, Double_t y);
   Double_t       Eval(Double_t x, Double_t y , Double_t z);
   Double_t       Eval(Double_t x, Double_t y , Double_t z , Double_t t );

   ClassDef(TFormula,1)
};
#endif