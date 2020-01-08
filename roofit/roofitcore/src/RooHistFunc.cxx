/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofit:$Id$
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (c) 2000-2005, Regents of the University of California          *
 *                          and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

/**
\file RooHistFunc.cxx
\class RooHistFunc
\ingroup Roofitcore

RooHistFunc implements a real-valued function sampled from a 
multidimensional histogram. The histogram can have an arbitrary number of real or 
discrete dimensions and may have negative values
**/

#include "RooFit.h"
#include "Riostream.h"

#include "RooHistFunc.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooWorkspace.h"

#include "TError.h"

using namespace std;

ClassImp(RooHistFunc)
;



////////////////////////////////////////////////////////////////////////////////
/// Default constructor

RooHistFunc::RooHistFunc() :
  _dataHist(0),
  _intOrder(0),
  _cdfBoundaries(kFALSE),
  _totVolume(0),
  _unitNorm(kFALSE)
{
  TRACE_CREATE 

  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _depList.createIterator() ;
}


////////////////////////////////////////////////////////////////////////////////
/// Constructor from a RooDataHist. The variable listed in 'vars' control the dimensionality of the
/// function. Any additional dimensions present in 'dhist' will be projected out. RooDataHist dimensions
/// can be either real or discrete. See RooDataHist::RooDataHist for details on the binning.
/// RooHistFunc neither owns or clone 'dhist' and the user must ensure the input histogram exists
/// for the entire life span of this function.

RooHistFunc::RooHistFunc(const char *name, const char *title, const RooArgSet& vars, 
		       const RooDataHist& dhist, Int_t intOrder) :
  RooAbsReal(name,title), 
  _depList("depList","List of dependents",this),
  _dataHist((RooDataHist*)&dhist), 
  _codeReg(10),
  _intOrder(intOrder),
  _cdfBoundaries(kFALSE),
  _totVolume(0),
  _unitNorm(kFALSE)
{
  _histObsList.addClone(vars) ;
  _depList.add(vars) ;

  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _depList.createIterator() ;

  // Verify that vars and dhist.get() have identical contents
  const RooArgSet* dvars = dhist.get() ;
  if (vars.getSize()!=dvars->getSize()) {
    coutE(InputArguments) << "RooHistFunc::ctor(" << GetName() 
			  << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
    assert(0) ;
  }
  TIterator* iter = vars.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!dvars->find(arg->GetName())) {
      coutE(InputArguments) << "RooHistFunc::ctor(" << GetName() 
			    << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
      assert(0) ;
    }
  }
  delete iter ;
  TRACE_CREATE 
}



////////////////////////////////////////////////////////////////////////////////
/// Constructor from a RooDataHist. The variable listed in 'vars' control the dimensionality of the
/// function. Any additional dimensions present in 'dhist' will be projected out. RooDataHist dimensions
/// can be either real or discrete. See RooDataHist::RooDataHist for details on the binning.
/// RooHistFunc neither owns or clone 'dhist' and the user must ensure the input histogram exists
/// for the entire life span of this function.

RooHistFunc::RooHistFunc(const char *name, const char *title, const RooArgList& funcObs, const RooArgList& histObs, 
  		       const RooDataHist& dhist, Int_t intOrder) :
  RooAbsReal(name,title), 
  _depList("depList","List of dependents",this),
  _dataHist((RooDataHist*)&dhist), 
  _codeReg(10),
  _intOrder(intOrder),
  _cdfBoundaries(kFALSE),
  _totVolume(0),
  _unitNorm(kFALSE)
{
  _histObsList.addClone(histObs) ;
  _depList.add(funcObs) ;

  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _depList.createIterator() ;

  // Verify that vars and dhist.get() have identical contents
  const RooArgSet* dvars = dhist.get() ;
  if (histObs.getSize()!=dvars->getSize()) {
    coutE(InputArguments) << "RooHistFunc::ctor(" << GetName() 
			  << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
    assert(0) ;
  }
  TIterator* iter = histObs.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    if (!dvars->find(arg->GetName())) {
      coutE(InputArguments) << "RooHistFunc::ctor(" << GetName() 
			    << ") ERROR variable list and RooDataHist must contain the same variables." << endl ;
      assert(0) ;
    }
  }
  delete iter ;
  TRACE_CREATE 
}



////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

RooHistFunc::RooHistFunc(const RooHistFunc& other, const char* name) :
  RooAbsReal(other,name), 
  _depList("depList",this,other._depList),
  _dataHist(other._dataHist),
  _codeReg(other._codeReg),
  _intOrder(other._intOrder),
  _cdfBoundaries(other._cdfBoundaries),
  _totVolume(other._totVolume),
  _unitNorm(other._unitNorm)
{
  TRACE_CREATE 

  _histObsList.addClone(other._histObsList) ;

  _histObsIter = _histObsList.createIterator() ;
  _pdfObsIter = _depList.createIterator() ;
}



////////////////////////////////////////////////////////////////////////////////

RooHistFunc::~RooHistFunc() 
{ 
  TRACE_DESTROY 

  delete _histObsIter ;
  delete _pdfObsIter ;
}




////////////////////////////////////////////////////////////////////////////////
/// Return the current value: The value of the bin enclosing the current coordinates
/// of the dependents, normalized by the histograms contents. Interpolation
/// is applied if the RooHistFunc is configured to do that

Double_t RooHistFunc::evaluate() const
{
  // Transfer values from   
  if (_depList.getSize()>0) {
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    RooAbsArg* harg, *parg ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
      parg = (RooAbsArg*)_pdfObsIter->Next() ;
      if (harg != parg) {
	parg->syncCache() ;
	harg->copyCache(parg,kTRUE) ;
	if (!harg->inRange(0)) {
	  return 0 ;
	}
      }
    }
  }

  //Double_t ret =  _dataHist->weight(_histObsList,_intOrder,kFALSE,_cdfBoundaries) ;  
  Double_t ret =  _dataHist->weight(_histObsList,_intOrder,kTRUE,_cdfBoundaries) ;  
  return ret ;
}

////////////////////////////////////////////////////////////////////////////////
/// Only handle case of maximum in all variables

Int_t RooHistFunc::getMaxVal(const RooArgSet& vars) const 
{
  RooAbsCollection* common = _depList.selectCommon(vars) ;
  if (common->getSize()==_depList.getSize()) {
    delete common ;
    return 1;
  }
  delete common ;
  return 0 ;
}

////////////////////////////////////////////////////////////////////////////////

Double_t RooHistFunc::maxVal(Int_t code) const 
{
  R__ASSERT(code==1) ;

  Double_t max(-1) ;
  for (Int_t i=0 ; i<_dataHist->numEntries() ; i++) {
    _dataHist->get(i) ;
    Double_t wgt = _dataHist->weight() ;
    if (wgt>max) max=wgt ;
  }

  return max*1.05 ;
}

////////////////////////////////////////////////////////////////////////////////
/// Return the total volume spanned by the observables of the RooDataHist

Double_t RooHistFunc::totVolume() const
{
  // Return previously calculated value, if any
  if (_totVolume>0) {
    return _totVolume ;
  }
  _totVolume = 1. ;
  TIterator* iter = _depList.createIterator() ;
  RooAbsArg* arg ;
  while((arg=(RooAbsArg*)iter->Next())) {
    RooRealVar* real = dynamic_cast<RooRealVar*>(arg) ;
    if (real) {
      _totVolume *= (real->getMax()-real->getMin()) ;
    } else {
      RooCategory* cat = dynamic_cast<RooCategory*>(arg) ;
      if (cat) {
	_totVolume *= cat->numTypes() ;
      }
    }
  }
  delete iter ;
  return _totVolume ;
}

namespace {
    bool fullRange(const RooAbsArg& x, const RooAbsArg& y ,const char* range)
    {
      const RooAbsRealLValue *_x = dynamic_cast<const RooAbsRealLValue*>(&x);
      const RooAbsRealLValue *_y = dynamic_cast<const RooAbsRealLValue*>(&y);
      if (!_x || !_y) return false;
      if (!range || !strlen(range) || !_x->hasRange(range) ||
	  _x->getBinningPtr(range)->isParameterized()) {
	// parameterized ranges may be full range now, but that might change,
	// so return false
	if (range && strlen(range) && _x->getBinningPtr(range)->isParameterized())
	    return false;
	return (_x->getMin() == _y->getMin() && _x->getMax() == _y->getMax());
      }
      return (_x->getMin(range) == _y->getMin() && _x->getMax(range) == _y->getMax());
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Determine integration scenario. If no interpolation is used,
/// RooHistFunc can perform all integrals over its dependents
/// analytically via partial or complete summation of the input
/// histogram. If interpolation is used, only the integral
/// over all RooHistPdf observables is implemented.

Int_t RooHistFunc::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{

  // First make list of pdf observables to histogram observables
  // and select only those for which the integral is over the full range

  RooFIter it = _depList.fwdIterator();
  RooFIter jt = _histObsList.fwdIterator();
  Int_t code = 0, frcode = 0, n = 0;
  for (RooAbsArg *pa = 0, *ha = 0; (pa = it.next()) && (ha = jt.next()); ++n) {
    if (allVars.find(*pa)) {
      code |= 2 << n;
      analVars.add(*pa);
      if (fullRange(*pa, *ha, rangeName)) {
	frcode |= 2 << n;
      }
    }
  }

  if (code == frcode) {
    // integrate over full range of all observables - use bit 0 to indicate
    // full range integration over all observables
    code |= 1;
  }
  // Disable partial analytical integrals if interpolation is used, and we
  // integrate over sub-ranges, but leave them enabled when we integrate over
  // the full range of one or several variables
  if (_intOrder > 1 && !(code & 1)) {
    analVars.removeAll();
    return 0;
  }
  return (code >= 2) ? code : 0;
}



////////////////////////////////////////////////////////////////////////////////
/// Return integral identified by 'code'. The actual integration
/// is deferred to RooDataHist::sum() which implements partial
/// or complete summation over the histograms contents

Double_t RooHistFunc::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  // Simplest scenario, full-range integration over all dependents
  if (((2 << _histObsList.getSize()) - 1) == code) {
    return _dataHist->sum(kFALSE);
  }

  // Partial integration scenario, retrieve set of variables, calculate partial
  // sum, figure out integration ranges (if needed)
  RooArgSet intSet;
  std::map<const RooAbsArg*, std::pair<Double_t, Double_t> > ranges;
  RooFIter it = _depList.fwdIterator();
  RooFIter jt = _histObsList.fwdIterator();
  Int_t n(0);
  for (RooAbsArg *pa = 0, *ha = 0; (pa = it.next()) && (ha = jt.next()); ++n) {
    if (code & (2 << n)) {
      intSet.add(*ha);
    }
    if (!(code & 1)) {
      RooAbsRealLValue* rlv = dynamic_cast<RooAbsRealLValue*>(pa);
      if (rlv) {
	const RooAbsBinning* binning = rlv->getBinningPtr(rangeName);
	if (rangeName && rlv->hasRange(rangeName)) {
	  ranges[ha] = std::make_pair(
	      rlv->getMin(rangeName), rlv->getMax(rangeName));
	} else if (binning) {
	  if (!binning->isParameterized()) {
	    ranges[ha] = std::make_pair(
		binning->lowBound(), binning->highBound());
	  } else {
	    ranges[ha] = std::make_pair(
		binning->lowBoundFunc()->getVal(), binning->highBoundFunc()->getVal());
	  }
	}
      }
    }
    // WVE must sync hist slice list values to pdf slice list
    // Transfer values from
    if (ha != pa) {
      pa->syncCache();
      ha->copyCache(pa,kTRUE);
    }
  }

  Double_t ret = (code & 1) ?
    _dataHist->sum(intSet,_histObsList,kFALSE,kFALSE) :
    _dataHist->sum(intSet,_histObsList,kFALSE,kFALSE, ranges);

    //_dataHist->sum(intSet,_histObsList,kTRUE,kTRUE) :
    //_dataHist->sum(intSet,_histObsList,kFALSE,kTRUE, ranges);

  
  //    cout << "intSet = " << intSet << endl ;
  //    cout << "slice position = " << endl ;
  //    _histObsList.Print("v") ;
  //    cout << "RooHistPdf::ai(" << GetName() << ") code = " << code << " ret = " << ret << endl ;
  
  return ret ;
}



////////////////////////////////////////////////////////////////////////////////
/// Return sampling hint for making curves of (projections) of this function
/// as the recursive division strategy of RooCurve cannot deal efficiently
/// with the vertical lines that occur in a non-interpolated histogram

list<Double_t>* RooHistFunc::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const
{
  // No hints are required when interpolation is used
  if (_intOrder>1) {
    return 0 ;
  }


  // Find histogram observable corresponding to pdf observable
  RooAbsArg* hobs(0) ;
  _histObsIter->Reset() ;
  _pdfObsIter->Reset() ;
  RooAbsArg* harg, *parg ;
  while((harg=(RooAbsArg*)_histObsIter->Next())) {
    parg = (RooAbsArg*)_pdfObsIter->Next() ;
    if (string(parg->GetName())==obs.GetName()) {
      hobs=harg ; 
    }
  }
  if (!hobs) {
    return 0 ;
  }

  // Check that observable is in dataset, if not no hint is generated
  RooAbsLValue* lvarg = dynamic_cast<RooAbsLValue*>(_dataHist->get()->find(hobs->GetName())) ;
  if (!lvarg) {
    return 0 ;
  }

  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lvarg->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;

  list<Double_t>* hint = new list<Double_t> ;

  // Widen range slighty
  xlo = xlo - 0.01*(xhi-xlo) ;
  xhi = xhi + 0.01*(xhi-xlo) ;

  Double_t delta = (xhi-xlo)*1e-8 ;
 
  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      hint->push_back(boundaries[i]-delta) ;
      hint->push_back(boundaries[i]+delta) ;
    }
  }

  return hint ;
}


////////////////////////////////////////////////////////////////////////////////
/// Return sampling hint for making curves of (projections) of this function
/// as the recursive division strategy of RooCurve cannot deal efficiently
/// with the vertical lines that occur in a non-interpolated histogram

std::list<Double_t>* RooHistFunc::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const 
{
  // No hints are required when interpolation is used
  if (_intOrder>1) {
    return 0 ;
  }

  // Find histogram observable corresponding to pdf observable
  RooAbsArg* hobs(0) ;
  _histObsIter->Reset() ;
  _pdfObsIter->Reset() ;
  RooAbsArg* harg, *parg ;
  while((harg=(RooAbsArg*)_histObsIter->Next())) {
    parg = (RooAbsArg*)_pdfObsIter->Next() ;
    if (string(parg->GetName())==obs.GetName()) {
      hobs=harg ; 
    }
  }

  // cout << "RooHistFunc::bb(" << GetName() << ") histObs = " << _histObsList << endl ;
  // cout << "RooHistFunc::bb(" << GetName() << ") pdfObs = " << _depList << endl ;

  RooAbsRealLValue* transform(0) ;
  if (!hobs) {

    // Considering alternate: input observable is histogram observable and pdf observable is transformation in terms of it
    RooAbsArg* pobs(0) ;
    _histObsIter->Reset() ;
    _pdfObsIter->Reset() ;
    while((harg=(RooAbsArg*)_histObsIter->Next())) {
      parg = (RooAbsArg*)_pdfObsIter->Next() ;
      if (string(harg->GetName())==obs.GetName()) {
	pobs=parg ; 
	hobs=harg ;
      }
    }

    // Not found, or check that matching pdf observable is an l-value dependent on histogram observable fails
    if (!hobs || !(pobs->dependsOn(obs) && dynamic_cast<RooAbsRealLValue*>(pobs))) {
      cout << "RooHistFunc::binBoundaries(" << GetName() << ") obs = " << obs.GetName() << " hobs is not found, returning null" << endl ;
      return 0 ;
    }

    // Now we are in business - we are in a situation where the pdf observable LV(x), mapping to a histogram observable x
    // We can return bin boundaries by mapping the histogram boundaties through the inverse of the LV(x) transformation
    transform = dynamic_cast<RooAbsRealLValue*>(pobs) ;
  }


  // cout << "hobs = " << hobs->GetName() << endl ;
  // cout << "transform = " << (transform?transform->GetName():"<none>") << endl ;

  // Check that observable is in dataset, if not no hint is generated
  RooAbsArg* xtmp = _dataHist->get()->find(hobs->GetName()) ;
  if (!xtmp) {
    cout << "RooHistFunc::binBoundaries(" << GetName() << ") hobs = " << hobs->GetName() << " is not found in dataset?" << endl ;
    _dataHist->get()->Print("v") ;
    return 0 ;
  }
  RooAbsLValue* lvarg = dynamic_cast<RooAbsLValue*>(_dataHist->get()->find(hobs->GetName())) ;
  if (!lvarg) {
    cout << "RooHistFunc::binBoundaries(" << GetName() << ") hobs = " << hobs->GetName() << " but is not an LV, returning null" << endl ;
    return 0 ;
  }

  // Retrieve position of all bin boundaries
  const RooAbsBinning* binning = lvarg->getBinningPtr(0) ;
  Double_t* boundaries = binning->array() ;

  list<Double_t>* hint = new list<Double_t> ;

  // Construct array with pairs of points positioned epsilon to the left and
  // right of the bin boundaries
  for (Int_t i=0 ; i<binning->numBoundaries() ; i++) {
    if (boundaries[i]>=xlo && boundaries[i]<=xhi) {
      
      Double_t boundary = boundaries[i] ;
      if (transform) {
	transform->setVal(boundary) ;
	//cout << "transform bound " << boundary << " using " << transform->GetName() << " result " << obs.getVal() << endl ;
	hint->push_back(obs.getVal()) ;
      } else {	
	hint->push_back(boundary) ;
      }
    }
  }

  return hint ;
}



////////////////////////////////////////////////////////////////////////////////
/// Check if our datahist is already in the workspace

Bool_t RooHistFunc::importWorkspaceHook(RooWorkspace& ws) 
{  
  std::list<RooAbsData*> allData = ws.allEmbeddedData() ;
  std::list<RooAbsData*>::const_iterator iter ;
  for (iter = allData.begin() ; iter != allData.end() ; ++iter) {
    // If your dataset is already in this workspace nothing needs to be done
    if (*iter == _dataHist) {
      return kFALSE ;
    }
  }

  // Check if dataset with given name already exists
  RooAbsData* wsdata = ws.embeddedData(_dataHist->GetName()) ;

  if (wsdata) {

    // Yes it exists - now check if it is identical to our internal histogram 
    if (wsdata->InheritsFrom(RooDataHist::Class())) {

      // Check if histograms are identical
      if (areIdentical((RooDataHist&)*wsdata,*_dataHist)) {

	// Exists and is of correct type, and identical -- adjust internal pointer to WS copy
	_dataHist = (RooDataHist*) wsdata ;
      } else {

	// not identical, clone rename and import
	TString uniqueName = Form("%s_%s",_dataHist->GetName(),GetName()) ;
	Bool_t flag = ws.import(*_dataHist,RooFit::Rename(uniqueName.Data()),RooFit::Embedded()) ;
	if (flag) {
	  coutE(ObjectHandling) << " RooHistPdf::importWorkspaceHook(" << GetName() << ") unable to import clone of underlying RooDataHist with unique name " << uniqueName << ", abort" << endl ;
	  return kTRUE ;
	}
	_dataHist = (RooDataHist*) ws.embeddedData(uniqueName.Data()) ;
      }

    } else {

      // Exists and is NOT of correct type: clone rename and import
      TString uniqueName = Form("%s_%s",_dataHist->GetName(),GetName()) ;
      Bool_t flag = ws.import(*_dataHist,RooFit::Rename(uniqueName.Data()),RooFit::Embedded()) ;
      if (flag) {
	coutE(ObjectHandling) << " RooHistPdf::importWorkspaceHook(" << GetName() << ") unable to import clone of underlying RooDataHist with unique name " << uniqueName << ", abort" << endl ;
	return kTRUE ;
      }
      _dataHist = (RooDataHist*) ws.embeddedData(uniqueName.Data()) ;
      
    }
    return kFALSE ;
  }
  
  // We need to import our datahist into the workspace
  ws.import(*_dataHist,RooFit::Embedded()) ;
  
  // Redirect our internal pointer to the copy in the workspace
  _dataHist = (RooDataHist*) ws.embeddedData(_dataHist->GetName()) ;
  return kFALSE ;
}


////////////////////////////////////////////////////////////////////////////////

Bool_t RooHistFunc::areIdentical(const RooDataHist& dh1, const RooDataHist& dh2) 
{
  if (fabs(dh1.sumEntries()-dh2.sumEntries())>1e-8) return kFALSE ;
  if (dh1.numEntries() != dh2.numEntries()) return kFALSE ;
  for (int i=0 ; i < dh1.numEntries() ; i++) {
    dh1.get(i) ;
    dh2.get(i) ;
    if (fabs(dh1.weight()-dh2.weight())>1e-8) return kFALSE ;
  }
  if (!(RooNameSet(*dh1.get())==RooNameSet(*dh2.get()))) return kFALSE ;
  return kTRUE ;
}



////////////////////////////////////////////////////////////////////////////////
/// Stream an object of class RooHistFunc.

void RooHistFunc::Streamer(TBuffer &R__b)
{
   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RooHistFunc::Class(),this);
      // WVE - interim solution - fix proxies here
      //_proxyList.Clear() ; // DANILO
      //registerProxy(_depList) ; // DANILO
   } else {
      R__b.WriteClassBuffer(RooHistFunc::Class(),this);
   }
}


////////////////////////////////////////////////////////////////////////////////
/// Schema evolution: if histObsList wasn't filled from persistence (v1)
/// then fill it here. Can't be done in regular schema evolution in LinkDef
/// as _depList content is not guaranteed to be initialized there

void RooHistFunc::ioStreamerPass2() 
{
  if (_histObsList.getSize()==0) {
    _histObsList.addClone(_depList) ;
  }
}

