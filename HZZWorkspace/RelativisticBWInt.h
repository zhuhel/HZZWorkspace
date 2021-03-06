/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef RELATIVISTIC_BWInt_H
#define RELATIVISTIC_BWInt_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RelativisticBWInt : public RooAbsPdf {
public:
  RelativisticBWInt() {} ;
  RelativisticBWInt(const char *name, const char *title,
	      RooAbsReal& _x,
                 RooAbsReal& _mH,
                 RooAbsReal& _w_fr,
                 RooAbsReal& _a0acc,
                 RooAbsReal& _a1acc,
                 RooAbsReal& _a2acc,
                 RooAbsReal& _a3acc,
                 RooAbsReal& _a0int,
                 RooAbsReal& _a1int,
                 RooAbsReal& _a2int,
                 RooAbsReal& _a3int,
                 RooAbsReal& _a4int,
                 RooAbsReal& _b0int,
                 RooAbsReal& _b1int,
                 RooAbsReal& _b2int,
                 RooAbsReal& _b3int,
                 RooAbsReal& _b4int,
                 RooAbsReal& _kappa);
  RelativisticBWInt(const RelativisticBWInt& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RelativisticBWInt(*this,newname); }
  inline virtual ~RelativisticBWInt() { }

protected:

  RooRealProxy x ;
  RooRealProxy mH ;
  RooRealProxy w_fr ;
  RooRealProxy a0acc ;
  RooRealProxy a1acc ;
  RooRealProxy a2acc ;
  RooRealProxy a3acc ;
  RooRealProxy a0int ;
  RooRealProxy a1int ;
  RooRealProxy a2int ;
  RooRealProxy a3int ;
  RooRealProxy a4int ;
  RooRealProxy b0int ;
  RooRealProxy b1int ;
  RooRealProxy b2int ;
  RooRealProxy b3int ;
  RooRealProxy b4int ;
  RooRealProxy kappa ;

  Double_t evaluate() const ;
  Double_t L_gg(Double_t x) const ;

private:

  ClassDef(RelativisticBWInt,1) // Your description goes here...
};

#endif
