/* -*- mode: c++ -*- *********************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * authors:                                                                  *
 *  Lydia Brenner (lbrenner@cern.ch), Carsten Burgard (cburgard@cern.ch)     *
 *  Katharina Ecker (kecker@cern.ch), Adam Kaluza      (akaluza@cern.ch)     *
 *****************************************************************************/
/*****************************************************************************
 *  Katharina Ecker (kecker@cern.ch): Copied (3rd May 2017 ) SVN VERSION  645*
 *****************************************************************************/

#ifndef ROO_LAGRANGIAN_MORPH_FUNC
#define ROO_LAGRANGIAN_MORPH_FUNC

#include "RooAbsPdf.h"
#include "RooRealSumPdf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooProduct.h"
#include "TMatrixD.h"

class RooParamHistFunc;
class TPair;

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

class RooLagrangianMorphFunc : public RooAbsPdf {
 public:
  typedef std::map<const std::string,double> ParamSet;
  static bool gAllowExceptions;
  
  RooLagrangianMorphFunc();
  template<class T> RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const std::vector<T>& vertices, const char* basefolder, const RooArgList& folders);
  template<class T> RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const std::vector<T>& vertices, const RooArgList& folders);
  template<class T> RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const std::vector<T>& vertices);
  RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings, const char* basefolder, const RooArgList& folders);
  RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings, const RooArgList& folders);
  RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings);
  RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const char* basefolder, const RooArgList& folders);
  RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooArgList& folders);
  RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName);
  RooLagrangianMorphFunc(const RooLagrangianMorphFunc& other, const char *newName);
 
  virtual ~RooLagrangianMorphFunc();

  virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const override;
  virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& /*obs*/, Double_t /*xlo*/, Double_t /*xhi*/) const override;
  virtual Bool_t isBinnedDistribution(const RooArgSet& obs) const override;
  virtual Double_t evaluate() const override;
  virtual RooAbsPdf::ExtendMode extendMode() const override;
  virtual Double_t expectedEvents(const RooArgSet* nset) const override;
  virtual Double_t expectedEvents(const RooArgSet& nset) const override;
  virtual TObject* clone(const char* newname) const override;
  virtual Double_t getValV(const RooArgSet* set=0) const override;
  virtual Bool_t selfNormalized() const override;


  void setParameters(const char* foldername);
  void setParameters(TH1* paramhist);
  void setParameter(const char* name, double value);
  void setParameters(const ParamSet& params);
  void setParameters(const RooArgList* list);
  double getParameterValue(const char* name) const;
  RooRealVar* getParameter(const char* name) const;
  bool hasParameter(const char* paramname) const;
  bool isParameterUsed(const char* paramname) const;
  bool isParameterConstant(const char* paramname) const;
  void setParameterConstant(const char* paramname, bool constant) const;
  void setParameter(const char* name, double value, double min, double max);
  void setParameter(const char* name, double value, double min, double max, double error);
  void randomizeParameters(double z);
  RooArgList* getParameterSet() const;
  using RooAbsArg::getParameters;
  ParamSet getParameters(const char* foldername) const;
  ParamSet getParameters() const;

  int nParameters() const;
  int nSamples() const;
  int nPolynomials() const;

  bool isCouplingUsed(const char* couplname) const;
  RooArgList* getCouplingSet() const;
  ParamSet getCouplings() const;

  TMatrixD getMatrix() const;
  TMatrixD getInvertedMatrix() const;
  double getCondition() const;
     
  RooRealVar* getObservable() const;
  RooRealSumPdf* getPdf() const;
  RooRealSumPdf* clonePdf() const;
 
  void printEvaluation() const;
  void printCouplings() const;
  void printParameters() const;
  void printParameters(const char* samplename) const;
  void printSamples() const;
  void printFormulas() const;
  void printPhysics() const;

  RooProduct* getSumElement(const char* name) const;
  
  const std::vector<std::string>& getSamples() const;
  
  Double_t expectedEvents() const;
  double expectedUncertainty() const;
  TH1* createTH1(const std::string& name, RooFitResult* r = NULL);
  TH1* createTH1(const std::string& name, bool correlateErrors, RooFitResult* r = NULL);
  
 protected:

  void printAuthors() const;

  void setup(const RooArgSet& operators, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings, bool ownParams = true);
  template<class T> void setup(const RooArgSet& operators, const std::vector<T>& vertices, bool ownParams = true);

  bool _ownParameters = false; 
  
  mutable RooObjCacheManager _cacheMgr; //! The cache manager
  class CacheElem;
 
  void addFolders(const RooArgList& folders);
  
  bool hasCache() const;
  RooLagrangianMorphFunc::CacheElem* getCache(const RooArgSet* nset) const;
 public:
  
  bool updateCoefficients();
  bool useCoefficients(const TMatrixD& inverse);
  bool useCoefficients(const char* filename);
  bool writeCoefficients(const char* filename);
  bool writeFormulas(const char* filename);
  bool writePhysics(const char* filename);
  
  int countContributingFormulas() const;
  RooParamHistFunc* getBaseTemplate();
  RooAbsReal* getSampleWeight(const char* name);
  void printSampleWeights();
  
 protected:
  std::string _fileName;
  std::string _obsName;
  std::string _baseFolder;
  std::vector<std::string>  _folders;
  RooListProxy _operators;
  RooListProxy _observables ; //!
  std::vector<RooListProxy*> _vertices;

  mutable const RooArgSet* _curNormSet ; //! 

 public:
  // these are a coupleof helper functions for use with the Higgs Characterization (HC) Model
  // arXiv: 1306.6464
  static RooArgSet makeHCggFCouplings(RooAbsCollection& kappas);
  static RooArgSet makeHCVBFCouplings(RooAbsCollection& kappas);
  static RooArgSet makeHCHWWCouplings(RooAbsCollection& kappas);
  static RooArgSet makeHCHZZCouplings(RooAbsCollection& kappas);
  static RooArgSet makeHCHllCouplings(RooAbsCollection& kappas);

  static void writeMatrixToFile(const TMatrixD& matrix, const char* fname);
  static void writeMatrixToStream(const TMatrixD& matrix, std::ostream& stream);
  static TMatrixD readMatrixFromFile(const char* fname);
  static TMatrixD readMatrixFromStream(std::istream& stream);

  static RooDataHist* makeDataHistogram(TH1* hist, RooRealVar* observable, const char* histname = NULL);
  static void setDataHistogram(TH1* hist, RooRealVar* observable, RooDataHist* dh);
  static void printDataHistogram(RooDataHist* hist, RooRealVar* obs);

  static int countSamples(std::vector<RooArgList*>& vertices);
  static int countSamples(int nprod, int ndec, int nboth);

  static TPair* makeCrosssectionContainer(double xs, double unc);

  ClassDefOverride(RooLagrangianMorphFunc,1) 

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // RooLagrangianMorphFunc
    //
    // The RooLagrangianMorphFunc is a type of RooAbsPdf that allows to morph
    // different input EFT samples to some arbitrary output EFT
    // sample, as long as the desired set of output parameters lie
    // within the realm spanned by the input samples. More
    // specifically, it expects as an input a TFile (or TDirectory)
    // with the following layout:
    //
    // TDirectory 
    //  |-sample1
    //  | |-param_card     // a TH1 which encodes the EFT parameter values used for this sample
    //  | | histogram      // a TH1 with a distribution of some physical observable 
    //  | |-subfolder1     // a subfolder (optional)
    //  | | |-histogram1   // another TH1 with a distribution of some physical observable 
    //  | | |-histogram2   // another TH1 with a distribution of some physical observalbe 
    //  | | ...            // more of these
    //  |-sample2
    //  | |-param_card     // a TH1 which encodes the EFT parameter values used for this sample
    //  | | histogram      // a TH1 with a distribution of the same physical observable as above
    //  | | ...
    //  | ...
    //
    // The RooLagrangianMorphFunc operates on this structure, extracts data
    // and meta-data and produces a morphing result as a RooRealSumPdf
    // consisting of the input histograms with appropriate prefactors.
    //
    // The histograms to be morphed can be accessed via their paths in
    // the respective sample, e.g. using
    //    "histogram"
    // or "subfolder1/histogram1"
    // or "some/deep/path/to/some/subfolder/histname"
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////

};

#define ROOLAGRANGIANMORPHFUNC(CLASSNAME) public:                              \
  CLASSNAME(const char *name, const char *title, const char* fileName, const char* obsName, const char* basefolder, const RooArgList& folders): RooLagrangianMorphFunc(name,title,fileName,obsName,basefolder,folders) {this->makeCouplings();} \
CLASSNAME(const char *name, const char *title, const char* fileName, const char* obsName, const RooArgList& folders) : RooLagrangianMorphFunc(name,title,fileName,obsName,folders){this->makeCouplings();} \
CLASSNAME(const char *name, const char *title, const char* fileName, const char* obsName) : RooLagrangianMorphFunc(name,title,fileName,obsName){this->makeCouplings();} \
CLASSNAME(const CLASSNAME& other, const char* newname) : RooLagrangianMorphFunc(other,newname){ }; \
  CLASSNAME(){ };                                                       \
  virtual ~CLASSNAME(){};                                                \
  virtual TObject* clone(const char* newname) const override { return new CLASSNAME(*this,newname); }

////////////////////////////////////////////////////////////////////////////////////////////////
// DERIVED CLASSES to implement specific PHYSICS ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

class RooHCggfWWMorphFunc : public RooLagrangianMorphFunc {
  ROOLAGRANGIANMORPHFUNC(RooHCggfWWMorphFunc)
  ClassDefOverride(RooHCggfWWMorphFunc,1)
 protected:
  void makeCouplings(){
    RooArgSet kappas("ggfWW");
    this->setup(kappas,this->makeHCggFCouplings(kappas),this->makeHCHWWCouplings(kappas),true);
  }
};

class RooHCvbfWWMorphFunc : public RooLagrangianMorphFunc {
  ROOLAGRANGIANMORPHFUNC(RooHCvbfWWMorphFunc)
  ClassDefOverride(RooHCvbfWWMorphFunc,1)
 protected:
  void makeCouplings(){
    RooArgSet kappas("vbfWW");
    this->setup(kappas,this->makeHCVBFCouplings(kappas),this->makeHCHWWCouplings(kappas),true);
  }
};

class RooHCggfZZMorphFunc : public RooLagrangianMorphFunc {
  ROOLAGRANGIANMORPHFUNC(RooHCggfZZMorphFunc)
  ClassDefOverride(RooHCggfZZMorphFunc,1) 
 protected:
  void makeCouplings(){
    RooArgSet kappas("ggfZZ");
    this->setup(kappas,this->makeHCggFCouplings(kappas),this->makeHCHZZCouplings(kappas),true);
  }
};

class RooHCvbfZZMorphFunc : public RooLagrangianMorphFunc {
  ROOLAGRANGIANMORPHFUNC(RooHCvbfZZMorphFunc)
  ClassDefOverride(RooHCvbfZZMorphFunc,1)
 protected:
  void makeCouplings(){
    RooArgSet kappas("vbfZZ");
    this->setup(kappas,this->makeHCVBFCouplings(kappas),this->makeHCHZZCouplings(kappas),true);
  }  
};

class RooHCvbfMuMuMorphFunc : public RooLagrangianMorphFunc {
  ROOLAGRANGIANMORPHFUNC(RooHCvbfMuMuMorphFunc)
  ClassDefOverride(RooHCvbfMuMuMorphFunc,1)
 protected:
  void makeCouplings(){
    RooArgSet kappas("vbfMuMu");
    this->setup(kappas,this->makeHCVBFCouplings(kappas),this->makeHCHllCouplings(kappas),true);
  }    
};

ClassImp(RooLagrangianMorphFunc)
ClassImp(RooHCggfWWMorphFunc)
ClassImp(RooHCvbfWWMorphFunc) 
ClassImp(RooHCggfZZMorphFunc) 
ClassImp(RooHCvbfZZMorphFunc) 
ClassImp(RooHCvbfMuMuMorphFunc) 

#endif
