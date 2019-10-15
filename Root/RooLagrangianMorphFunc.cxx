// #define USE_UBLAS 1
#undef USE_UBLAS
/*****************************************************************************
 * Project: RooLagrangianMorphFunc                                                  *
 *                                                                           *
 * authors:                                                                  *
 *  Lydia Brenner (lbrenner@cern.ch), Carsten Burgard (cburgard@cern.ch)     *
 *  Katharina Ecker (kecker@cern.ch), Adam Kaluza      (akaluza@cern.ch)     *
 *****************************************************************************/
/*****************************************************************************
 *  Katharina Ecker (kecker@cern.ch): Copied (3rd May 2017 ) SVN VERSION  645*
 *****************************************************************************/
#include "HZZWorkspace/RooLagrangianMorphFunc.h"

#include "Riostream.h"

// RooFit includes
#include "RooDataHist.h"
#include "RooHistFunc.h"
#include "RooParamHistFunc.h"
#include "RooRealSumPdf.h"
#include "RooAbsArg.h"
#include "RooStringVar.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooHistFunc.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooFitResult.h"
#include "RooHistConstraint.h"
#include "RooUniformBinning.h"
#include "RooBinning.h"

// plain ROOT includes
#include "TH1.h"
#include "TParameter.h"
#include "TFile.h"
#include "TKey.h"
#include "TFolder.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TMatrixD.h"

// stl includes
#include <map>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstddef>
#include <cmath>
#include <iostream>
#include <limits>

#include <typeinfo>

// #define _DEBUG_

// various preprocessor helpers
#define NaN std::numeric_limits<double>::quiet_NaN()
#define NODEBUG(arg) std::cout << arg << std::endl;
#ifdef _DEBUG_
#define DEBUG(arg) std::cout << arg << std::endl;
//cxcoutD(Eval) << arg << std::endl
#else
#define DEBUG(arg)
#endif
bool RooLagrangianMorphFunc::gAllowExceptions = true;
#define ERROR(arg){                                                     \
  if(RooLagrangianMorphFunc::gAllowExceptions){                                \
    std::stringstream err; err << arg << std::endl; throw(std::runtime_error(err.str())); \
  } else {                                                              \
    std::cerr << arg << std::endl;                                      \
  }}
#define INFO(arg) std::cout << arg << std::endl;

typedef std::vector<std::vector<std::vector<std::vector<bool>>>> opsarray;

template<class MatrixT>
inline size_t size(const MatrixT& matrix);
template <> inline size_t size<TMatrixD> (const TMatrixD& mat){
  // retrieve the size of a square matrix
  return mat.GetNrows();
}
using namespace std;

#ifdef USE_UBLAS
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <boost/numeric/ublas/symmetric.hpp> //inc diag
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/numeric/ublas/lu.hpp>
#pragma GCC diagnostic pop
typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100> > SuperFloat;
typedef std::numeric_limits<boost::multiprecision::cpp_dec_float<100> > SuperFloatPrecision;
typedef boost::numeric::ublas::matrix<SuperFloat> Matrix;
#else
typedef std::numeric_limits<double> SuperFloatPrecision;
#endif

template<class MatrixT>
inline void writeMatrixToStreamT(const MatrixT& matrix, std::ostream& stream){
  // write a matrix to a stream
  for(size_t i=0; i<size(matrix); ++i){
    for(size_t j=0; j<size(matrix); ++j){
#ifdef USE_UBLAS
      stream << std::setprecision(SuperFloatPrecision::max_digits10) << matrix(i,j) << "\t";
#else
      stream << matrix(i,j) << "\t";
#endif
    }
    stream << std::endl;
  }
}
template<class MatrixT>
inline void writeMatrixToFileT(const MatrixT& matrix, const char* fname){
  // write a matrix to a text file
  std::ofstream of(fname);
  if(!of.good()){
    ERROR("unable to read file '"<<fname<<"'!");
  }
  writeMatrixToStreamT(matrix,of);
  of.close();
}

// boost includes
#ifdef USE_UBLAS
inline void printMatrix(const Matrix& mat){
  // write a matrix
  for(size_t i=0; i<mat.size1(); ++i){
    for(size_t j=0; j<mat.size2(); ++j){
      std::cout << std::setprecision(SuperFloatPrecision::max_digits10) << mat(i,j) << " ,\t";
    }
    std::cout << std::endl;
  }
}
template <> inline size_t size<Matrix> (const Matrix& matrix){
  // retrieve the size of a square matrix
  return matrix.size1();
}
inline Matrix diagMatrix(size_t n){
  // create a new diagonal matrix of size n
  return boost::numeric::ublas::identity_matrix<SuperFloat>(n);
}
inline TMatrixD makeRootMatrix(const Matrix& in){
  // convert a matrix into a TMatrixD
  size_t n = size(in);
  TMatrixD mat(n,n);
  for(size_t i=0; i<n; ++i){
    for(size_t j=0; j<n; ++j){
      mat(i,j) = (double)(in(i,j));
    }
  }
  return mat;
}
inline Matrix makeSuperMatrix(const TMatrixD& in){
  // convert a TMatrixD into a Matrix
  size_t n = in.GetNrows();
  Matrix mat(n,n);
  for(size_t i=0; i<n; ++i){
    for(size_t j=0; j<n; ++j){
      mat(i,j) = in(i,j);
    }
  }
  return mat;
}
inline SuperFloat invertMatrix(const Matrix& matrix, Matrix& inverse){
  // calculate the inverse of a matrix, returning the condition
  boost::numeric::ublas::permutation_matrix<size_t> pm(size(matrix));
  SuperFloat mnorm = norm_inf(matrix);
  Matrix lu(matrix);
  try {
    int res = lu_factorize(lu,pm);
    if( res != 0 ){
      printMatrix(matrix);
      ERROR("Error: matrix is not invertable!");
    }
    // backsubstitute to get the inverse
    lu_substitute(lu, pm, inverse);
  } catch (boost::numeric::ublas::internal_logic& error){
    ERROR("Error: matrix is not invertable!");
  }
  SuperFloat inorm = norm_inf(inverse);
  SuperFloat condition = mnorm * inorm;
  return condition;
}
inline Matrix operator* (const Matrix&m, const Matrix& otherM){
  return prod(m,otherM);
}
#else
typedef double SuperFloat;
#include "TDecompLU.h"
typedef TMatrixD Matrix;
inline TMatrixD makeRootMatrix(const Matrix& in){
  // convert a matrix into a TMatrixD
  return TMatrixD(in);
}
inline Matrix makeSuperMatrix(const TMatrixD& in){
  // convert a TMatrixD into a Matrix
  return in;
}
inline Matrix diagMatrix(size_t n){
  // create a new diagonal matrix of size n
  TMatrixD mat(n,n);
  mat.UnitMatrix();
  return mat;
}
inline void printMatrix(const TMatrixD& mat){
  // write a matrix
  writeMatrixToStreamT(mat,std::cout);
}
inline double invertMatrix(const Matrix& matrix, Matrix& inverse){
  // calculate the inverse of a matrix, returning the condition
  TDecompLU lu(matrix);
  bool status = lu.Invert(inverse);
  // check if the matrix is invertible
  if(!status){
    printMatrix(matrix);
    ERROR("Error: matrix is not invertable!");
  }
  double condition = lu.GetCondition();
  const size_t n = size(inverse);
  // sanitize numeric problems
  for(size_t i= 0; i<n; ++i)
    for(size_t j=0; j<n; ++j)
      if(fabs(inverse(i,j)) < 1e-9) inverse(i,j) = 0;
  return condition;
}
#endif

RooDataHist* RooLagrangianMorphFunc::makeDataHistogram(TH1* hist, RooRealVar* observable, const char* histname){
  // convert a TH1 into a RooDataHist
  RooDataHist* dh = new RooDataHist(histname ? histname : hist->GetName(),histname ? histname : hist->GetName(),*observable);
  setDataHistogram(hist,observable,dh);
  return dh;
}

void RooLagrangianMorphFunc::setDataHistogram(TH1* hist, RooRealVar* observable, RooDataHist* dh){
  // set the values of a RooDataHist to those of a TH1
  int nrBins = observable->getBins();
  for (int i=0;i<nrBins;i++) {
    observable->setBin(i);
    dh->set(*observable,hist->GetBinContent(i+1),hist->GetBinError(i+1));
    dh->get(i);
    DEBUG("dh = " << dh->weight() << " +/- " << sqrt(dh->weightSquared()) << ", hist=" <<  hist->GetBinContent(i+1) << " +/- " << hist->GetBinError(i+1));
  }
}


void RooLagrangianMorphFunc::printDataHistogram(RooDataHist* hist, RooRealVar* obs){
  // print the contents of a RooDataHist
  for(Int_t i=0; i<obs->getBins(); ++i){
    hist->get(i);
    obs->setBin(i);
    std::cout << hist->weight() << " +/- " << hist->weightSquared() << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////
// LOCAL FUNCTIONS AND DEFINITIONS ////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// anonymous namespace to prohibit use of these functions outside the class itself
namespace {

  ///////////////////////////////////////////////////////////////////////////////
  // HELPERS ////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  typedef std::map<const std::string,RooLagrangianMorphFunc::ParamSet > ParamMap;
  typedef std::vector<std::vector<bool> > VertexMap;
  typedef std::vector<std::vector<int> > MorphFuncPattern;
  typedef std::map<int,RooFormulaVar*> FormulaList;

  //_____________________________________________________________________________

  inline bool begins_with(const std::string& input, const std::string& match){
    // check if a std::string begins with the given character set
    return input.size() >= match.size()
      && equal(match.begin(), match.end(), input.begin());
  }

  //_____________________________________________________________________________

  template<class A,class B> inline void assignElement(A& a,const B& b){
    // this is a workaround for the missing implicit conversion from SuperFloat<>double
    a = static_cast<A>(b);
  }
  template<class MatrixT>
  inline MatrixT readMatrixFromStreamT(std::istream& stream){
    // read a matrix from a stream
    std::vector<std::vector<SuperFloat> > matrix;
    std::vector<SuperFloat> line;
    while(!stream.eof()){
      if(stream.peek() == '\n'){
        stream.get();
        stream.peek();
        continue;
      }
      SuperFloat val;
      stream>>val;
      line.push_back(val);
      while(stream.peek() == ' ' || stream.peek() == '\t'){
        stream.get();
      }
      if(stream.peek() == '\n'){
        matrix.push_back(line);
        line.clear();
      }
    }
    MatrixT retval(matrix.size(),matrix.size());
    for(size_t i=0; i<matrix.size(); ++i){
      if(matrix[i].size() != matrix.size()){
        ERROR("matrix read from stream doesn't seem to be square!");
      }
      for(size_t j=0; j<matrix[i].size(); ++j){
        assignElement(retval(i,j),matrix[i][j]);
      }
    }
    return retval;
  }
  template<class MatrixT>
  inline MatrixT readMatrixFromFileT(const char* fname){
    // read a matrix from a text file
    std::ifstream in(fname);
    if(!in.good()){
      ERROR("unable to read file '"<<fname<<"'!");
    }
    MatrixT mat = readMatrixFromStreamT<MatrixT>(in);
    in.close();
    return mat;
  }

  //_____________________________________________________________________________

  inline RooLagrangianMorphFunc::ParamSet readParamSet(TH1* h_pc){
    // convert a TH1* param hist into the corresponding ParamSet object
    RooLagrangianMorphFunc::ParamSet point;
    // loop over all bins of the param_card histogram
    for(int ibx = 1; ibx <= h_pc->GetNbinsX(); ++ibx){
      // read the value of one parameter
      const std::string s_coup(h_pc->GetXaxis()->GetBinLabel(ibx));
      double coup_val = h_pc->GetBinContent(ibx);
      // add it to the map
      if(!s_coup.empty()){
        point[s_coup] = coup_val;
      }
    }
    return point;
  }

  //_____________________________________________________________________________

  inline TH1F* getParamHist(TDirectory* file, const std::string& name){
    // retrieve a param_hist from a certain subfolder 'name' of the file
    TFolder* f_tmp = dynamic_cast<TFolder*>(file->Get(name.c_str()));
    if(!f_tmp) ERROR("unable to retrieve folder '"<<name<<"' from file '"<<file->GetName()<<"'!");
    // retrieve the histogram param_card which should live directly in the folder
    TH1F* h_pc = dynamic_cast<TH1F*>(f_tmp->FindObject("param_card"));
    if(!h_pc){
      ERROR("unable to retrieve param_card histo from folder '"<<name<<"'");
      return NULL;
    }
    DEBUG("found param_card for '" << name << "'");
    return h_pc;
  }

  //_____________________________________________________________________________

  inline RooLagrangianMorphFunc::ParamSet readParamSet(TDirectory* file, const std::string& name){
    TH1F* h_pc = getParamHist(file,name);
    return readParamSet(h_pc);
  }

  //_____________________________________________________________________________

  inline ParamMap readParamSets(TDirectory* f, const std::vector<std::string>& names){
    // retrieve the param_hists file and return a map of the parameter values
    // by providing a list of names, only the param_hists of those subfolders are read
    // leaving the list empty is interpreted as meaning 'read everyting'

    ParamMap inputFiles;
    // if the list of names is empty, we assume that this means 'all'
    // loop over all folders in the file
    for(size_t i=0; i<names.size(); i++){
      const std::string name(names[i]);
      // actually read an individual param_hist
      DEBUG("reading param hist '" << name << "'!");
      inputFiles[name] = readParamSet(f,name);
    }

    // return the map of all parameter values found for all samples
    return inputFiles;
  }

  //_____________________________________________________________________________

  inline TDirectory* openFile(const std::string& filename){
    // open the file and return a file pointer
    if(filename.empty()){
      return gDirectory;
    } else {
      DEBUG("opening file '" << filename << "'");
      TFile *file= TFile::Open(filename.c_str(),"READ");
      if (!file|| !file->IsOpen()) {
        if(file) delete file;
        ERROR("could not open file '"<<filename<<"'!");
      }
      return file;
    }
  }

  //_____________________________________________________________________________

  inline void closeFile(TDirectory*& d){
    // open the file and return a file pointer
    TFile* f = dynamic_cast<TFile*>(d);
    if(f){
      f->Close();
      delete f;
      d=NULL;
    }
  }

  //_____________________________________________________________________________

  inline ParamMap readParamSets(const std::string& name, const std::vector<std::string>& names){
    // retrieve the param_hists file and return a map of the parameter values
    // by providing a list of names, only the param_hists of those subfolders are read
    // leaving the list empty is interpreted as meaning 'read everyting'
    TDirectory* f = openFile(name.c_str());
    if(!f) ERROR("unable to open file '"<<name<<"'!");
    ParamMap retval = readParamSets(f,names);
    closeFile(f);
    return retval;
  }


  //_____________________________________________________________________________

  inline RooLagrangianMorphFunc::ParamSet readParamSet(const std::string& fname, const std::string& name){
    // build the ParamSet for a sample in a file
    TDirectory* f = openFile(fname.c_str());
    if(!f) ERROR("unable to open file '"<<fname<<"'!");
    RooLagrangianMorphFunc::ParamSet retval = readParamSet(f,name);
    closeFile(f);
    return retval;
  }

  //_____________________________________________________________________________

  inline void extractOperators(const RooAbsCollection& couplings, RooAbsCollection& operators){
    // extract the operators from a list of couplings
    DEBUG("extracting operators from '" << couplings.GetName() << "'...");
    RooAbsArg* obj;
    RooFIter itr(couplings.fwdIterator());
    while((obj = itr.next())){
      RooFormulaVar* formulacoupling = dynamic_cast<RooFormulaVar*>(obj);
      if(formulacoupling){
        RooAbsCollection* c = formulacoupling->getVariables();
        DEBUG("looking at sub-formula '"<< formulacoupling->GetName() << "' @" << formulacoupling << " with components @" << c);
#ifdef _DEBUG_
        c->Print();
#endif
        extractOperators(*c,operators);
      } else {
        operators.add(*obj);
      }
    }
  }

  //_____________________________________________________________________________

  inline void extractCouplings(const RooAbsCollection& inCouplings, RooAbsCollection& outCouplings){
    // extract the couplings from a given set and copy them to a new one
    RooAbsArg* obj;
    RooFIter itr(inCouplings.fwdIterator());
    while((obj = itr.next())){
      if(!outCouplings.find(obj->GetName())){
        DEBUG("adding parameter " << obj->GetName());
        outCouplings.add(*obj);
      }
    }
  }

  //_____________________________________________________________________________

  inline RooAbsArg& get(RooAbsCollection& operators, const char* name, double defaultval=0){
    // find and, if necessary, create a parameter from a list
    RooAbsArg* kappa = operators.find(name);
    if(kappa) return *kappa;
    RooRealVar* newKappa = new RooRealVar(name,name,defaultval);
    double minVal = 0.9*defaultval;
    double maxVal = 1.1*defaultval;
    newKappa->setRange(std::min(minVal,maxVal),std::max(minVal,maxVal));
    newKappa->setConstant(false);
    operators.add(*newKappa);
    return *newKappa;
  }
  inline RooAbsArg& get(RooAbsCollection& operators, const std::string& name, double defaultval=0){
    // find and, if necessary, create a parameter from a list
    return get(operators,name.c_str(),defaultval);
  }

  //_____________________________________________________________________________

  inline bool setParam(RooRealVar* p, double val, bool force){
    bool ok = true;
    if(val > p->getMax()){
      if(force){
        p->setMax(val);
      } else {
        std::cerr << "ERROR: parameter " << p->GetName() << " out of bounds: " << val << " > " << p->getMax() << std::endl;
        ok=false;
      }
    } else if(val < p->getMin()){
      if(force){
        p->setMin(val);
      } else {
        std::cerr << "ERROR: parameter " << p->GetName() << " out of bounds: " << val << " < " << p->getMin() << std::endl;
        ok=false;
      }
    }
    if(ok) p->setVal(val);
    return ok;
  }

  //_____________________________________________________________________________

  inline bool setParams(const RooLagrangianMorphFunc::ParamSet& point,const RooAbsCollection& args,bool force=false){
    // set all parameters to the values in the param_card histogram
    bool ok = true;
    for(auto paramit=point.begin(); paramit!=point.end(); ++paramit){
      // loop over all the parameters
      const std::string param(paramit->first);
      // retrieve them from the map
      RooRealVar* p = dynamic_cast<RooRealVar*>(args.find(param.c_str()));
      if(!p) continue;
      // set them to their nominal value
      ok = setParam(p,paramit->second,force) && ok;
    }
    return ok;
  }

  //_____________________________________________________________________________

  inline bool setParams(TH1* hist,const RooAbsCollection& args,bool force=false){
    // set all parameters to the values in the param_card histogram
    bool ok = true;
    TAxis* ax = hist->GetXaxis();
    for(int i=1; i<=ax->GetNbins(); ++i){
      // loop over all the parameters
      RooRealVar* p = dynamic_cast<RooRealVar*>(args.find(ax->GetBinLabel(i)));
      if(!p) continue;
      // set them to their nominal value
      ok = setParam(p,hist->GetBinContent(i),force) && ok;
    }
    return ok;
  }

  //_____________________________________________________________________________

  inline RooLagrangianMorphFunc::ParamSet getParams(const RooAbsCollection& parameters){
    // create a set of parameters
    RooFIter itr(parameters.fwdIterator());
    TObject* obj;
    RooLagrangianMorphFunc::ParamSet retval;
    while((obj = itr.next())){
      RooRealVar* param = dynamic_cast<RooRealVar*>(obj);
      if(!param) continue;
      retval[param->GetName()] = param->getVal();
    }
    return retval;
  }

  //_____________________________________________________________________________

  inline void adjustParamRanges(const ParamMap& input, RooArgList& args){
    // build the set of parameters
    DEBUG("adjusting parameter set");
    std::map<std::string,bool> isZero;
    for(Int_t i=0; i<args.getSize(); i++){
      const std::string parname(args.at(i)->GetName());
      isZero[parname] = true;
    }
    for(auto sampleit : input){
      auto point(sampleit.second);
      for(auto it=point.begin(); it!=point.end(); ++it){
        const std::string parname = it->first;
        RooRealVar* param = dynamic_cast<RooRealVar*>(args.find(parname.c_str()));
        if(!param) continue;
        double val(fabs(it->second));
        double max(param->getMax());
        double min(param->getMin());
        if(val != 0){
          isZero[parname] = false;
          if(parname[0] == 'k' || parname[0] == 'g'){
            if( val > 0.5*  max )   param->setMax( 2*val);
            if( val > 0.5*(-min))   param->setMin(-2*val);
            param->setConstant(0);
            param->setError(0.01);
          } else if(begins_with(parname,"cos") || begins_with(parname,"sin")){
            param->setMin( -1);
            param->setMax(  1);
            param->setConstant(0);
            param->setError(0.01);
          } else {
            if( val > 0.9*  max )   param->setMax( 1.1*val);
            if( val < 1.1*  min )   param->setMin( 0.9*val);
            param->setConstant(0);
            param->setError(0.01);
          }
        }
      }
    }
    for(Int_t i=0; i<args.getSize(); i++){
      RooRealVar* param = dynamic_cast<RooRealVar*>(args.at(i));
      if(!param) continue;
      const std::string parname(param->GetName());
      if(isZero[parname]){
        DEBUG("setting parameter to zero: " << param->GetName());
        param->setConstant(1);
      }
    }
  }


  //_____________________________________________________________________________

  TClass* findClass(TFolder* example, const char* name){
    // find the class of an object contained in a TFolder
    if(!example){
      ERROR("unable to find class of object from sample that is NULL");
      return NULL;
    }
    TObject* obj = example->FindObject(name);
    if(!obj){
      ERROR("unable to locate object '"<<name<<"' in folder '" << example->GetName() << "'!");
      return NULL;
    }
    return TClass::GetClass(obj->ClassName());
  }

  //_____________________________________________________________________________

  void collectHistograms(const char* name,TDirectory* file, RooArgList& list_hf, RooRealVar& var, const std::string& varname, const std::string& /*basefolder*/, const ParamMap& inputFiles) {
    // collect the histograms from the input file and convert them to RooStats objects
    DEBUG("building list of histogram functions");
    bool binningOK = false;
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      TFolder* folder = dynamic_cast<TFolder*>(file->Get(sample.c_str()));
      if(!folder){
        ERROR("Error: unable to access data from folder '" << sample << "'!");
        continue;
      }
      TH1* hist = dynamic_cast<TH1*>(folder->FindObject(varname.c_str()));
      if(!hist){
        std::stringstream errstr;
        errstr << "Error: unable to retrieve histogram '" << varname << "' from folder '" << sample << "'. contents are:";
        for(TIter itr = folder->GetListOfFolders()->begin(); itr != folder->GetListOfFolders()->end(); ++itr){
          errstr << " " << itr.Next()->GetName();
        }
        ERROR(errstr.str());
      }
      
      TString histname(sample);
      TString constraintname(sample);
      TString prodname(sample);
      prodname.Append("_");
      prodname.Append(name);
      TString funcname(prodname);
      histname.Prepend("dh_");
      funcname.Prepend("phys_");
      constraintname.Prepend("hc_");
      prodname.Prepend("hc_");

      // TODO: replace the base template with a ParamSetfunc - this seems to be broken currently
      //      if(sample == basefolder){

      //        RooParamHistFunc* hf = new RooParamHistFunc(funcname,funcname,var,*dh);
      //        int nrBins = var.getBins();
      //        var.getBinning().Print();
      //        for (int i=0;i<nrBins;i++) {
      //          var.setBin(i);
      //          std::cout << "phf = " << hf->getVal() << ", dh = " << dh->weight(var) << ", hist=" <<  hist->GetBinContent(i+1) << std::endl;
      //        }
      //        list_hf.add(*hf);
      //      } else

      RooHistFunc* hf = dynamic_cast<RooHistFunc*>(list_hf.find(funcname));
      if(hf){
        hf->setValueDirty(); 
        RooDataHist* dh = &(hf->dataHist());
        RooLagrangianMorphFunc::setDataHistogram(hist,&var,dh);
      } else {
        if(!binningOK){
//          binningOK=true;
//          TVirtualPad* oldPad = gPad;
//          TCanvas* c = new TCanvas();
//          hist->Draw();
//          delete c;
//          gPad=oldPad;
          int n = hist->GetNbinsX();
//           double max = hist->GetXaxis()->GetXmax();
//           double min = hist->GetXaxis()->GetXmin();
//           var.setBinning(RooUniformBinning(min,max,n));
	  std::vector<double> bins;
	  for(int i =1 ; i < n+1 ; ++i){
	    bins.push_back(hist->GetBinLowEdge(i));
	  }
	  bins.push_back(hist->GetBinLowEdge(n)+hist->GetBinWidth(n));
	  var.setBinning(RooBinning(n,&(bins[0])));
          //	  var.getBinning().Print();
        }

        // generate the mean value
        RooDataHist* dh = RooLagrangianMorphFunc::makeDataHistogram(hist,&var,histname.Data());
        hf = new RooHistFunc(funcname,funcname,var,*dh);
        // add it to the list
        list_hf.add(*hf);
      }
      DEBUG("found histogram " << hist->GetName() << " with integral " << hist->Integral());
    }
  }

  template<class T>
  void collectCrosssections(const char* name, TDirectory* file, RooArgList& list_xs, const std::string& varname, const std::string& /*basefolder*/, const ParamMap& inputFiles) {
    // collect the histograms from the input file and convert them to RooStats objects
    DEBUG("building list of histogram functions");
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      TFolder* folder = dynamic_cast<TFolder*>(file->Get(sample.c_str()));
      if(!folder) ERROR("unable to access data from folder '" << sample << "'!");
      TObject* obj = folder->FindObject(varname.c_str());
      TParameter<T>* xsection = NULL;
      TParameter<T>* error = NULL;
      TParameter<T>* p = dynamic_cast<TParameter<T>*>(obj);
      if(p){
        xsection = p;
      }
      TPair* pair = dynamic_cast<TPair*>(obj);
      if(pair){
        xsection = dynamic_cast<TParameter<T>*>(pair->Key());
        error = dynamic_cast<TParameter<T>*>(pair->Value());
      }
      if(!xsection){
        std::stringstream errstr;
        errstr << "Error: unable to retrieve cross section '" << varname << "' from folder '" << sample << "'. contents are:";
        for(TIter itr = folder->GetListOfFolders()->begin(); itr != folder->GetListOfFolders()->end(); ++itr){
          errstr << " " << itr.Next()->GetName();
        }
        ERROR(errstr.str());
      }
      TString objname(sample);
      objname.Prepend("phys_");
      objname.Append("_");
      objname.Append(name);
      RooRealVar* xs = dynamic_cast<RooRealVar*>(list_xs.find(objname));
      if(xs){
        xs->setVal(xsection->GetVal());
      } else {
        xs = new RooRealVar(objname,objname,xsection->GetVal());
        xs->setConstant(true);
        list_xs.add(*xs);
      }
      if(error) xs->setError(error->GetVal());
    }
  }


  ///////////////////////////////////////////////////////////////////////////////
  // formula calculation ////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  void calculateFunctionHelper(const VertexMap& vertexmap, MorphFuncPattern& morphfunc, std::vector<int>& term, int vertexid, bool first){
    // recursive function to determine polynomials
    if(vertexid > 0){
      for(size_t i=0; i<vertexmap[vertexid-1].size(); ++i){
        if(!vertexmap[vertexid-1][i]) continue;
        std::vector<int> newterm(term);
        newterm[i]++;
        if(first){
          calculateFunctionHelper(vertexmap,morphfunc,newterm,vertexid,false);
        } else {
          calculateFunctionHelper(vertexmap,morphfunc,newterm,vertexid-1,true);
        }
      }
    } else {
      bool found = false;
      for(size_t i=0; i<morphfunc.size(); ++i){
        bool thisfound = true;
        for(size_t j=0; j<morphfunc[i].size(); ++j){
          if(morphfunc[i][j] != term[j]){
            thisfound=false;
            break;
          }
        }
        if(thisfound) {
          found = true;
          break;
        }
      }
      if(!found){
        morphfunc.push_back(term);
      }
    }
  }

  MorphFuncPattern calculateFunction(const VertexMap& vertexmap){
    // calculate the morphing function pattern based on a vertex map
    MorphFuncPattern morphfunc; // empty, later result
    int nvtx(vertexmap.size());
    std::vector<int> term(vertexmap[0].size(),0);

    calculateFunctionHelper(vertexmap,morphfunc,term,nvtx,true);
    return morphfunc;
  }

  template<class List>
  inline VertexMap buildVertexMap(const std::vector<List*>& vertices,RooArgList& couplings){
    // build a vertex map based on vertices and couplings appearing
    const int ncouplings = couplings.getSize();
    VertexMap vertexmap;
    for(size_t i=0; i<vertices.size(); ++i){
      const RooAbsCollection* vertex = vertices[i];
      RooFIter citr = couplings.fwdIterator();
      RooAbsReal* coupling;
      std::vector<bool> vertexCouplings(ncouplings,false);
      int idx = -1;
      while((coupling = dynamic_cast<RooAbsReal*>(citr.next()))){
        idx++;
        if(!coupling){
          ERROR("encountered invalid list of couplings in vertex!");
        }
        if(vertex->find(coupling->GetName())){
          vertexCouplings[idx] = true;
        }
      }
      vertexmap.push_back(vertexCouplings);
    }
    return vertexmap;
  }

  template<class MatrixT>
  inline MatrixT buildMatrixT(const ParamMap& inputFiles, const FormulaList& formulas, RooAbsCollection& args){
    // fill the matrix of coefficients
    const size_t dim = inputFiles.size();
    MatrixT matrix(dim,dim);
    int row = 0;
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      // set all vars to value stored in input file
      if(!setParams(sampleit->second,args,true)){
        ERROR("unable to set parameters for sample "<<sample<<"!");
      }
      // loop over all the formulas
      int col = 0;
      for(auto formulait=formulas.begin(); formulait!=formulas.end(); ++formulait){
        RooFormulaVar* formula = formulait->second;
        if(!formula){
          ERROR("Error: invalid formula encountered!");
        }
        matrix(row,col) = formula->getVal();
        col++;
      }
      row++;
    }
    return matrix;
  }

  inline void checkMatrix(const ParamMap& inputFiles, const FormulaList& formulas){
    // check if the matrix is square
    if(inputFiles.size() != formulas.size()){
      std::stringstream ss;
      ss << "ERROR: matrix is not square, consistency check failed: " <<
        inputFiles.size() << " samples, " <<
      formulas.size() << " expressions:" << std::endl;
      ss << "formulas: " << std::endl;
      for(auto formula : formulas){
        ss << formula.second->GetTitle() << std::endl;
      }
      ss << "samples: " << std::endl;
      for(auto sample : inputFiles){
        ss << sample.first << std::endl;
      }
      ERROR(ss.str());
    }
  }

  inline void inverseSanity(const Matrix& matrix, const Matrix& inverse, double& unityDeviation, double& largestWeight){
    DEBUG("multiplying for sanity check");
    Matrix unity(inverse * matrix);
    DEBUG("matrix operations done");

    // check if the entries in the inverted matrix are sensible
    unityDeviation = 0.;
    largestWeight = 0.;
    const size_t dim = size(unity);
    for(size_t i=0; i<dim; ++i){
      for(size_t j=0; j<dim; ++j){
        if(inverse(i,j) > largestWeight){
          largestWeight = (double)inverse(i,j);
        }
        if(fabs(unity(i,j) - (i==j)) > unityDeviation){
          unityDeviation = fabs((double)unity(i,j)) - (i==j);
        }
      }
    }
    DEBUG("found deviation of " << unityDeviation << ", largest weight is " << largestWeight << ".");
  }
 
  template<class List>
  inline void checkNameConflict(const ParamMap& inputFiles, List& args){
    // check for name conflicts between the input samples and an argument set
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      RooAbsArg* arg = args.find(sample.c_str());
      if(arg){
        ERROR("detected name conflict: cannot use sample '" << sample << "' - a parameter with the same name of type '" << arg->ClassName() << "' is present in set '" << args.GetName() << "'!");
      }
    }
  }


  inline FormulaList buildFormulas(const ParamMap& inputFiles, const MorphFuncPattern& morphfunc, RooArgList& couplings, RooArgList& operators){
    // build the formulas corresponding to the given set of input files and the physics process
    // example vbf hww:
    //                        Operators kSM,  kHww, kAww, kHdwR,kHzz, kAzz
    // std::vector<bool> vertexProd  = {true, true, true, true, true, true };
    // std::vector<bool> vertexDecay = {true, true, true, true, false,false};
    // vertexmap.push_back(vertexProd);
    // vertexmap.push_back(vertexDecay);
    const int ncouplings = couplings.getSize();
    std::vector<bool> couplingsZero(ncouplings,true);

    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      if(!setParams(sampleit->second,operators,true)){
        ERROR("unable to set parameters for sample '"<<sample<<"'!");
      }

      RooAbsReal* obj;
      int idx = 0;
      RooFIter itr(couplings.fwdIterator());
      while((obj = dynamic_cast<RooAbsReal*>(itr.next()))){
        if(obj->getVal() != 0){
          if(couplingsZero[idx]){
            DEBUG(obj->GetName() << " is non-zero for sample " << sample << " (idx=" << idx << ")!");
            couplingsZero[idx] = false;
          }
        }
        idx++;
      }
    }

    #ifdef _DEBUG_
    {
      int idx = 0;
      RooAbsReal* obj;
      RooFIter itr(couplings.fwdIterator());
      while((obj = dynamic_cast<RooAbsReal*>(itr.next()))){
        if(couplingsZero[idx]){
          DEBUG(obj->GetName() << " is zero (idx=" << idx << ")");
        } else {
          DEBUG(obj->GetName() << " is non-zero (idx=" << idx << ")");
        }
        idx++;
      }
    }
    #endif

    FormulaList formulas;
    for(size_t i=0; i<morphfunc.size(); ++i){
      std::stringstream ss;
      bool first = true;
      bool isZero = false;
      std::string nullname;
      for(size_t j=0; j<morphfunc[i].size(); ++j){
        const int exponent = morphfunc[i][j];
        if(exponent == 0) continue;
        std::string cname(couplings.at(j)->GetName());
        if(!first) ss << "*";
        if(exponent > 1)
          ss << "TMath::Power(" << cname << "," << exponent << ")";
        else
          ss << cname;
        first = false;
        if(!isZero && couplingsZero[j]){
          isZero = true;
          nullname=cname;
        }
      }
      if(!isZero){
        // build the name
        const TString name = TString::Format("_a%lu",i);
        const std::string formula(ss.str());
        DEBUG("creating formula " << name << ": " << formula);
        // create and add the formula
        formulas[i] = new RooFormulaVar(name.Data(),formula.c_str(),couplings);
      } else {
        DEBUG("killing formula " << ss.str() << " because " << nullname << " is zero");
      }
    }
    return formulas;
  }

  // fancy numerics ahead
#ifdef USE_UBLAS
  //in test code...
  class LinearCombination : public RooAbsReal {
    RooListProxy _actualVars ;
    std::vector<SuperFloat> _coefficients;
    mutable RooArgSet* _nset ;

    mutable SuperFloat _result;
    mutable SuperFloat _tmp;

  public:
    LinearCombination(const char* name) :
      RooAbsReal(name,name),
      _actualVars("actualVars","Variables used by formula expression",this),
      _nset(0)
    {
      // constructor
    }
    LinearCombination(const LinearCombination& other, const char* name) :
      RooAbsReal(other, name),
      _actualVars("actualVars",this,other._actualVars),
      _coefficients(other._coefficients),
      _nset(0)
    {
    }

    ~LinearCombination(){
      // destructor
    }

    virtual TObject* clone(const char* newname) const override {
      return new LinearCombination(*this,newname);
    }

    void add(SuperFloat c,RooAbsReal* t){
      // add a new term
      _actualVars.add(*t);
      _coefficients.push_back(c);
    }

    void setCoefficient(size_t idx,SuperFloat c){
      // set the coefficient with the given index
      this->_coefficients[idx]=c;
    }

    SuperFloat getCoefficient(size_t idx){
      // get the coefficient with the given index
      return this->_coefficients[idx];
    }

    virtual Double_t evaluate() const override {
      // call the evaluation
      this->_result.assign(0.);
      const std::size_t n(this->_actualVars.getSize());
      for(std::size_t i=0;i<n; ++i){
        this->_tmp.assign( (static_cast<const RooAbsReal*>(this->_actualVars.at(i)))->getVal() );
        this->_result += this->_coefficients[i] * this->_tmp;
      }
      return this->_result.convert_to<double>();
    }

    virtual std::list<Double_t>* binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const override{
      // Forward the plot sampling hint from the p.d.f. that defines the observable obs
      RooFIter iter = this->_actualVars.fwdIterator();
      RooAbsReal* func;
      while((func=(RooAbsReal*)iter.next())) {
        list<Double_t>* binb = func->binBoundaries(obs,xlo,xhi);
        if (binb) {
          return binb;
        }
      }
      return 0;
    }
    virtual std::list<Double_t>* plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
      // Forward the plot sampling hint from the p.d.f. that defines the observable obs
      RooFIter iter = this->_actualVars.fwdIterator();
      RooAbsReal* func;
      while((func=(RooAbsReal*)iter.next())){
        list<Double_t>* hint = func->plotSamplingHint(obs,xlo,xhi);
        if (hint) {
          return hint;
        }
      }
      return 0;
    }
  };
#endif

///////////////////////////////////////////////////////////////////////////////

}


RooArgSet RooLagrangianMorphFunc::makeHCggFCouplings(RooAbsCollection& operators) {
  // create the couplings needed for ggF vertices
  DEBUG("creating ggf couplings");
  RooArgSet prodCouplings("ggf");
  DEBUG("adding cosa");
  RooAbsArg& cosa = get(operators,"cosa",1);
  DEBUG("adding Hgg");
  prodCouplings.add(*(new RooFormulaVar("_gHgg" ,"cosa*kHgg",                       RooArgList(cosa,get(operators,"kHgg")))));
  DEBUG("adding Agg");
  prodCouplings.add(*(new RooFormulaVar("_gAgg" ,"sqrt(1-(cosa*cosa))*kAgg",        RooArgList(cosa,get(operators,"kAgg")))));
  return prodCouplings;
}

RooArgSet RooLagrangianMorphFunc::makeHCVBFCouplings(RooAbsCollection& operators) {
  // create the couplings needed for VBF vertices
  RooArgSet prodCouplings("vbf");
  RooAbsArg& cosa = get(operators,"cosa",1);
  RooAbsArg& lambda = get(operators,"Lambda",1000);
  prodCouplings.add(*(new RooFormulaVar("_gSM"  ,"cosa*kSM",                        RooArgList(cosa,get(operators,"kSM")))));
  prodCouplings.add(*(new RooFormulaVar("_gHaa" ,"cosa*kHaa",                       RooArgList(cosa,get(operators,"kHaa")))));
  prodCouplings.add(*(new RooFormulaVar("_gAaa" ,"sqrt(1-(cosa*cosa))*kAaa",        RooArgList(cosa,get(operators,"kAaa")))));
  prodCouplings.add(*(new RooFormulaVar("_gHza" ,"cosa*kHza",                       RooArgList(cosa,get(operators,"kHza")))));
  prodCouplings.add(*(new RooFormulaVar("_gAza" ,"sqrt(1-(cosa*cosa))*kAza",        RooArgList(cosa,get(operators,"kAza")))));
  prodCouplings.add(*(new RooFormulaVar("_gHzz" ,"cosa*kHzz/Lambda",                RooArgList(cosa,get(operators,"kHzz"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gAzz" ,"sqrt(1-(cosa*cosa))*kAzz/Lambda", RooArgList(cosa,get(operators,"kAzz"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gHdz","cosa*kHdz/Lambda",                 RooArgList(cosa,get(operators,"kHdz"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gHww" ,"cosa*kHww/Lambda",                RooArgList(cosa,get(operators,"kHww"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gAww" ,"sqrt(1-(cosa*cosa))*kAww/Lambda", RooArgList(cosa,get(operators,"kAww"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gHdwR","cosa*kHdwR/Lambda",               RooArgList(cosa,get(operators,"kHdwR"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gHdwI","cosa*kHdwI/Lambda",               RooArgList(cosa,get(operators,"kHdwI"),lambda))));
  prodCouplings.add(*(new RooFormulaVar("_gHda","cosa*kHda/Lambda",                 RooArgList(cosa,get(operators,"kHda"),lambda))));
  return prodCouplings;
}

RooArgSet RooLagrangianMorphFunc::makeHCHWWCouplings(RooAbsCollection& operators) {
  // create the couplings needed for HWW vertices
  DEBUG("creating HWW couplings");
  RooArgSet decCouplings("HWW");
  RooAbsArg& cosa = get(operators,"cosa",1);
  RooAbsArg& lambda = get(operators,"Lambda",1000);
  decCouplings.add(*(new RooFormulaVar("_gSM"  ,"cosa*kSM",                        RooArgList(cosa,get(operators,"kSM")))));
  decCouplings.add(*(new RooFormulaVar("_gHww" ,"cosa*kHww/Lambda",                RooArgList(cosa,get(operators,"kHww"),lambda))));
  decCouplings.add(*(new RooFormulaVar("_gAww" ,"sqrt(1-(cosa*cosa))*kAww/Lambda", RooArgList(cosa,get(operators,"kAww"),lambda))));
  decCouplings.add(*(new RooFormulaVar("_gHdwR","cosa*kHdwR/Lambda",               RooArgList(cosa,get(operators,"kHdwR"),lambda))));
  decCouplings.add(*(new RooFormulaVar("_gHdwI","cosa*kHdwI/Lambda",               RooArgList(cosa,get(operators,"kHdwI"),lambda))));
  return decCouplings;
}
RooArgSet RooLagrangianMorphFunc::makeHCHZZCouplings(RooAbsCollection& operators) {
  // create the couplings needed for HZZ vertices
  RooArgSet decCouplings("HZZ");
  RooAbsArg& cosa = get(operators,"cosa",1);
  RooAbsArg& lambda = get(operators,"Lambda",1000);
  decCouplings.add(*(new RooFormulaVar("_gSM"  ,"cosa*kSM",                        RooArgList(cosa,get(operators,"kSM")))));
  decCouplings.add(*(new RooFormulaVar("_gHzz" ,"cosa*kHzz/Lambda",                RooArgList(cosa,get(operators,"kHzz"),lambda))));
  decCouplings.add(*(new RooFormulaVar("_gAzz" ,"sqrt(1-(cosa*cosa))*kAzz/Lambda", RooArgList(cosa,get(operators,"kAzz"),lambda))));
  decCouplings.add(*(new RooFormulaVar("_gHdz","cosa*kHdz/Lambda",                 RooArgList(cosa,get(operators,"kHdz"),lambda))));
  decCouplings.add(*(new RooFormulaVar("_gHaa" ,"cosa*kHaa",                       RooArgList(cosa,get(operators,"kHaa")))));
  decCouplings.add(*(new RooFormulaVar("_gAaa" ,"sqrt(1-(cosa*cosa))*kAaa",        RooArgList(cosa,get(operators,"kAaa")))));
  decCouplings.add(*(new RooFormulaVar("_gHza" ,"cosa*kHza",                       RooArgList(cosa,get(operators,"kHza")))));
  decCouplings.add(*(new RooFormulaVar("_gAza" ,"sqrt(1-(cosa*cosa))*kAza",        RooArgList(cosa,get(operators,"kAza")))));
  decCouplings.add(*(new RooFormulaVar("_gHda","cosa*kHda/Lambda",                 RooArgList(cosa,get(operators,"kHda"),lambda))));
  return decCouplings;
}

RooArgSet RooLagrangianMorphFunc::makeHCHllCouplings(RooAbsCollection& operators) {
  // create the couplings needed for Hll vertices
  RooArgSet decCouplings("Hmumu");
  RooAbsArg& cosa = get(operators,"cosa",1);
  decCouplings.add(*(new RooFormulaVar("_gHll" ,"cosa*kHll",                RooArgList(cosa,get(operators,"kHll")))));
  //  decCouplings.add(*(new RooFormulaVar("_gAll" ,"sqrt(1-(cosa*cosa))*kAll", RooArgList(cosa,get(operators,"kAll")))));
  return decCouplings;
}

///////////////////////////////////////////////////////////////////////////////
// CacheElem magic ////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#define MAXTERMSFORMULA 100

class RooLagrangianMorphFunc::CacheElem : public RooAbsCacheElement {
  public:
    CacheElem(){};
    void operModeHook(RooAbsArg::OperMode) {};
    virtual ~CacheElem();
    virtual RooArgList containedArgs(Action);

    RooRealSumPdf* _sumFunc = 0 ;
    RooArgList* _params = 0 ;
    RooArgList* _couplings = 0 ;
    RooRealVar* _observable = 0 ;
    bool _ownObs = true ;

    FormulaList _formulas;
    RooArgList _phys;
    RooArgList _weights;

    Matrix _matrix;
    Matrix _inverse;
    double _condition;

  public:

#ifndef USE_UBLAS
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
  inline void updateSampleWeights(RooLagrangianMorphFunc* func,const ParamMap& inputFiles){
#ifdef USE_UBLAS
    int sampleidx = 0;
    const size_t n(size(this->_inverse));
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      // build the formula with the correct normalization
      LinearCombination* sampleformula = dynamic_cast<LinearCombination*>(func->getSampleWeight(sample.c_str()));
      if(!sampleformula){
        throw std::runtime_error(TString::Format("unable to access formula for sample '%s'!",sample.c_str()).Data());
      }
      DEBUG("updating formula for sample '" << sample << "'");
      for(size_t formulaidx = 0; formulaidx<n; ++formulaidx){
        const SuperFloat val(this->_inverse(formulaidx,sampleidx));
        DEBUG("   " << formulaidx << ":" << sampleformula->getCoefficient(formulaidx) << " -> " << val);
        sampleformula->setCoefficient(formulaidx,val);
        assert(sampleformula->getCoefficient(formulaidx) == val);
      }
      ++sampleidx;
    }
#else
    throw std::runtime_error("updating sample weights currently not possible without boost!");
#endif
  }
#ifndef USE_UBLAS
#pragma GCC diagnostic pop
#endif

  inline void buildSampleWeights(const char* fname,const ParamMap& inputFiles){
    // build the sampleformulas
    int sampleidx = 0;
    //printMatrix(this->_inverse);
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      const std::string sample(sampleit->first);
      DEBUG("building formula for sample '" << sample << "'");
      TString name_full(sample.c_str());
      name_full.Append("_");
      name_full.Append(fname);
      name_full.Prepend("w_");
      int formulaidx = 0;
      // build the formula with the correct normalization
#ifdef USE_UBLAS
      LinearCombination* sampleformula = new LinearCombination(name_full.Data());
      for(auto formulait=this->_formulas.begin(); formulait!=this->_formulas.end(); ++formulait){
        const SuperFloat val(this->_inverse(formulaidx,sampleidx));
        RooFormulaVar* formula = formulait->second;
        sampleformula->add(val,formula);
        formulaidx++;
      }
      this->_weights.add(*sampleformula);
#else
      std::stringstream ss,ss_full;
      RooArgList sampleformulas;
      bool isNew (true);
      bool isNew_full (true);
      for(auto formulait=this->_formulas.begin(); formulait!=this->_formulas.end(); ++formulait){
        const double val(this->_inverse(formulaidx,sampleidx));
        if(val != 0){
          RooFormulaVar* formula = formulait->second;
          if(!isNew) ss << " + ";
          isNew = false;
          ss << val << "*" << "(" << formula->GetTitle() << ")";
        }
        formulaidx++;
        if((formulaidx % MAXTERMSFORMULA == 0 || formulait == (--this->_formulas.end())) && ss.str().size() != 0){
          TString name(sample.c_str());
          name.Prepend("w_"+std::to_string(formulaidx)+"_");
	  name.Append("_");
          name.Append(fname);
          RooFormulaVar* sampleformula = new RooFormulaVar(name.Data(),ss.str().c_str(),*(this->_couplings));
#ifdef _DEBUG_
          setParams(sampleit->second,*(this->_params),true);
          std::cerr << sample << "_" << formulaidx << ": " << sampleformula->getVal() << " = " << sampleformula->GetTitle() << std::endl;;
#endif
          sampleformulas.add(*sampleformula);
          isNew = true;
          ss.str(std::string());

          if(!isNew_full) ss_full << " + ";
          isNew_full = false;
          ss_full << name;
        }
      }
      if(ss_full.str().size() == 0){
        ERROR("formula received for sample "<<sample<<" is empty!");
      }
      RooFormulaVar* sampleformula_full = new RooFormulaVar(name_full.Data(),ss_full.str().c_str(),sampleformulas);
#ifdef _DEBUG_
      setParams(sampleit->second,*(this->_params),true);
      std::cerr << sample << ": " << sampleformula_full->getVal() << " = " << sampleformula_full->GetTitle() << std::endl;;
#endif
      DEBUG("adding formula '" << sample << "'");
      this->_weights.add(*sampleformula_full);
#endif
      sampleidx++;
    }

    // TODO: for large matrices this takes alot of time, need to implement an option to skip check
    // // check if the results of the sample formulas is sensible
    // for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
    //   const std::string sample(sampleit->first);
    //   RooFIter itr(list_sampleformulas.fwdIterator());
    //   TObject* obj;
    //   setParams(sampleit->second,this->_params);
    //   TString fname(sample.c_str());
    //   fname.Prepend("w_");

    //   std::cout<<"test checking formula: "<<fname<<std::endl;

    //   while((obj = itr.next())){
    //     RooAbsReal* formula = dynamic_cast<RooAbsReal*>(obj);
    //     if(!formula) continue;
    //     std::stringstream ss;
    //     ss<<"checking formula: "<<fname<<" vs. "<<formula->GetName()<<"\t\t\t";
    //     if(fname.CompareTo(formula->GetName()) == 0){
    //       if(fabs(formula->getVal() - 1) > 0.1){
    //         ss << "Warning: weight " << formula->GetName() << "==" << formula->getVal() << " is not close to unity for its own parameter set.";
    //         std::cerr << ss.str() << std::endl;
    //       }
    //     } else {
    //       if(fabs(formula->getVal()) > 0.1){
    //         ss << "Warning: weight " << formula->GetName() << "==" << formula->getVal() << " is not close to zero in the parameter set of " << sample << ".";
    //         std::cerr << ss.str() << std::endl;
    //       }
    //     }
    //   }
    // }
    DEBUG("done building sample weights");
  }

  inline void createBasicObjects(const ParamMap& inputFiles,const char* name, const RooListProxy& operators,const std::vector<RooListProxy*>& vertices,bool adjustParameters){
    // create the basic objects required for the morphing
    TString vname(name);
    vname.ReplaceAll("/","_");
    TString varname(vname);
    varname.Prepend("var_");
    this->_couplings = new RooArgList();
    if (!this->_observable) {
      this->_observable = new RooRealVar(vname,vname,0,1);
    } else {
      DEBUG("createBasicObjects: recycling existing observable object " << this->_observable);
      this->_ownObs = false ;
      if (vname != this->_observable->GetName()) {
        std::cerr << "WARNING: name of existing observable " << this->_observable->GetName() << " does not match expected name " << vname << std::endl ;
      }
    }
    this->_params = new RooArgList(operators);
    if(adjustParameters){
      DEBUG("adjusting parameter ranges");
      adjustParamRanges(inputFiles,*this->_params);
    }
    DEBUG("collecting couplings");
    for(auto vertex : vertices){
      extractCouplings(*vertex,*this->_couplings);
    }
    DEBUG("building vertex map");
    VertexMap vertexmap(buildVertexMap<RooListProxy>(vertices,*this->_couplings));
    DEBUG("calculating pattern for vertexmap of size " << vertexmap.size());
    MorphFuncPattern morphfuncpattern(calculateFunction(vertexmap));
    DEBUG("building formulas");
    this->_formulas = buildFormulas(inputFiles,morphfuncpattern,*this->_couplings,*this->_params);
    DEBUG("checking matrix consistency");
    checkMatrix(inputFiles,this->_formulas);
  }


  //_____________________________________________________________________________

  inline void buildMatrix(const ParamMap& inputFiles){
    // build and invert the morphing matrix
    DEBUG("filling matrix");
    Matrix matrix(buildMatrixT<Matrix>(inputFiles,this->_formulas,*this->_params));
    if(size(matrix) < 1 ){
      ERROR("input matrix is empty, please provide suitable input samples!");
    }
    Matrix inverse(diagMatrix(size(matrix)));
#ifdef _DEBUG_
    printMatrix(matrix);
#endif
    DEBUG("inverting matrix");
    double condition = (double)(invertMatrix(matrix,inverse));
    DEBUG("inverse matrix (condition " << condition << ") is:");
#ifdef _DEBUG_
    printMatrix(inverse);
#endif
    
    double unityDeviation, largestWeight;
    inverseSanity(matrix, inverse, unityDeviation, largestWeight);
    bool weightwarning(largestWeight > 10e7 ? true : false);
    bool unitywarning(unityDeviation > 10e-6 ? true : false);

    // if(unitywarning || weightwarning){
    if(false){
      if(unitywarning){
        std::cerr << "Warning: The matrix inversion seems to be unstable. This can be a result to input samples that are not sufficiently different to provide any morphing power." << std::endl;
      } else if(weightwarning){
        std::cerr << "Warning: Some weights are excessively large. This can be a result to input samples that are not sufficiently different to provide any morphing power." << std::endl;
      }
      std::cerr << "         Please consider the couplings encoded in your samples to cross-check:" << std::endl;
      for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
        const std::string sample(sampleit->first);
        std::cerr << "         " << sample << ": ";
        // set all vars to value stored in input file
        setParams(sampleit->second,*this->_params,true);
        RooFIter itr(this->_couplings->fwdIterator());
        bool first = true;
        RooAbsReal* obj;
        while((obj = dynamic_cast<RooAbsReal*>(itr.next()))){
          if(!first) std::cerr << ", ";
          std::cerr << obj->GetName() << "=" << obj->getVal();
          first = false;
        }
        std::cerr << std::endl;
      }
    }
#ifndef USE_UBLAS
    this->_matrix.ResizeTo(matrix.GetNrows(),matrix.GetNrows());
    this->_inverse.ResizeTo(matrix.GetNrows(),matrix.GetNrows());
#endif
    this->_matrix  = matrix;
    this->_inverse = inverse;
    this->_condition=condition;
  }

  //_____________________________________________________________________________

  inline bool updateMorphingFunction(RooLagrangianMorphFunc* func,TDirectory* file,const ParamMap& inputFiles,
                                     const char* obsName, const std::vector<std::string>& folders,const char* baseFolder){
    // update the morphing function with a new inverse matrix

    // first, receive the new physics objects
    TClass* mode = findClass(dynamic_cast<TFolder*>(file->Get(folders[0].c_str())),obsName);
    if(mode->InheritsFrom(TH1::Class())){
      collectHistograms(func->GetName(), file, this->_phys, *this->_observable, obsName, baseFolder, inputFiles);
    } else if(mode->InheritsFrom(TParameter<double>::Class())){
      collectCrosssections<double>(func->GetName(), file, this->_phys, obsName, baseFolder, inputFiles);
    } else if(mode->InheritsFrom(TParameter<float>::Class())){
      collectCrosssections<float>(func->GetName(), file, this->_phys, obsName, baseFolder, inputFiles);
    } else if(mode->InheritsFrom(TPair::Class())){
      TFolder* folder = (dynamic_cast<TFolder*>(file->Get(folders[0].c_str())));
      TPair* pair = dynamic_cast<TPair*>(folder->FindObject(obsName));
      TParameter<double>* xsec_double = dynamic_cast<TParameter<double>*>(pair->Key());
      if(xsec_double){
        collectCrosssections<double>(func->GetName(), file, this->_phys, obsName, baseFolder, inputFiles);
      } else {
        TParameter<float>* xsec_float = dynamic_cast<TParameter<float>*>(pair->Key());
        if(xsec_float) {
          collectCrosssections<float>(func->GetName(), file, this->_phys, obsName, baseFolder, inputFiles);
        }
        else {
          ERROR("cannot morph objects of class 'TPair' if parameter is not double or float!");
        }
      }
    } else {
      ERROR("cannot morph objects of class '"<<mode->GetName()<<"'!");
      return false;
    }

    // then, update the weights in the morphing function
    updateSampleWeights(func,inputFiles);
    return true;
  }

  //_____________________________________________________________________________

  inline void buildMorphingFunction(const char* name,TDirectory* file,const ParamMap& inputFiles,
                                    const char* obsName, const std::vector<std::string>& folders,const char* baseFolder){
    // build the final morphing function

    // retrieve the physics inputs
    DEBUG("initializing physics inputs");
    TClass* mode = findClass(dynamic_cast<TFolder*>(file->Get(folders[0].c_str())),obsName);
    TString sbw = TString::Format("binWidth_%s",obsName);
    sbw.ReplaceAll("/","_");
    RooRealVar* binWidth = new RooRealVar(sbw.Data(),sbw.Data(),1.);
    if(mode->InheritsFrom(TH1::Class())){
      collectHistograms(name, file, this->_phys, *this->_observable, obsName, baseFolder, inputFiles);
      double bw = this->_observable->numBins()/(this->_observable->getMax() - this->_observable->getMin());
      binWidth->setVal(bw);
    } else if(mode->InheritsFrom(TParameter<double>::Class())){
      collectCrosssections<double>(name, file, this->_phys, obsName, baseFolder, inputFiles);
    } else if(mode->InheritsFrom(TParameter<float>::Class())){
      collectCrosssections<float>(name, file, this->_phys, obsName, baseFolder, inputFiles);
    } else if(mode->InheritsFrom(TPair::Class())){
      TFolder* folder = (dynamic_cast<TFolder*>(file->Get(folders[0].c_str())));
      TPair* pair = dynamic_cast<TPair*>(folder->FindObject(obsName));
      TParameter<double>* xsec_double = dynamic_cast<TParameter<double>*>(pair->Key());
      if(xsec_double){
        collectCrosssections<double>(name, file, this->_phys, obsName, baseFolder, inputFiles);
      } else {
        TParameter<float>* xsec_float = dynamic_cast<TParameter<float>*>(pair->Key());
        if(xsec_float) {
          collectCrosssections<float>(name, file, this->_phys, obsName, baseFolder, inputFiles);
        }
        else {
          ERROR("cannot morph objects of class 'TPair' if parameter is not double or float!");
        }
      }
    } else {
      ERROR("cannot morph objects of class '"<<mode->GetName()<<"'!");
    }
    binWidth->setConstant(true);

    // retrieve the weights
    DEBUG("creating Sample Weights");
    buildSampleWeights(name,inputFiles);

    DEBUG("creating RooProducts");
    // build the products of element and weight for each sample
    size_t i=0;
    RooArgList sumElements;
    RooArgList scaleElements;
    for(auto sampleit=inputFiles.begin(); sampleit!=inputFiles.end(); ++sampleit){
      // for now, we assume all the lists are nicely ordered
      TString prodname (sampleit->first.c_str());
      prodname.Append("_");
      prodname.Append(name);
//       RooProduct* prod = new RooProduct(sampleit->first.c_str(),sampleit->first.c_str(),RooArgList(*(list_weights.at(i)),*(list_phys.at(i))));
      RooProduct* prod = new RooProduct(prodname,prodname,RooArgList(*(this->_weights.at(i)),*(this->_phys.at(i))));
      sumElements.add(*prod);
      scaleElements.add(*binWidth);
      i++;
    }


    // put everything together
    DEBUG("creating RooRealSumPdf");
    RooRealSumPdf* morphfunc = new RooRealSumPdf(TString::Format("%s_morphfunc",name), name,sumElements,scaleElements,true);

    // WVE check this
    RooArgList observables(*this->_observable);
    morphfunc->addServerList(observables);
    morphfunc->addServerList(*this->_params);
    morphfunc->addOwnedComponents(this->_weights);
    morphfunc->addOwnedComponents(this->_phys);
    morphfunc->addOwnedComponents(sumElements);
    morphfunc->addOwnedComponents(scaleElements);


#ifdef USE_UBLAS    
    std::cout.precision(std::numeric_limits<double>::digits);
#endif 
#ifdef DEBUG_
    morphfunc->Print();
#endif
    DEBUG("successfully created morphing function");

    // fill the this
    this->_sumFunc = morphfunc;
  }

  static RooLagrangianMorphFunc::CacheElem* createCache(const RooLagrangianMorphFunc* func) {
    DEBUG("creating cache for basePdf" << func);
    RooLagrangianMorphFunc::ParamSet values = getParams(func->_operators);
    // create all the temporary objects required by the class

    RooLagrangianMorphFunc::CacheElem* cache = new RooLagrangianMorphFunc::CacheElem();
    TDirectory* file = openFile(func->_fileName);
    if(!file) ERROR("unable to open file '"<<func->_fileName<<"'!");
    DEBUG("reading parameter sets.");
    ParamMap inputFiles = readParamSets(file,func->_folders);
    checkNameConflict(inputFiles,func->_operators);

    DEBUG("creating basic objects");

    // Recycle existing observable, if defined
    Bool_t obsExists(kFALSE) ;
    if (func->_observables.at(0)!=0) {
      cache->_observable = (RooRealVar*)func->_observables.at(0) ;
      obsExists = kTRUE ;
    }
    cache->createBasicObjects(inputFiles,func->_obsName.c_str(),func->_operators,func->_vertices,func->_ownParameters);

    // WV register observable with class
    if (!obsExists) {
      const_cast<RooLagrangianMorphFunc*>(func)->_observables.add(*cache->_observable) ;
    }

    DEBUG("performing matrix operations");
    cache->buildMatrix(inputFiles);
    if(func->_obsName.size() == 0){
      ERROR("Matrix inversion succeeded, but no observable was supplied. quitting...");
      return cache;
    }
    DEBUG("building morphing function");

    cache->buildMorphingFunction(func->GetName(),file,inputFiles,func->_obsName.c_str(),func->_folders,func->_baseFolder.c_str());
    DEBUG("closing file");
    closeFile(file);
    setParams(values,func->_operators,true);
    return cache;
  }

  static RooLagrangianMorphFunc::CacheElem* createCache(const RooLagrangianMorphFunc* func, const Matrix& inverse) {
    // create all the temporary objects required by the class
    // function variant with precomputed inverse matrix
    DEBUG("creating cache for basePdf = " << func << " with matrix");
    RooLagrangianMorphFunc::ParamSet values = getParams(func->_operators);
    RooLagrangianMorphFunc::CacheElem* cache = new RooLagrangianMorphFunc::CacheElem();
    TDirectory* file = openFile(func->_fileName);
    if(!file) ERROR("unable to open file '"<<func->_fileName<<"'!");
    ParamMap inputFiles = readParamSets(file,func->_folders);
    checkNameConflict(inputFiles,func->_operators);
    cache->createBasicObjects(inputFiles,func->_obsName.c_str(),func->_operators,func->_vertices,func->_ownParameters);
#ifndef USE_UBLAS
    cache->_inverse.ResizeTo(inverse.GetNrows(),inverse.GetNrows());
#endif
    cache->_inverse = inverse;
    cache->_condition = NaN;
    cache->buildMorphingFunction(func->GetName(),file,inputFiles,func->_obsName.c_str(),func->_folders,func->_baseFolder.c_str());
    closeFile(file);
    setParams(values,func->_operators,true);
    return cache;
  }
};

///////////////////////////////////////////////////////////////////////////////
// Class Implementation ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// general static I/O utils
void RooLagrangianMorphFunc::writeMatrixToFile(const TMatrixD& matrix, const char* fname){
  // write a matrix to a file
  writeMatrixToFileT(matrix,fname);
}

void RooLagrangianMorphFunc::writeMatrixToStream(const TMatrixD& matrix, std::ostream& stream){
  // write a matrix to a stream
  writeMatrixToStreamT(matrix,stream);
}

TMatrixD RooLagrangianMorphFunc::readMatrixFromFile(const char* fname){
  // read a matrix from a text file
  return readMatrixFromFileT<TMatrixD>(fname);
}

TMatrixD RooLagrangianMorphFunc::readMatrixFromStream(std::istream& stream){
  // read a matrix from a stream
  return readMatrixFromStreamT<TMatrixD>(stream);
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::addFolders(const RooArgList& folders){
  // convert the RooArgList folders into a simple vector of std::string
  RooFIter folderItr = folders.fwdIterator();
  RooAbsArg* folder;
  bool foundBase = false;
  while((folder = (RooAbsArg*)(folderItr.next()))){
    RooStringVar* var = dynamic_cast<RooStringVar*>(folder);
    const std::string sample(var ? var->getVal() : var->GetName());
    if(sample.size() == 0) continue;
    DEBUG("adding sample: '" << sample << "'");
    this->_folders.push_back(sample);
    if(sample == this->_baseFolder){
      foundBase = true;
    }
  }
  if(this->_folders.size() > 0){
    if(!foundBase){
      if(this->_baseFolder.size() > 0){
        this->_folders.insert(this->_folders.begin(),this->_baseFolder);
      } else {
        this->_baseFolder= _folders[0];
      }
    }
  } else {
    TDirectory* file = openFile(this->_fileName.c_str());
    TIter next(file->GetList());
    TObject *obj = NULL;
    while ((obj = (TObject*)next())) {
      TFolder * f = dynamic_cast<TFolder*>(file->Get(obj->GetName()));
      if(!f) continue;
      std::string name(f->GetName());
      if(name.size() == 0) continue;
      // was 
      //if(this->_baseFolder.size() == 0) this->_baseFolder == name;
      // GB I presume this was intented ?
      if(this->_baseFolder.size() == 0) this->_baseFolder = name;
      if(this->_baseFolder == name){
        this->_folders.insert(this->_folders.begin(),name);
      } else {
        this->_folders.push_back(name);
      }
    }
    closeFile(file);
  }
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title,
                                 const char* fileName, const char* obsName, const char* basefolder, const RooArgList& folders) :
  RooAbsPdf(name,title),
  _cacheMgr(this,10,kTRUE,kTRUE),
  _fileName(fileName),
  _obsName(obsName),
  _baseFolder(basefolder),
  _operators("operators","set of operators",this,kTRUE,kFALSE),
  _observables("observables","set of observables",this,kTRUE,kFALSE),
  _curNormSet(0)
{
  // protected constructor with proper arguments
  this->printAuthors();
  this->addFolders(folders);
}

//_____________________________________________________________________________
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title,
                                 const char* fileName, const char* obsName, const RooArgList& folders) :
  RooLagrangianMorphFunc(name,title,fileName,obsName,"",folders)
{
  // protected constructor with proper arguments
}

//_____________________________________________________________________________
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title,
                                 const char* fileName, const char* obsName) :
  RooLagrangianMorphFunc(name,title,fileName,obsName,"",RooArgList())
{
  // protected constructor with proper arguments
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings, const char* basefolder, const RooArgList& folders) :
  RooLagrangianMorphFunc(name,title,fileName,obsName,basefolder,folders)
{
  // constructor with proper arguments
  RooArgList operators;
  extractOperators(prodCouplings,operators);
  extractOperators(decCouplings,operators);
  this->setup(operators,prodCouplings,decCouplings, false);
}

//_____________________________________________________________________________
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings, const RooArgList& folders) :
  RooLagrangianMorphFunc(name,title,fileName,obsName,folders)
{
  // constructor with proper arguments
  RooArgList operators;
  extractOperators(prodCouplings,operators);
  extractOperators(decCouplings,operators);
  this->setup(operators,prodCouplings,decCouplings, false);
}

//_____________________________________________________________________________
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings) :
  RooLagrangianMorphFunc(name,title,fileName,obsName)
{
  // constructor with proper arguments
  RooArgList operators;
  extractOperators(prodCouplings,operators);
  extractOperators(decCouplings,operators);
  this->setup(operators,prodCouplings, decCouplings, false);
}

//_____________________________________________________________________________

template<class T>
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const std::vector<T>& vertices, const char* basefolder, const RooArgList& folders) :
  RooLagrangianMorphFunc(name,title,fileName,obsName,basefolder,folders)
{
  // constructor with proper arguments
  RooArgList operators;
  for(const auto& v:vertices){
    extractOperators(v,operators);
  }
  this->setup(operators,vertices, false);
}

// Are the lines below needed ? 
//template RooLagrangianMorphFunc::RooLagrangianMorphFunc<RooArgList>(const char *, const char *, const char*, const char*, const std::vector<RooArgList>&, const char*, const RooArgList&){};
//template RooLagrangianMorphFunc::RooLagrangianMorphFunc<RooArgSet> (const char *, const char *, const char*, const char*, const std::vector<RooArgSet> &, const char*, const RooArgList&){};

//_____________________________________________________________________________
template<class T>
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const std::vector<T>& vertices, const RooArgList& folders) :
  RooLagrangianMorphFunc(name,title,fileName,obsName,folders)
{
  // constructor with proper arguments
  RooArgList operators;
  for(const auto& v:vertices){
    extractOperators(v,operators);
  }
  this->setup(operators,vertices, false);
}

// Are the lines below needed ? 
//template RooLagrangianMorphFunc::RooLagrangianMorphFunc<RooArgList>(const char *, const char *, const char*, const char*, const std::vector<RooArgList>&, const RooArgList&){}
//template RooLagrangianMorphFunc::RooLagrangianMorphFunc<RooArgSet> (const char *, const char *, const char*, const char*, const std::vector<RooArgSet> &, const RooArgList&){}

//_____________________________________________________________________________
template<class T>
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const char *name, const char *title, const char* fileName, const char* obsName, const std::vector<T>& vertices) :
  RooLagrangianMorphFunc(name,title,fileName,obsName)
{
  // constructor with proper arguments
  RooArgList operators;
  for(const auto& v:vertices){
    extractOperators(v,operators);
  }
  this->setup(operators,vertices, false);
}

// Are the lines below needed ? 
//template RooLagrangianMorphFunc::RooLagrangianMorphFunc<RooArgList>(const char *, const char *, const char*, const char*, const std::vector<RooArgList>&){};
//template RooLagrangianMorphFunc::RooLagrangianMorphFunc<RooArgSet> (const char *, const char *, const char*, const char*, const std::vector<RooArgSet> &){};

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setup(const RooArgSet& operators, const RooAbsCollection& prodCouplings, const RooAbsCollection& decCouplings, bool own){
  // setup this instance with the given set of operators and vertices
  // if own=true, the class will own the operators
  DEBUG("setup(ops,prod,decay,"<<own<<") called");
  this->_ownParameters = own;
  this->_vertices.push_back(new RooListProxy("!production","set of couplings in the production vertex",this,kTRUE,kFALSE));
  this->_vertices.push_back(new RooListProxy("!decay",     "set of couplings in the decay vertex",     this,kTRUE,kFALSE));
  if(own){
    DEBUG("adding own operators");
    this->_operators.addOwned(operators);
    this->_vertices[0]->addOwned(prodCouplings);
    this->_vertices[1]->addOwned(decCouplings);
  } else {
    DEBUG("adding non-own operators");
    this->_operators.add(operators);
    this->_vertices[0]->add(prodCouplings);
    this->_vertices[1]->add(decCouplings);
  }
}

//_____________________________________________________________________________

template<class T>
void RooLagrangianMorphFunc::setup(const RooArgSet& operators, const std::vector<T>& vertices, bool own){
  // setup this instance with the given set of operators and vertices
  // if own=true, the class will own the operators
  this->_ownParameters = own;
  if(own){
    this->_operators.addOwned(operators);
  } else {
    this->_operators.add(operators);
  }
  for(size_t i=0; i<vertices.size(); i++){
    std::stringstream name;
    name << "!vertex" << i;
    std::stringstream title;
    title << "set of couplings in the vertex " << i;
    this->_vertices.push_back(new RooListProxy(name.str().c_str(),title.str().c_str(),this,kTRUE,kFALSE));
    if(own){
      this->_vertices[i]->addOwned(vertices[i]);
    } else {
      this->_vertices[i]->add(vertices[i]);
    }
  }
}

//_____________________________________________________________________________
RooLagrangianMorphFunc::RooLagrangianMorphFunc(const RooLagrangianMorphFunc& other, const char* name) :
  RooAbsPdf(other,name),
  _cacheMgr(other._cacheMgr,this),
  _fileName(other._fileName),
  _obsName(other._obsName),
  _baseFolder(other._baseFolder),
  _folders(other._folders),
  _operators(other._operators.GetName(),this,other._operators),
  _observables(other._observables.GetName(),this,other._observables),
  _curNormSet(0)
{
  // copy constructor
  DEBUG("copy constructor called");
  for(size_t i=0; i<other._vertices.size(); ++i){
    this->_vertices.push_back(new RooListProxy(other._vertices[i]->GetName(),this,*(other._vertices[i])));
  }
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::RooLagrangianMorphFunc() :
  _observables("observables","set of observables",this,kTRUE,kFALSE)
{
  // default constructor
  DEBUG("default constructor called");
  this->printAuthors();
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::~RooLagrangianMorphFunc() {
  // default destructor
  DEBUG("destructor called");
  for(auto v:this->_vertices){
    delete v;
  }
}

//_____________________________________________________________________________

TObject* RooLagrangianMorphFunc::clone(const char* newname) const {
  // cloning method
  return new RooLagrangianMorphFunc(*this,newname);
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::printAuthors() const {
  // print the author information
  std::cout << "\033[1mRooLagrangianMorphFunc\033[0m: a RooStats class for morphing physics distributions between configurations. authors:" << std::endl;
  std::cout << "   " << "Lydia Brenner   (lbrenner@cern.ch)" << std::endl;
  std::cout << "   " << "Carsten Burgard (cburgard@cern.ch)" << std::endl;
  std::cout << "   " << "Katharina Ecker (kecker@cern.ch)" << std::endl;
  std::cout << "   " << "Adam Kaluza     (akaluza@cern.ch)" << std::endl;
  std::cout << "please feel free to contact with questions and suggestions." << std::endl;
}

//_____________________________________________________________________________

int RooLagrangianMorphFunc::countSamples(int nprod, int ndec, int nboth){
  // calculate the number of samples needed to morph a bivertex, 2-2 physics process
  VertexMap vertexmap;
  std::vector<bool> prod;
  std::vector<bool> dec;
  for(int i=0; i<nboth; ++i){
    prod.push_back(true);
    dec.push_back(true);
  }
  for(int i=0; i<nprod; ++i){
    prod.push_back(true);
    dec.push_back(false);
  }
  for(int i=0; i<ndec; ++i){
    prod.push_back(false);
    dec.push_back(true);
  }
  vertexmap.push_back(prod);
  vertexmap.push_back(dec);
  MorphFuncPattern morphfuncpattern(calculateFunction(vertexmap));
  return morphfuncpattern.size();
}

//_____________________________________________________________________________

int RooLagrangianMorphFunc::countSamples(std::vector<RooArgList*>& vertices){
  // calculate the number of samples needed to morph a certain physics process
  // usage:
  //   countSamples ( { RooLagrangianMorphFunc::makeHCggfCouplings(), RooLagrangianMorphFunc::makeHCHZZCouplings() } )

  RooArgList operators,couplings;
  for(auto vertex: vertices){
    extractOperators(*vertex,operators);
    extractCouplings(*vertex,couplings);
  }
  VertexMap vertexmap(buildVertexMap<RooArgList>(vertices,couplings));
  MorphFuncPattern morphfuncpattern(calculateFunction(vertexmap));
  return morphfuncpattern.size();
}

//_____________________________________________________________________________

TPair* RooLagrangianMorphFunc::makeCrosssectionContainer(double xs, double unc){
  TPair* v = new TPair(new TParameter<double>("xsection",xs),new TParameter<double>("uncertainty",unc));
  return v;
}

//_____________________________________________________________________________

RooParamHistFunc* RooLagrangianMorphFunc::getBaseTemplate(){
  // find the one component that is a ParamHistFunc
  RooRealSumPdf* mf = dynamic_cast<RooRealSumPdf*>(this->getPdf());
  if(!mf) ERROR("unable to retrieve morphing function");
  RooArgSet* args = mf->getComponents();
  RooFIter itr(args->fwdIterator());
  TObject* obj;
  while((obj = itr.next())){
    RooProduct* prod = dynamic_cast<RooProduct*>(obj);
    RooFIter subitr(prod->components().fwdIterator());
    TObject* subobj;
    while((subobj = itr.next())){
      RooParamHistFunc* p = dynamic_cast<RooParamHistFunc*>(obj);
      if(p){
        return p;
      }
    }
  }
  return NULL;
}

//_____________________________________________________________________________

RooProduct* RooLagrangianMorphFunc::getSumElement(const char* name) const {
  // return the RooProduct that is the element of the RooRealSumPdf corresponding to the given sample name
  RooRealSumPdf* mf = dynamic_cast<RooRealSumPdf*>(this->getPdf());
  if(!mf) ERROR("unable to retrieve morphing function");
  RooArgSet* args = mf->getComponents();
  RooFIter itr(args->fwdIterator());
  TObject* obj;
  TString prodname (name);
  prodname.Append("_");
  prodname.Append(this->GetName());
  while((obj = itr.next())){
    RooProduct* prod = dynamic_cast<RooProduct*>(obj);
    if(!prod) continue;
    TString sname(prod->GetName());
    if(sname.CompareTo(prodname) == 0){
      return prod;
    }
  }
  return NULL;
}

//_____________________________________________________________________________

const std::vector<std::string>& RooLagrangianMorphFunc::getSamples() const {
  // return the vector of sample names, used to build the morph func
  return this->_folders;
}

//_____________________________________________________________________________

RooAbsReal* RooLagrangianMorphFunc::getSampleWeight(const char* name){
  // retrieve the weight (prefactor) of a sample with the given name
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  TString wname(name);
  wname.Prepend("w_");
  wname.Append("_");
  wname.Append(this->GetName());
  return dynamic_cast<RooAbsReal*>(cache->_weights.find(wname));
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::printSampleWeights(){
  // print the current sample weights
  RooRealSumPdf* mf = dynamic_cast<RooRealSumPdf*>(this->getPdf());
  if(!mf) ERROR("unable to retrieve morphing function");
  RooFIter itr(mf->getComponents()->fwdIterator());
  TObject* w = NULL;
  mf->getComponents()->Print();
  while((w = itr.next())){
    TString wname(w->GetName());
    if(!wname.BeginsWith("w_")) continue;
    RooFormulaVar* weight = dynamic_cast<RooFormulaVar*>(w);
    if(!weight) continue;
    std::cout << weight->GetName() << " = " << weight->getVal() << std::endl;
  }
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::randomizeParameters(double z){
  // randomize the parameters a bit
  // useful to test and debug fitting
  RooFIter itr(_operators.fwdIterator());
  RooRealVar* obj;
  TRandom3 r;
  while((obj = dynamic_cast<RooRealVar*>(itr.next()))){
    double val = obj->getVal();
    if(obj->isConstant()) continue;
    double variation = r.Gaus(1,z);
    obj->setVal(val*variation);
  }
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::updateCoefficients(){
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  TDirectory* file = openFile(this->_fileName);
  if(!file){
    ERROR("unable to open file '"<<this->_fileName<<"'!");
    return false;
  }
  DEBUG("reading parameter sets.");
  ParamMap inputFiles = readParamSets(file,this->_folders);
  cache->buildMatrix(inputFiles);
  cache->updateMorphingFunction(this,file,inputFiles,this->_obsName.c_str(),this->_folders,this->_baseFolder.c_str());
  closeFile(file);
  return true;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::useCoefficients(const TMatrixD& inverse){
  // setup the morphing function with a predefined inverse matrix
  // call this function *before* any other after creating the object
  RooLagrangianMorphFunc::CacheElem* cache = (RooLagrangianMorphFunc::CacheElem*) _cacheMgr.getObj(0,(RooArgSet*)0);
  Matrix m = makeSuperMatrix(inverse);
  if (cache) {
#ifdef USE_UBLAS
    cache->_inverse = m;
    TDirectory* file = openFile(this->_fileName);
    if(!file) ERROR("unable to open file '"<<this->_fileName<<"'!");
    DEBUG("reading parameter sets.");
    ParamMap inputFiles = readParamSets(file,this->_folders);
    cache->updateMorphingFunction(this,file,inputFiles,this->_obsName.c_str(),this->_folders,this->_baseFolder.c_str());
#else
    return false;
#endif
  } else {
    cache = RooLagrangianMorphFunc::CacheElem::createCache(this,m);
    if(!cache) ERROR("unable to create cache!");
    this->_cacheMgr.setObj(0,0,cache,0) ;
  }
  return true;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::useCoefficients(const char* filename){
  // setup the morphing function with a predefined inverse matrix
  // call this function *before* any other after creating the object
  RooLagrangianMorphFunc::CacheElem* cache = (RooLagrangianMorphFunc::CacheElem*) _cacheMgr.getObj(0,(RooArgSet*)0);
  if (cache) {
    return false;
  }
  cache = RooLagrangianMorphFunc::CacheElem::createCache(this,readMatrixFromFileT<Matrix>(filename));
  if(!cache) ERROR("unable to create cache!");
  this->_cacheMgr.setObj(0,0,cache,0);
  return true;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::writeCoefficients(const char* filename){
  // write the inverse matrix to a file
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  if(!cache) return false;
  writeMatrixToFileT(cache->_inverse,filename);
  return true;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::writeFormulas(const char* filename){
    // write a matrix to a text file
  std::ofstream stream(filename);
  if(!stream.good()){
    ERROR("unable to read file '"<<filename<<"'!");
    return false;
  }
  
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  for(auto const & formula : cache->_formulas){
    stream << (formula.second)->GetTitle() << std::endl;
  }
  stream.close();  
  return true;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::writePhysics(const char* filename){
      // write a matrix to a text file
  std::ofstream stream(filename);
  if(!stream.good()){
    ERROR("unable to read file '"<<filename<<"'!");
    return false;
  }
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  for(int i=0; i<cache->_phys.getSize(); ++i){
    RooProduct* phys = (RooProduct*) cache->_phys.at(i);
    if(!phys) continue;
    stream << "# " << phys->GetName() << std::endl;
    RooAbsReal* integral = phys->createIntegral(*this->getObservable());
    if(integral) stream << integral->getVal() << std::endl;
  }
  stream.close();  
  return true;
}
//_____________________________________________________________________________

RooLagrangianMorphFunc::CacheElem* RooLagrangianMorphFunc::getCache(const RooArgSet* /*nset*/) const {
  // retrieve the cache object
  RooLagrangianMorphFunc::CacheElem* cache = (RooLagrangianMorphFunc::CacheElem*) _cacheMgr.getObj(0,(RooArgSet*)0);
  if (!cache) {
    cache = RooLagrangianMorphFunc::CacheElem::createCache(this);
    if(cache) this->_cacheMgr.setObj(0,0,cache,0);
    else ERROR("unable to create cache!");
  }
  return cache;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::hasCache() const {
  // return true if a cache object is present, false otherwise
  return (bool)(_cacheMgr.getObj(0,(RooArgSet*)0));
}

//_____________________________________________________________________________
RooArgList RooLagrangianMorphFunc::CacheElem::containedArgs(Action)
{
    // RooRealSumPdf* _sumFunc;
    // RooArgList* _params;
    // RooArgList* _couplings;
    // RooRealVar* _observable;

    // FormulaList _formulas;
    // RooArgList _phys;
    // RooArgList _weights;

  // retrieve the list of contained args
  return RooArgList(*_sumFunc,*_observable);
}

//_____________________________________________________________________________
RooLagrangianMorphFunc::CacheElem::~CacheElem()
{
  // default destructor
  delete _sumFunc; // the sumfunc owns all its contents
  if (_ownObs) {
    delete _observable; // the observable is a single object
  }
  delete _params; // the param list does _NOT_ own its contents
  delete _couplings; // the coupling list does _NOT_ own its contents
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameter(const char* name, double value){
  // set one parameter to a specific value
  RooRealVar* param = this->getParameter(name);
  if(!param){
    return;
  }
  if(value > param->getMax()) param->setMax(value);
  if(value < param->getMin()) param->setMin(value);
  param->setVal(value);
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameter(const char* name, double value, double min, double max){
  // set one parameter to a specific value and range
  RooRealVar* param = this->getParameter(name);
  if(!param){
    return;
  }
  param->setMin(min);
  param->setMax(max);
  param->setVal(value);
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameter(const char* name, double value, double min, double max, double error){
  // set one parameter to a specific value and range
  RooRealVar* param = this->getParameter(name);
  if(!param){
    return;
  }
  param->setMin(min);
  param->setMax(max);
  param->setVal(value);
  param->setError(error);
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::isParameterConstant(const char* name) const{
  // return true if the parameter with the given name is set constant, false otherwise
  RooRealVar* param = this->getParameter(name);
  if(param){
    return param->isConstant();
  }
  return true;
}

//_____________________________________________________________________________

RooRealVar* RooLagrangianMorphFunc::getParameter(const char* name) const{
  // retrieve the RooRealVar object incorporating the parameter with the given name
  RooAbsCollection* args = this->getParameterSet();
  if(!args){
    return NULL;
  }
  RooRealVar* param = dynamic_cast<RooRealVar*>(args->find(name));
  if(!param){
    return NULL;
  }
  return param;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::hasParameter(const char* name) const{
  // check if a parameter of the given name is contained in the list of known parameters
  RooRealVar* p = this->getParameter(name);
  if(p) return true;
  return false;
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameterConstant(const char* name, bool constant) const {
  // call setConstant with the boolean argument provided on the parameter with the given name
  RooRealVar* param = this->getParameter(name);
  if(param){
    return param->setConstant(constant);
  }
}

//_____________________________________________________________________________

double RooLagrangianMorphFunc::getParameterValue(const char* name) const{
  // set one parameter to a specific value
  RooRealVar* param = this->getParameter(name);
  if(param){
    return param->getVal();
  }
  return 0;
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameters(TH1* paramhist){
  // set the morphing parameters to those supplied in the given param hist
  setParams(paramhist,*(this->getParameterSet()),false);
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameters(const char* foldername){
  // set the morphing parameters to those supplied in the sample with the given name
  TDirectory* file = openFile(this->_fileName);
  TH1* paramhist = getParamHist(file,foldername);
  setParams(paramhist,*(this->getParameterSet()),false);
  closeFile(file);
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::ParamSet RooLagrangianMorphFunc::getParameters(const char* foldername) const {
  // retrieve the morphing parameters associated to the sample with the given name
  TDirectory* file = openFile(this->_fileName);
  RooLagrangianMorphFunc::ParamSet point = readParamSet(file,foldername);
  closeFile(file);
  return point;
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameters(const RooArgList* list){
  // set the morphing parameters to those supplied in the list with the given name
    RooFIter itr(list->fwdIterator());
    TObject* obj;
    while((obj = itr.next())){
      RooRealVar* param = dynamic_cast<RooRealVar*>(obj);
      if(!param) continue;
      this->setParameter(param->GetName(),param->getVal());

    }
}

//_____________________________________________________________________________

TH1* RooLagrangianMorphFunc::createTH1(const std::string& name, RooFitResult* r){
  // retrieve a histogram output of the current morphing settings
  return this->createTH1(name,false,r);
}

TH1* RooLagrangianMorphFunc::createTH1(const std::string& name, bool correlateErrors, RooFitResult* r){
  // retrieve a histogram output of the current morphing settings
  RooRealSumPdf* pdf = this->getPdf();
  RooRealVar* observable = this->getObservable();
//   var.getBinning().Print()
  
  const int nbins = observable->getBins();
  
  TH1* hist = new TH1F(name.c_str(),name.c_str(),nbins,observable->getBinning().array());
  
//   TH1* hist = new TH1F(name.c_str(),name.c_str(),nbins,observable->getMin(),observable->getMax());
  bool ownResult = !(bool)(r);
  // if(!r){
  //   r = RooFitResult::prefitResult(*(this->getParameterSet()));
  // }
  RooArgSet* args = pdf->getComponents();
  TObject* obj;
  for (int i=0; i<nbins; ++i) {
    observable->setBin(i);
    RooFIter itr(args->fwdIterator());
    double val = 0;
    double unc2 = 0;
    double unc = 0;
    while((obj = itr.next())){
      RooProduct* prod = dynamic_cast<RooProduct*>(obj);
      if(!prod) continue;
      RooAbsArg* phys = prod->components().find(TString::Format("phys_%s",prod->GetName()));
      RooHistFunc* hf = dynamic_cast<RooHistFunc*>(phys);
      const RooDataHist& dhist = hf->dataHist();
      dhist.get(i);
      RooAbsReal* formula = dynamic_cast<RooAbsReal*>(prod->components().find(TString::Format("w_%s",prod->GetName())));
      double weight = formula->getVal();
      unc2 += dhist.weightSquared()*weight*weight;
      unc += sqrt(dhist.weightSquared())*weight;
      val += dhist.weight()*weight;
    }
    hist->SetBinContent(i+1,val);
    hist->SetBinError(i+1,correlateErrors ? unc : sqrt(unc2));
  }
  if(ownResult) delete r;
  return hist;
}


//_____________________________________________________________________________

int RooLagrangianMorphFunc::countContributingFormulas() const{
  // count the number of formulas that correspond to the current parameter set
  int nFormulas = 0;
  RooRealSumPdf* mf = dynamic_cast<RooRealSumPdf*>(this->getPdf());
  if(!mf) ERROR("unable to retrieve morphing function");
  RooArgSet* args = mf->getComponents();
  RooFIter itr(args->fwdIterator());
  TObject* obj;
  while((obj = itr.next())){
    RooProduct* prod = dynamic_cast<RooProduct*>(obj);
    if(prod->getVal() != 0){
      nFormulas++;
    }
  }
  return nFormulas;
}

//_____________________________________________________________________________

bool RooLagrangianMorphFunc::isParameterUsed(const char* paramname) const {
  // check if there is any morphing power provided for the given parameter
  // morphing power is provided as soon as any two samples provide different, non-zero values for this parameter
  std::string pname(paramname);
  ParamMap inputFiles = readParamSets(this->_fileName,this->_folders);
  double val = 0;
  bool isUsed = false;
  for(auto sample : inputFiles){
    double thisval = sample.second[pname];
    if(thisval != val){
      if(val != 0) isUsed = true;
      val = thisval;
    }
  }
  return isUsed;
}




//_____________________________________________________________________________

bool RooLagrangianMorphFunc::isCouplingUsed(const char* couplname) const {
  // check if there is any morphing power provided for the given coupling
  // morphing power is provided as soon as any two samples provide different, non-zero values for this coupling
  std::string cname(couplname);
  RooArgList* args = this->getCouplingSet();
  RooAbsReal* coupling = dynamic_cast<RooAbsReal*>(args->find(couplname));
  if(!coupling) return false;
  RooLagrangianMorphFunc::ParamSet params = this->getParameters();
  ParamMap inputFiles = readParamSets(this->_fileName,this->_folders);
  RooLagrangianMorphFunc* self = const_cast<RooLagrangianMorphFunc*>(this);
  double val = 0;
  bool isUsed = false;
  for(auto sample : inputFiles){
    self->setParameters(sample.second);
    double thisval = coupling->getVal();
    if(thisval != val){
      if(val != 0) isUsed = true;
      val = thisval;
    }
  }
  self->setParameters(params);
  return isUsed;
}


//_____________________________________________________________________________

void RooLagrangianMorphFunc::printParameters(const char* samplename) const {
  // print all the parameters and their values in the given sample to the console
  RooLagrangianMorphFunc::ParamSet params = readParamSet(this->_fileName,samplename);
  for(auto param : params){
    if(this->hasParameter(param.first.c_str())){
      std::cout << param.first << " = " << param.second;
      if(this->isParameterConstant(param.first.c_str())) std::cout << " (const)";
      std::cout << std::endl;
    }
  }
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::printSamples() const {
  // print all the known samples to the console
  for(auto folder : this->_folders){
    std::cout << folder << std::endl;
    if(folder ==  this->_baseFolder) std::cout << "*" << std::endl;
  }
}

//_____________________________________________________________________________
void RooLagrangianMorphFunc::printFormulas() const {
  // print the morphing formulas
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  for(auto const & formula : cache->_formulas){
    std::cout << formula.first << " ";
    (formula.second)->Print();
    std::cout << std::endl;
  }
}

//_____________________________________________________________________________
void RooLagrangianMorphFunc::printPhysics() const {
  // print the current physics values
  RooLagrangianMorphFunc::CacheElem* cache = this->getCache(_curNormSet);
  for(int i=0; i<cache->_phys.getSize(); ++i){
    RooAbsArg* phys = cache->_phys.at(i);
    if(!phys) continue;
    phys->Print();
  }
}

//_____________________________________________________________________________

int RooLagrangianMorphFunc::nParameters() const {
  // return the number of parameters in this morphing function
  return this->getParameterSet()->getSize();
}

//_____________________________________________________________________________

int RooLagrangianMorphFunc::nSamples() const {
  // return the number of samples in this morphing function
  return this->_folders.size();
}

//_____________________________________________________________________________

int RooLagrangianMorphFunc::nPolynomials() const {
  // return the number of samples in this morphing function
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  return cache->_formulas.size();
}


//_____________________________________________________________________________

void RooLagrangianMorphFunc::printEvaluation() const {
  // print the contributing samples and their respective weights
  RooRealSumPdf* mf = dynamic_cast<RooRealSumPdf*>(this->getPdf());
  if(!mf){
    std::cerr << "Error: unable to retrieve morphing function" << std::endl;
    return;
  }
  RooArgSet* args = mf->getComponents();
  RooFIter itr(args->fwdIterator());
  TObject* obj;
  while((obj = itr.next())){
    RooFormulaVar* formula = dynamic_cast<RooFormulaVar*>(obj);
    if(formula){
      TString name(formula->GetName());
      name.Remove(0,2);
      name.Prepend("phys_");
      if(!args->find(name.Data())){
        continue;
      }
      double val = formula->getVal();
      if(val != 0){
        std::cout << formula->GetName() << ": " << val << " = " << formula->GetTitle() << std::endl;
      }
    }
  }
}

//_____________________________________________________________________________

RooArgList* RooLagrangianMorphFunc::getParameterSet() const {
  // get the set of parameters
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  return cache->_params;
}

//_____________________________________________________________________________

RooArgList* RooLagrangianMorphFunc::getCouplingSet() const {
  // get the set of couplings
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  return cache->_couplings;
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::ParamSet RooLagrangianMorphFunc::getCouplings() const {
  // retrieve a set of couplings
  RooFIter itr(this->getCouplingSet()->fwdIterator());
  RooAbsArg* obj;
  RooLagrangianMorphFunc::ParamSet couplings;
  while((obj = itr.next())){
    RooAbsReal* var = dynamic_cast<RooAbsReal*>(obj);
    if(!var) continue;
    const std::string name(var->GetName());
    double val = var->getVal();
    couplings[name] = val;
  }
  return couplings;
}

//_____________________________________________________________________________

RooLagrangianMorphFunc::ParamSet RooLagrangianMorphFunc::getParameters() const {
  // retrieve a set of couplings
  return getParams(*(this->getParameterSet()));
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::setParameters(const ParamSet& params) {
  // retrieve a set of couplings
  setParams(params,*(this->getParameterSet()),false);
}

//_____________________________________________________________________________

RooRealVar* RooLagrangianMorphFunc::getObservable() const {
  // get the observable
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  return cache->_observable;
}

//_____________________________________________________________________________

RooRealSumPdf* RooLagrangianMorphFunc::getPdf() const {
  // get the pdf
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  return cache->_sumFunc;
}

//_____________________________________________________________________________

RooRealSumPdf* RooLagrangianMorphFunc::clonePdf() const {
  // get a standalone clone of the pdf that does not depend on this object
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  RooRealSumPdf* orig = cache->_sumFunc;
  RooRealSumPdf* pdfclone = new RooRealSumPdf(orig->GetName(),orig->GetTitle(),orig->funcList(),orig->coefList(),"true");
  return pdfclone;
}

//_____________________________________________________________________________

RooAbsPdf::ExtendMode RooLagrangianMorphFunc::extendMode() const {
  // Return extended mode capabilities
  return this->getPdf()->extendMode();
}

//_____________________________________________________________________________

Double_t RooLagrangianMorphFunc::expectedEvents(const RooArgSet* nset) const {
    // Return expected number of events for extended likelihood calculation
    // which is the sum of all coefficients
  return this->getPdf()->expectedEvents(nset);
}

//_____________________________________________________________________________

Double_t RooLagrangianMorphFunc::expectedEvents() const {
  // return the number of expected events for the current parameter set
  RooArgSet set;
  set.add(*this->getObservable());
  return this->getPdf()->expectedEvents(set);
}

//_____________________________________________________________________________

Double_t RooLagrangianMorphFunc::expectedEvents(const RooArgSet& nset) const {
  // Return expected number of events for extended likelihood calculation
  // which is the sum of all coefficients
  return getPdf()->expectedEvents(&nset) ;
}

//_____________________________________________________________________________

double RooLagrangianMorphFunc::expectedUncertainty() const {
  // return the expected uncertainty for the current parameter set
  RooRealVar* observable = this->getObservable();
  CacheElem* cache = this->getCache(_curNormSet);
  double unc2 = 0;
  for(int i=0; i<cache->_phys.getSize(); ++i){
    RooAbsArg* phys = cache->_phys.at(i);
    RooAbsReal* weight = (RooAbsReal*)(cache->_weights.at(i));
    double newunc2 = 0;
    RooHistFunc* hf = dynamic_cast<RooHistFunc*>(phys);
    RooRealVar* rv = dynamic_cast<RooRealVar*>(phys);
    if(hf){
      const RooDataHist& hist = hf->dataHist();
      for(Int_t j=0; j<observable->getBins(); ++j){
        hist.get(j);
        newunc2 += hist.weightSquared();
      }
    } else if(rv){
      newunc2 = pow(rv->getError(),2);
    }
    double w = weight->getVal();
    unc2 += newunc2*w*w;
    // std::cout << phys->GetName() << " : " << weight->GetName() << " thisweight: " <<  w << " thisxsec2: " << newunc2 << " weight " << weight << std::endl;
  }
  return sqrt(unc2);
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::printParameters() const {
  // print the parameters and their current values
  RooArgList* parameters = this->getParameterSet();
  RooFIter itr(parameters->fwdIterator());
  TObject* obj;
  while((obj = itr.next())){
    RooRealVar* param = dynamic_cast<RooRealVar*>(obj);
    if(!param) continue;
    std::cout << param->GetName() << ": " << param->getVal();
    if(param->isConstant()) std::cout << " (const)";
    else {
      std::cout << " +" << param->getAsymErrorHi() << " -" << param->getAsymErrorLo();
      std::cout << " (" << param->getMin() << " - " << param->getMax() << ")";
    }
    std::cout << std::endl;
  }
}

//_____________________________________________________________________________

void RooLagrangianMorphFunc::printCouplings() const {
  // print a set of couplings
  RooLagrangianMorphFunc::ParamSet couplings = this->getCouplings();
  for(auto c : couplings){
    std::cout << c.first << ": " << c.second << std::endl;
  }
}

//_____________________________________________________________________________

std::list<Double_t>* RooLagrangianMorphFunc::binBoundaries(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  // retrieve the list of bin boundaries
  return this->getPdf()->binBoundaries(obs,xlo,xhi);
}

//_____________________________________________________________________________

std::list<Double_t>* RooLagrangianMorphFunc::plotSamplingHint(RooAbsRealLValue& obs, Double_t xlo, Double_t xhi) const {
  // retrieve the sampling hint
  return this->getPdf()->plotSamplingHint(obs,xlo,xhi);
}

//_____________________________________________________________________________

Double_t RooLagrangianMorphFunc::getValV(const RooArgSet* set) const
{
  //cout << "XX RooLagrangianMorphFunc::getValV(" << this << ") set = " << set << endl ;
  _curNormSet = set ;
  return RooAbsPdf::getValV(set) ;
}

//_____________________________________________________________________________

Bool_t RooLagrangianMorphFunc::selfNormalized() const {
  // P.d.f is self normalized
  return kTRUE ;
}

//_____________________________________________________________________________

Double_t RooLagrangianMorphFunc::evaluate() const {
  // call getVal on the Pdf
  RooAbsPdf* pdf = this->getPdf();
  if(pdf) return pdf->getVal(_curNormSet);
  else ERROR("unable to aquire in-built pdf!");
  return 0.;
}

//_____________________________________________________________________________

Bool_t  RooLagrangianMorphFunc::isBinnedDistribution(const RooArgSet& obs) const {
  // check if this PDF is a binned distribution in the given observable
  return this->getPdf()->isBinnedDistribution(obs);
}

//_____________________________________________________________________________

TMatrixD RooLagrangianMorphFunc::getMatrix() const {
  // retrieve the matrix of coefficients before inversion
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  if(!cache) ERROR("unable to retrieve cache!");
  return makeRootMatrix(cache->_matrix);
}

//_____________________________________________________________________________

TMatrixD RooLagrangianMorphFunc::getInvertedMatrix() const {
  // retrieve the matrix of coefficients after inversion
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  if(!cache) ERROR("unable to retrieve cache!");
  return makeRootMatrix(cache->_inverse);
}

//_____________________________________________________________________________


double RooLagrangianMorphFunc::getCondition() const {
  // retrieve the condition of the coefficient matrix
  RooLagrangianMorphFunc::CacheElem* cache = getCache(_curNormSet);
  if(!cache) ERROR("unable to retrieve cache!");
  return cache->_condition;
}
