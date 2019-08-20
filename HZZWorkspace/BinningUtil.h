#ifndef _HZZWS_BINNINGUTIL_H_
#define _HZZWS_BINNINGUTIL_H_
#include "TString.h"
class TTree;
class TH2F;
class TH2D;
class TH2;
class TF1;
class RooBinning;
class RooArgList;
class RooDataHist;
class RooCurve;
class RooAbsPdf;
class RooArgSet;

#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooBinning.h"
#include "RooFitResult.h"
#include <vector>

namespace BinningUtil 
{
  void getBinning( RooBinning& binning, RooRealVar& x, RooAbsReal& func );
  double nextstep(TF1* wf, double x0, double x1);
  
  void randomizeBinned( RooDataHist& myCopy );

  std::vector<RooDataHist> dataVariations(RooArgList& parList, const RooDataHist& input, const double& scale=1.0, const double& cutoff=-1.0 ); 
  std::vector<RooDataHist> histVariations(RooRealVar& x, const std::vector<RooDataHist>& inputVec, const double& rho=1.0, 
                                          const TString& options="am", const Int_t& nSigma=3, 
                                          const Bool_t& rotation=kTRUE, const Bool_t& sortInput=kTRUE);

  void pdfVariations(RooArgList&, RooRealVar& x, const std::vector<RooDataHist>& inputVec);

  void createPdf( RooArgList& parList, RooArgList& pdfList, RooRealVar& x, const RooDataHist& input, 
		  const TString& options, const double& rho, const Int_t& nSigma, const Bool_t& rotation, const Bool_t& sortInput,
		  const double& scale, const double& cutoff );

  RooCurve* plotOnWithErrorBand( RooPlot* frame, Double_t Z, RooAbsPdf& pdf, RooArgSet& obs, const int& nCurves=-1, const int& nEvents=-1, 
				 const char* binning=0, RooDataHist* dataHist=0 ); 

  RooDataHist* makeAsimov1D( const RooAbsPdf& pdf, RooRealVar& obs, const RooBinning& binning, const char* binningName, TH1* ehist=0 );

  RooDataHist* makeAsimov1DScale( const RooAbsPdf& pdf, RooRealVar& obs, const RooBinning& binning, const char* binningName, const double& scale, 
				  const Bool_t& overwrite, RooRealVar& overwriteObs );

  void setAverageErrors( RooDataHist& data, const RooDataHist& input1, const RooDataHist& input2);

  RooDataHist* makeAsimov1DError( const RooAbsPdf& pdf, RooRealVar& obs, const RooBinning& binning, const char* binningName, const RooFitResult& result );
  
  void copyScaleAndErrors( RooDataHist& data, const RooDataHist& input) ;

}
#endif
