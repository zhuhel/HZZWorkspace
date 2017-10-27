#include "Hzzws/BinningUtil.h"

#include <iostream>

#include "TMath.h"
#include "TF1.h"
#include "TH2D.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph2D.h"
#include "TString.h"
#include "RooBinning.h"
#include "RooNDKeysPdf.h"
#include "RooArgList.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include "RooCurve.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooMsgService.h"

using namespace std;

void 
BinningUtil::getBinning( RooBinning& binning, RooRealVar& x, RooAbsReal& wfunc )
{
  double xbegin = x.getMin();
  double xend = x.getMax();

  binning.setRange(xbegin,xend);

  double x0 = xbegin;
  x.setVal( xbegin );

  // initial step size
  double wstep = wfunc.getVal();

  //cout << x.getVal() << " " << wfunc.getVal() << endl;

  TF1* wf = wfunc.asTF( x, RooArgList(), RooArgSet());

  //cout << "Staring wstep=" << wstep << endl;

  // Force at least one iteration to avoid
  // putting everything in one bin when the width at xbegin is wider
  // than xend-xbegin
  while (true) {
  //while ( x0+wstep<xend ) {

    wstep = BinningUtil::nextstep(wf,x0,x0+wstep);

    //cout << "x0 = " << x0 << "  step = " << wstep << " xend = " << xend << endl;
    x0 = x0+wstep;
    if (x0>xend) break;
    //cout << "x1 = " << x0 << " wstep " << wstep << endl;

    binning.addBoundary(x0);

    // for next iteration:
    wstep = wf->Eval(x0);
  } 
}


Double_t 
BinningUtil::nextstep(TF1* wf, Double_t x0, Double_t x1)
{
  //double w0 = wf->Eval(x0);
  //double w1 = wf->Eval(x1);

   // If "guessed" range is less than wmin, use that
  Double_t wmin = wf->GetMinimum(x0, x1);
  // Sometimes it ends up that (wmin - ((x0+wmin)-x0))==5e-15
  // so the code goes into an infinite loop
  // the second condition here is to avoid that
  // Seems a bit odd as the epsilon is 2e-16
  if (wmin>=(x1-x0) || fabs(wmin-(x1-x0))<1e-10) { 
     //cout << "Returning x1-x0 : " << x1 << " " << x0 << endl;
     return (x1-x0); } 

  // If distance between start and point at which width is minimised
  // is less than minimum width, use that
  Double_t xmin = wf->GetMinimumX(x0, x1);
  if (((xmin-x0)<=wmin) || fabs((xmin-x0)-wmin)<1e-10) { 
     //cout << "Returning wmin " << wmin << endl;
     return wmin; }

  Double_t xend(x1);
  // Update guess to the maximum of (xmin-wmin) [ x position of wmin - wmin]
  // and (x0+wmin) [start + wmin]
  if (xmin-wmin>x0+wmin) { xend=xmin-wmin; }
  else { xend= x0+wmin; }  
  //cout << "xend: " << xend << " xmin " << xmin << " wmin " << wmin << endl;

  return nextstep(wf,x0,xend);
}


void
BinningUtil::randomizeBinned( RooDataHist& myCopy )
{
  for (int i=0; i<myCopy.numEntries(); ++i) {
    myCopy.get(i);
    double w  = myCopy.weight();
    double we = myCopy.weightError( RooAbsData::SumW2 );

    double g = RooRandom::gaussian();
    double val = w + g*we;
    if (val<0) { val=0.; }

    myCopy.set( val );
  }
}


std::vector<RooDataHist> 
BinningUtil::dataVariations(RooArgList& parList, const RooDataHist& input, const double& scale, const double& cutoff)
{
  RooDataHist myCopy(input);

  std::vector<RooDataHist> upDownVec;
  double w(0), we(0);

  for (int i=0; i<myCopy.numEntries(); ++i) {
    myCopy.get(i);
    w  = myCopy.weight();
    we = scale * myCopy.weightError( RooAbsData::SumW2 );

    if ((we/w) < cutoff) { continue; }

    myCopy.set( w+we );
    upDownVec.push_back( myCopy );

    double min( w-we > 0. ? w-we : 0. );
    myCopy.set( min );
    upDownVec.push_back( myCopy );

    RooRealVar alpha( Form("alpha%d",i), Form("alpha%d",i), 0., -6., 6.);
    alpha.setError(1);
    parList.add ( alpha );

    // and reset for next iteration
    myCopy.set( w );
  }

  upDownVec.push_back( input );  

  return upDownVec;
}


std::vector<RooDataHist> 
BinningUtil::histVariations(RooRealVar& x, const std::vector<RooDataHist>& inputVec, const double& rho, const TString& options, const Int_t& nSigma, 
                            const Bool_t& rotation, const Bool_t& sortInput )
{
  std::vector<RooDataHist> histVec;

  RooNDKeysPdf* keys(0);
  TH1F* keyshist(0);

  for (unsigned int i=0; i<inputVec.size(); ++i) {
    keys = new RooNDKeysPdf( "keys", "keys", RooArgSet(x), inputVec[i], options.Data(), rho, nSigma, rotation, sortInput );
    keyshist = (TH1F*) keys->createHistogram( Form("smooth%d",i), x, RooFit::Binning(100) );

    histVec.push_back( RooDataHist( Form("smoothhist%d",i), Form("smoothhist%d",i), x, keyshist ) );

    delete keyshist;
    delete keys;
  } 

  return histVec;
}


void
BinningUtil::pdfVariations(RooArgList& pdfList, RooRealVar& x, const std::vector<RooDataHist>& inputVec)
{
  for (unsigned int i=0; i<inputVec.size(); ++i) {
    RooHistPdf histPdf( Form("histpdf%d",i), Form("histpdf%d",i), x, inputVec[i], 1 );
    pdfList.add( histPdf );
  } 
}


void 
BinningUtil::createPdf( RooArgList& parList, RooArgList& pdfList, RooRealVar& x, const RooDataHist& input, 
                        const TString& options, const double& rho, const Int_t& nSigma, const Bool_t& rotation, const Bool_t& sortInput,
                        const double& scale, const double& cutoff )
{
  // up & down histogram variations
  std::vector<RooDataHist> dataVec = BinningUtil::dataVariations(parList, input, scale, cutoff);
  // smoothed histogram variations
  std::vector<RooDataHist> histVec = BinningUtil::histVariations(x, dataVec, rho, options, nSigma, rotation, sortInput);
  // histpdf variations
  pdfVariations( pdfList, x, histVec );
   
}


RooCurve*
BinningUtil::plotOnWithErrorBand( RooPlot* frame, Double_t Z, RooAbsPdf& pdf, RooArgSet& obs, const int& nCurves, const int& nEvents, 
			       const char* binning, RooDataHist* dataHist ) 
{
  RooCurve* cenCurve = frame->getCurve() ;
  frame->remove(0,kFALSE) ;

  RooCurve* band(0) ;

  // *** Interval method ***
  //
  // Make N variations of parameters samples from V and visualize N% central interval where N% is defined from Z
  
  // Generate 100 random parameter points distributed according to fit result covariance matrix
  Int_t n = Int_t(100./TMath::Erfc(Z/sqrt(2.))) ;
  if (n<100) n=100 ;
  if (nCurves>0) n=nCurves;
  
  // Generate variation curves with above set of parameter values
  Double_t ymin = frame->GetMinimum() ;
  Double_t ymax = frame->GetMaximum() ;

  vector<RooCurve*> cvec ;

  for (int i=0 ; i<n ; i++) {
    //std::cout << "toyset: " << i << std::endl;

    RooDataHist* binned(0);
    RooDataSet* unbinned(0); 

    if (dataHist!=0) {
      binned = new RooDataHist(*dataHist);
      BinningUtil::randomizeBinned(*binned);
    } else {
      unbinned = pdf.generate( obs, nEvents );//, RooFit::Binning(binning) );
      binned = new RooDataHist("binnedData", "binned data", obs, "adaptive");  
      binned->add( *unbinned );
    }
    RooNDKeysPdf* tmpkeys = new RooNDKeysPdf("tmpkeys","tmpkeys", obs, *binned, "3am", 0.2, 3, false, false);
    tmpkeys->plotOn( frame ); //, RooFit::Normalization(1.0,RooAbsReal::RelativeExpected) ) ;
    cvec.push_back( frame->getCurve() ) ;

    // reset for next iteration
    frame->remove(0,kFALSE) ;
    if (binned!=0) { delete binned; }
    if (unbinned!=0) { delete unbinned; }
    delete tmpkeys;
  }

  frame->SetMinimum(ymin) ;
  frame->SetMaximum(ymax) ;
  
  // Generate upper and lower curve points from 68% interval around each point of central curve
  band = cenCurve->makeErrorBand(cvec,Z) ;
  
  // Cleanup 
  for (vector<RooCurve*>::iterator i=cvec.begin() ; i!=cvec.end() ; i++) {
    delete (*i) ;
  }

  delete cenCurve ;

  //if (!band) return frame;  
  return band;
}



RooDataHist* 
BinningUtil::makeAsimov1D( const RooAbsPdf& pdf, RooRealVar& obs, const RooBinning& binning, const char* binningName, TH1* ehist )
{
  RooDataHist* binnedAsimov = new RooDataHist("binnedAsimov", "adaptive binned data Asimov", obs, binningName );    

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  //double obsmin = obs.getMin();
  //double obsmax = obs.getMax();

  // obs.Print();
  assert (binning.numBins() > 0); //
  assert (obs.getMin() == binning.binLow(0)) ; //
  assert ((obs.getMax()-binning.binHigh( binning.numBins()-1 ))<1e-10) ; //

  for (int bin=0; bin<binning.numBins(); ++bin) {

    double xlo = binning.binLow(bin);
    double xhi = binning.binHigh(bin);
    //double width = binning.binWidth(bin); 
    //double center = binning.binCenter(bin);

    TString rangeName = Form("%s%d",binningName,bin);
    obs.setRange( rangeName.Data(), xlo, xhi );

    RooAbsReal* igx = pdf.createIntegral( obs, RooFit::NormSet(obs), RooFit::Range(rangeName.Data()) ) ;
    double integral = igx->getVal();
    delete igx;

    double fval = integral;// / width;
    //cout << "bin = " << bin << " range = " << rangeName << " xlo = " << xlo << " xhi = " << xhi << " integral " << integral << endl;
    if (ehist) {
       Double_t hint, herr;
       hint = ehist->IntegralAndError(ehist->FindBin(xlo), ehist->FindBin(xhi), herr);
       double werr=0;
       if (hint>0) {
         werr = fval*herr/hint;
       }
       //cout << bin << " xlo = " << xlo << " xhi = " << xhi << " w " << fval << " werr " << werr << endl;
       binnedAsimov->get(bin);
       binnedAsimov->set(fval, werr);
    } else {
       binnedAsimov->get(bin);
       binnedAsimov->set(fval);
    }
  }

  return binnedAsimov;
}



RooDataHist* 
BinningUtil::makeAsimov1DScale( const RooAbsPdf& pdf, RooRealVar& obs, const RooBinning& binning, const char* binningName,
				const double& scale, const Bool_t& overwrite, RooRealVar& overwriteObs )
{
  RooDataHist* binnedAsimov(0);
  if (!overwrite) {
    binnedAsimov = new RooDataHist("binnedAsimov", "adaptive binned data Asimov", obs, binningName );    
  } else {
    binnedAsimov = new RooDataHist("binnedAsimov", "adaptive binned data Asimov", overwriteObs, binningName ); 
  }

  //double obsmin = obs.getMin();
  //double obsmax = obs.getMax();

  //assert (binning.numBins() > 0); //
  //assert (obs.getMin() == binning.binLow(0)) ; //
  //assert (obs.getMax() == binning.binHigh( binning.numBins()-1 )) ; //

  for (int bin=0; bin<binning.numBins(); ++bin) {

    double xlo = binning.binLow(bin);
    double xhi = binning.binHigh(bin);
    //double width = binning.binWidth(bin); 
    //double center = binning.binCenter(bin);

    //cout << "bin = " << bin << " xlo = " << xlo << " xhi = " << xhi << endl;
    TString rangeName = Form("%s%d",binningName,bin);
    obs.setRange( rangeName.Data(), xlo, xhi );

    RooAbsReal* igx = pdf.createIntegral( obs, RooFit::NormSet(obs), RooFit::Range(rangeName.Data()) ) ;
    double integral = igx->getVal();
    delete igx;

    double fval = integral ; // / width;

    binnedAsimov->get(bin);
    binnedAsimov->set(fval * scale, 0.0 );
  }

  return binnedAsimov;
}


void 
BinningUtil::setAverageErrors( RooDataHist& data, const RooDataHist& input1, const RooDataHist& input2)
{
   //cout << "BinningUtil::setAverageErrors "<< data.GetName() << " " << input1.GetName() << " " << input2.GetName() << endl;
for (int bin=0; bin<data.numEntries(); ++bin) {

 input1.get(bin);
 double ew1 = input1.weightError(RooAbsData::SumW2);
 double w1 = input1.weight();
 //double few1 = (w1==0) ? 0 : ew1/w1;

 input2.get(bin);
 double ew2 = input2.weightError(RooAbsData::SumW2);
 double w2 = input2.weight();
 //double few2 = (w2==0) ? 0 : ew2/w2;
 
 data.get(bin);
 double w = data.weight();
 // Assume each input pdf contributes equally (an assumption) so average sq weights
 // Then rescale to the actual weight of the morphed pdf
 //double werr = data.weight() / (0.5*(w1+w2)) * 0.5* TMath::Sqrt(ew1*ew1 + ew2*ew2);
 double werr = 0;
 if (w1+w2>0) {
   werr = w/((w1+w2)) * TMath::Sqrt(ew1*ew1 + ew2*ew2);
 } else {
   werr = w;
 }

 // Gives fraction of w1, (1-f) w2
 // w = f*w1 + (1-f)*w2
 //double f = (w-w2) / (w1-w2);
 //double f=-1;
 //double werr = TMath::Sqrt( f*f * ew1*ew1 + (1-f)*(1-f) * ew2*ew2);

 data.set(w, werr);
 //cout << bin << " " << input1.weight() << " " << ew1 << " " << input2.weight() << " " << ew2 << " " << f << " " << w << " " << werr << endl;
}
}


void 
BinningUtil::copyScaleAndErrors( RooDataHist& data, const RooDataHist& input)
{
  double scale = input.sumEntries();

  for (int bin=0; bin<data.numEntries(); ++bin) {
    input.get(bin);
    double ew = input.weightError(RooAbsData::SumW2);

    data.get(bin);
    double w = data.weight();
    data.set(scale * w, ew );
  }
}


RooDataHist* 
BinningUtil::makeAsimov1DError( const RooAbsPdf& pdf, RooRealVar& obs, const RooBinning& binning, const char* binningName, const RooFitResult& result )
{
  RooDataHist* binnedAsimov = new RooDataHist("binnedAsimov", "adaptive binned data Asimov", obs, binningName );    

  //double obsmin = obs.getMin();
  //double obsmax = obs.getMax();

  assert (binning.numBins() > 0); //
  assert (obs.getMin() == binning.binLow(0)) ; //
  assert (obs.getMax() == binning.binHigh( binning.numBins()-1 )) ; //

  for (int bin=0; bin<binning.numBins(); ++bin) {

    double xlo = binning.binLow(bin);
    double xhi = binning.binHigh(bin);
    //double width = binning.binWidth(bin); 
    //double center = binning.binCenter(bin);

    //cout << "bin = " << bin << " xlo = " << xlo << " xhi = " << xhi << endl;
    TString rangeName = Form("%s_bin%d",binningName,bin);
    obs.setRange( rangeName.Data(), xlo, xhi );

    RooAbsReal* igx  = pdf.createIntegral( obs, RooFit::NormSet(obs), RooFit::Range(rangeName.Data()) ) ;
    double integral  = igx->getVal();
    double fvalError = igx->getPropagatedError(result);
    delete igx;

    double fval = integral; // / width;
    //fvalError /= width;

    binnedAsimov->get(bin);
    binnedAsimov->set(fval,fvalError);
  }

  return binnedAsimov;
}
