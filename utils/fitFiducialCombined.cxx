/* Simple example for fitting and creating a profile for the
 * fiducial cross section
 * Gaetano Barone <gaetano.barone@cern.ch>
 * Hannah Herde   <hannah.herde@cern.ch>
 */
// ------------------------------------------------------------------------------
// Processes total XS workspaces in HIGGS COMBINATION to run nll scans
//
// Operates on a combined HZZ & HGam workspace
//
// Produced 2016; not updated since! Use with caution
// Particular usage uncertain... Use with caution. Based on fitFiducial
//
// First argument: workspace file name
// Space-separated arguments follow:
// - asimov: use Asimov data (instead of obsData)
// - sys: use Systematics (instead of stat-only)
//
// EG:
// fitFiducialCombined {path}/combinedHZZHGamWorkspace.root no no no
// fitFiducialCombined {path}/workspace.root asimov sys
//
// Assumes:
// - RooWorkspace  name is "combined"
// - Model configuration name is "modelConfig"
// ------------------------------------------------------------------------------

#include "RooStats/ProfileInspector.h"
#include "TFile.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TList.h"
#include "TMath.h"
#include "TSystem.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TF1.h"
#include "TVectorD.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooMinuit.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooHistPdf.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ModelConfig.h"
#include "RooCategory.h"
#include <vector>
#include "stdio.h"
#include "time.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

TF1 *GetNLL( string name,ModelConfig *model,RooAbsData *data,RooRealVar *mainPoi, vector <double> &vals,vector <double> &nps,vector <double> &npErrLow,vector <double> &npErrHigh,vector <string> &names){

//::: Calculate the profile likelihood given the data and the model
  ProfileLikelihoodCalculator *calc = new ProfileLikelihoodCalculator( *data, *model );
  calc->SetConfidenceLevel( 0.683 );
  //  calc->SetNullParameters(
  LikelihoodInterval *interval = calc->GetInterval();
  interval->Print();
  double lower = 0.5*mainPoi->getVal(); double upper = 1.5*mainPoi->getVal();
  double central = ( ( RooRealVar* ) interval->GetBestFitParameters()->find( *mainPoi ) )->getVal();

  RooRealVar *npReR = NULL;
  TIter npIt = model->GetNuisanceParameters()->createIterator();
  vector< string > processed;
  while( ( npReR = ( RooRealVar* ) npIt.Next() ) ){
    nps.push_back( npReR->getVal() );
    npErrLow.push_back( npReR->getErrorLo() );
    npErrHigh.push_back( npReR->getErrorHi() );
    names.push_back( npReR->GetName() );
  }

  interval->FindLimits( *mainPoi, lower, upper );
  std::cout << "Central " << central << " Lower value: " << lower << " upper " << upper << std::endl;
  vals.push_back( central ); vals.push_back( lower ); vals.push_back( upper );
  //vals.push_back( interval->GetLikelihoodRatio()->getVal() );
  /*
  LikelihoodIntervalPlot *pll_frac = new LikelihoodIntervalPlot( interval );
  //pll_frac->SetNPoints( 2 );
  pll_frac->SetNPoints( 50 );
  RooArgSet tmps;
  tmps.add( *mainPoi ) ;

  TCanvas *LiklehoodC = new TCanvas( string( "LiklehoodC" ).c_str(), "Liklehood profile vs parameter of interest" );{
    LiklehoodC->cd();
    //pll_frac->SetPlotParameters( &tmps );
    //pll_frac->Draw( "tf1" );
    pll_frac->Draw();
  }

  LiklehoodC->ls();

  TF1 *nll1 = ( TF1* ) LiklehoodC->FindObject( "_PLL_xsec_tot_bin0" );
  TF1 *nllC = ( TF1* ) nll1->Clone( name.c_str() );
  delete LiklehoodC;;
  delete pll_frac;
  delete interval;
  */

  RooAbsReal *Nll=NULL;
  RooMinuit *msys=NULL;
  Nll=(model->GetPdf())->createNLL(*data,NumCPU(2));
  msys=new RooMinuit(*Nll);
  msys->migrad();
  msys->hesse();
  msys->minos(*mainPoi);
  RooAbsReal *pll_frac=Nll->createProfile(*mainPoi);
  //TF1 *nllC=Nll->asTF(*mainPoi);
  //TF1 *nllC=pll_frac->asTF(*mainPoi);
  double var=2.0;

  vector <double> x;
  vector <double> y;
  int npoints=50;
  for( int i=0 ; i<(int)npoints; i++){
    var+=(200.- 2.0)/(double)npoints;
    mainPoi->setVal(var);
    y.push_back(pll_frac->getVal()) ;
    x.push_back(mainPoi->getVal());
    cout<<"Value is "<<x.back()<<" @ "<<y.back()<<endl;
  }
  double xmin=x[0];
  double xmax=x.back()-1;

  TGraph * g = new TGraph((int)x.size()-1, &x[0], &y[0]);
  TF1 * nll1 = new TF1("f",[&](double *xp, double *p){ return p[0]*g->Eval(xp[0]); }, xmin, xmax, 1);
  TF1 *nllC=(TF1*)nll1->Clone(name.c_str());

  //delete msys;
  //delete pll_frac;
  //delete  msys;
  g->Write(("graph"+name).c_str());
  delete g;
  delete calc;
  return nllC;
}

int main( int argc, char *argv[] ) {
  string fname = argv[ 1 ];
  string options;
  int nCPU=1;
  if( argc > 1 )
    for( int i = 2; i < ( int ) argc; i++ )  options += argv[ i ];

  cout << "File will be " << fname << endl;
  int isAsimov = options.find( "asimov" ) != string::npos ? 1:0;
  int doSyst = options.find( "sys" ) != string::npos ? 1:0;
  int doTotal=true;

  FitOptions( SumW2Error( false ), Minos( true ), Hesse( false ), Extended(), Save( true ), Strategy(2));
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer( "Minuit" ); //Used to be Minuit2
  string snapshotname = isAsimov > 0 ? "postFitObs":"postFitExp";
  string nllname;
  if( doSyst > 0 ) nllname = isAsimov > 0 ? "expSys":"obsSys";
  else nllname = isAsimov > 0 ? "expStat":"obsStat";

  std::cout<<"The options are "<<options<<std::endl;
  std::cout<<"The output will be "<<nllname<<std::endl;

  //::: Fetch the files - one for the statistical workspace and one for the stat+sys workspace
  TFile *f = new TFile( fname.c_str(), "READ" );

  //::: Open the workspaces
  // w: includes all uncertainties "sys"; wO: statistical uncertainty only "stat"

  // sys WS - (1) Open the workspace
  RooWorkspace *w = ( RooWorkspace* ) f->Get( "combined" );
  // sys WS - (2) Retrieve modelConfig
  ModelConfig *model = ( ModelConfig* ) w->obj( "ModelConfig" );
  //const RooArgSet *preFit = model->GetSnapshot();
  // sys WS - (3) Load the datasets
  string DataSet = isAsimov  > 0 ? "asimovData":"combData";
  RooAbsData *data = w->data( DataSet.c_str() );
  model->Print();
  data->Print();

  // sys WS - (4) Retrieve the PDF components
  const RooArgSet *poi = model->GetParametersOfInterest();
  // const RooArgSet *nps = model->GetNuisanceParameters();
  //const RooRealVar *m4l = ( RooRealVar* ) model->GetObservables()->first();
  // RooCategory *cat = w->cat( "channelCat" );
  //cout<<"********** Categroy" <<endl;
  //cat->Print();

  if( doTotal > 0 ){
    RooRealVar *bzz = ( RooRealVar* ) w->var( "BZZ_hzz" );
    bzz->setVal( 0.1251 ); //do total XS in pb instead of fb
  }

  // sys WS - (5) Assign the primary parameter of interest (inclusive cross section)
  // Iterate over the parameters of interest ("poi") in the sys WS. Pull out the cross section. Allow it to vary.
  TIter itr_poi = poi->createIterator();
  RooRealVar *mainPoi = NULL;
  RooRealVar *this_poi = NULL;
  while ( ( this_poi = ( RooRealVar* ) itr_poi.Next() ) ){
    this_poi->Print(); //print out each parameter of interest to the Command Line
    if( string ( this_poi->GetName() ).compare( "xsec_tot_bin0" ) == 0 ) { //Identify the XS_all, the inclusive cross section
      this_poi->setConstant( false ); // Allow it to vary in the sys workspace
      mainPoi = this_poi; // Assign the XS_all POI as the primary parameter of interest in the sys WS
    }
  }

  mainPoi->setVal(48.0);
  mainPoi->setRange(0,300);
  RooRealVar *iv=(RooRealVar*)model->GetNuisanceParameters()->find("alpha_ATLAS_MUON_EFF_TrigSystUncertainty");
  if(iv!=NULL)iv->setConstant(true);

  RooRealVar *iv2=(RooRealVar*)model->GetNuisanceParameters()->find("alpha_ATLAS_MUON_EFF_TrigStatUncertainty");
  if(iv2!=NULL)iv2->setConstant(true);

  RooRealVar *bzz = ( RooRealVar* ) w->var( "BZZ_hzz" );
   bzz->setVal( 0.1251 ); //do Total in pb instead of fb

  RooRealVar *bgamma = ( RooRealVar* ) w->var( "BR_hgamma" );
  //bgamma->setVal(2.246 ); //do Total in pb instead of fb
  bgamma->setConstant(true);

  RooRealVar *lumigamma = ( RooRealVar* ) w->var( "lumi_hgamma" );
  //lumigamma->setVal(lumigamma->getVal()/1e3);
  lumigamma->setConstant(true);
  lumigamma->Print();

  w->var("A_bin0_hgamma")->setConstant(kTRUE);
  w->var("C_bin0_hgamma")->setConstant(kTRUE);
  w->var("Uncert_EnRes_EnRes")->setVal(0.0);
  w->var("Uncert_EnScale_EnScale")->setVal(0.0);
  w->var("Uncert_EnRes_EnRes")->setError(0.0);
  w->var("Uncert_EnScale_EnScale")->setError(0.0);
  w->var("higgs_mass_hgamma")->setConstant(kTRUE);

  /*
  //::: Iterate over the nuisance paramters ("nps")
  if( doSyst ){
    TIter itr_nps = nps->createIterator();
    RooRealVar *this_np = NULL;
    while ( ( this_np = ( RooRealVar* ) itr_nps.Next() ) ){
      if( (string(this_np->GetName())).find("BR_")!=string::npos)
        this_np->setConstant(true);

      else if( (string(this_np->GetName())).find("XS_ch_")!=string::npos)
        this_np->setConstant(true);

      else if( (string(this_np->GetName())).find("mu_")!=string::npos)
        this_np->setConstant(true);

      else {
        this_np->setConstant( false );
        this_np->setRange( -3, 3 );
      }
    }
  }
  */

  //::: Create RooPlot of m4l
  //::: Dump information about the sys WS
  cout << "All PDFS" << endl;
  RooArgSet pdfSet = w->allPdfs();
  pdfSet.Print();
  cout << endl;

  cout << "All Vars" << endl;
  RooArgSet varSet = w->allVars();
  varSet.Print();
  cout << endl;

  cout << "All GenericObjects" << endl;
  RooArgSet genSet = w->allVars();
  genSet.Print();
  cout << endl;

  cout << "All Functions" << endl;
  RooArgSet funcSet = w->allFunctions();
  funcSet.Print();
  cout << endl;

  vector< double > npPreFitError;
  vector< double > npPreFitVal;
  RooRealVar *npReR = NULL;
  TIter npIt = model->GetNuisanceParameters()->createIterator();
  vector< string > processed;
  while( ( npReR = ( RooRealVar* ) npIt.Next() ) ){
    double v = npReR->getVal();
    double r = npReR->getError();

    npPreFitError.push_back( npReR->getError() );
    npPreFitVal.push_back( npReR->getVal() );
    npReR->Print();
    cout << " v " << v << " error " << r << endl;

    cout << endl;
  }

  //::: DO THE FIT
  //::: Calculate the profile likelihood fit, given the data and the model
  // sys WS - data
  if( doTotal > 0 ){
    if( isAsimov > 0 ) mainPoi->setVal( 50.6 );
    else mainPoi->setVal( 50.6 );
    mainPoi->setRange( -0.5, 250 ); //totalXS in pb
  }
  else{
    if( isAsimov > 0 ) mainPoi->setVal( 50.6 );
    else mainPoi->setVal( 50.6 );
    mainPoi->setRange( -0.5, 250 ); // fiducial XS in fb
  }


  vector <double> results;
  vector <double> npVals;
  vector <double> npEl;
  vector <double> npEh;
  vector <string> npNames;

  RooAbsReal *NLL=(model->GetPdf())->createNLL(*data,NumCPU(nCPU));
  RooMinuit minuit(*NLL);
  minuit.migrad();
  minuit.minos(*mainPoi);
  //minuit.minos();
  results.push_back(mainPoi->getVal());
  results.push_back(mainPoi->getVal()-fabs(mainPoi->getErrorLo()));
  results.push_back(mainPoi->getVal()+fabs(mainPoi->getErrorHi()));


  RooArgSet allPars;
  allPars.add(*poi);
  allPars.add(*model->GetNuisanceParameters());
  w->saveSnapshot("floatingNP",allPars);


  TIter it14 = model->GetNuisanceParameters()->createIterator();
  RooRealVar *t14 = NULL;
  while ( ( t14 = ( RooRealVar* ) it14.Next() ) ){
    t14->setConstant( true );
  }

  w->saveSnapshot("fixedNP",allPars);
  w->loadSnapshot("fixedNP");
  model->LoadSnapshot();
  minuit.minos(*mainPoi);

  results.push_back(mainPoi->getVal());
  results.push_back(mainPoi->getVal()-fabs(mainPoi->getErrorLo()));
  results.push_back(mainPoi->getVal()+fabs(mainPoi->getErrorHi()));

  w->loadSnapshot("floatingNP");
  model->LoadSnapshot();

  //TIter it15 = model->GetNuisanceParameters()->createIterator();
  //RooRealVar *t15 = NULL;
  //while ( ( t15 = ( RooRealVar* ) it15.Next() ) ){
  //t15->setConstant( false );
  //}


  TFile *out=new TFile( ("HZZHGammaComb_"+nllname+ ".root").c_str(),"RECREATE");
  out->cd();

  TF1 *observed = GetNLL(nllname, model, data, mainPoi, results, npVals, npEl, npEh, npNames );

  cout << "This TF1 - " << nllname << " - it's never gonna give you up, never gonna let you down" << endl;


  cout << "Begin recording results..." << endl;
  TVectorD resultsV( ( int ) results.size(), &results[ 0 ] );
  resultsV.Write( ( "vals_" + nllname ).c_str());
  for( int i=0; i < ( int ) npVals.size(); i++)
    std::cout << npNames.at( i ) << " " << npVals.at( i )
              << " + " << npEh.at( i ) << " - " << npEl.at( i ) << std::endl;

  TVectorD npValsD( ( int ) npVals.size(), &npVals[ 0 ] );
  npValsD.Write( ( "npVals_" + nllname ).c_str() );

  TVectorD npElD( ( int ) npEl.size (), &npEl[ 0 ] );
  npElD.Write( ( "npEl_" + nllname ).c_str() );

  TVectorD npEhD( ( int ) npEh.size(), &npEh[ 0 ] );
  npEhD.Write( ( "npEh_" + nllname ).c_str() );


  //Save the results
  cout << "Writing to file: " << nllname << endl;
  observed->Write(nllname.c_str(),TFile::kOverwrite);

  doSyst=false;
  if( doSyst > 0 ) nllname = isAsimov > 0 ? "expSys":"obsSys";
  else nllname = isAsimov > 0 ? "expStat":"obsStat";


  w->loadSnapshot("fixedNP");

  TF1 *stat = GetNLL(nllname, model, data, mainPoi, results, npVals, npEl, npEh, npNames );
  cout << "This TF1 - " << nllname << " - it's never gonna give you up, never gonna let you down" << endl;

  TVectorD resultsVS( ( int ) results.size(), &results[ 0 ] );
  resultsVS.Write( ( "vals_" + nllname ).c_str());

  TVectorD npValsDS( ( int ) npVals.size(), &npVals[ 0 ] );
  npValsDS.Write( ( "npVals_" + nllname ).c_str() );

  TVectorD npElDS( ( int ) npEl.size (), &npEl[ 0 ] );
  npElDS.Write( ( "npEl_" + nllname ).c_str() );

  TVectorD npEhDS( ( int ) npEh.size(), &npEh[ 0 ] );
  npEhDS.Write( ( "npEh_" + nllname ).c_str() );

  //Save the results
  cout << "Writing to file: " << nllname << endl;
  stat->Write(nllname.c_str(),TFile::kOverwrite);

  mainPoi->Print();
  out->Close();
  f->Close();
  cout << "Have a nice day!!" << endl;
}
