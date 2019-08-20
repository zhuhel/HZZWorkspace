// a set of tools used for calculating p0-values or setting limits
//
#ifndef _HZZWS_ROOSTATSHELPER_H_
#define _HZZWS_ROOSTATSHELPER_H_


#include <utility>
#include <map>
#include <memory>

#include "RooStats/ModelConfig.h"
#include "RooNLLVar.h"
#include "RooWorkspace.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
#include "TFile.h"

using namespace std;

namespace RooStatsHelper{
    void setDefaultMinimize();
    void setVarfixed(RooWorkspace* ws, const char* varName, double imass);
    void setVarFree(RooWorkspace* combined, const char* varName);
    pair<double,double> getVarVal(const RooWorkspace& w, const char* var);
    // MG: old signature - keeping a wrapper to the new one around with a warning, to avoid
    // unexpected results due to implicit RooWorkspace*  --> bool conversion in signature below 
    RooFitResult* minimize(RooNLLVar* nll, RooWorkspace* ws);
    // MG: This is the version to use
    RooFitResult* minimize(RooNLLVar* nll, bool save = true, const RooArgSet* minosSet = NULL);
    RooNLLVar* createNLL(RooAbsData* data, RooStats::ModelConfig* mc);
    // Make asimov data
    void unfoldConstraints(RooArgSet& initial, RooArgSet& final, RooArgSet& obs, RooArgSet& nuis, int& counter);
    RooDataSet* makeAsimovData(RooWorkspace* combined, 
            double muval, 
            double profileMu,  // used when fit data
            const char* muName, // name of POI
            const char* mcname, // name of ModelConfig
            const char* dataname, // name of observed Data
            bool doprofile    // profile to data?
            );
    RooDataSet* makeUnconditionalAsimov(RooWorkspace* combined, RooStats::ModelConfig* mc, const char * dataname);
    // get p0-value
    double getPvalue(RooWorkspace* combined, 
            RooStats::ModelConfig* mc, 
            RooAbsData* data, 
            const char* muName,
            bool isRatioLogLikelihood = false);
    // sqrt(2* ((s+b)ln(1+s/b) - b ))
    double getRoughSig(double s, double b);

    // fit the data and store the profiled NPs
    void fitData(RooWorkspace*, const char* mcName, const char* dataName, const char* muName, double poival, map<string,double>& result);
    void fitData(RooWorkspace*, const char* mcName, RooAbsData* data, const char* muName, double poival, map<string,double>& result);

    // generate pseudo-data and perform the fit
    void generateToy(RooWorkspace* w,
            const char* poi_name,
            double poi_value,
            int seed,
            map<string, double>& res);

    //generate pseudodata
    RooAbsData* generatePseudoData(RooWorkspace* w,
            const char* poi_name, int seed);

    void randomizeSet(RooAbsPdf* pdf, RooArgSet* globs, int seed);
    void SetRooArgSetConst(RooArgSet& argset, bool flag = true);
    void SetRooArgSetValue(RooArgSet& argset, float value);

    // Scan POI
    bool ScanPOI(RooWorkspace* ws,
            const string& data_name,
            const string& poi_name,
            int total, double low, double hi,
            TTree* tree);
    void PrintExpEvts(RooAbsPdf* simPdf,
            RooRealVar* mu, const RooArgSet* observables, RooAbsData* obs = NULL);
    double GetObsNevtsOfSignal(RooSimultaneous* simPdf,
            RooRealVar* mu, const RooArgSet* observables, bool subrange);

    // Fit some variables in the workspace
    bool fixTermsWithPattern(RooStats::ModelConfig* mc, const char* pat) ;
    void fixVariables(RooWorkspace* ws, const string& inputs, RooStats::ModelConfig* mc = NULL);
    void setOtherPOIs(RooArgSet* pois, const string& poi_name, float other_poi_value=1, bool is_other_poi_const=true);

    bool CheckNuisPdfConstraint(const RooArgSet* nuis, const RooArgSet* pdfConstraint);

    double calculateSignificance(double obs_nll_min, double obs_nll_bkg);

    bool compare_TObject_byName(TObject* a, TObject* b);
}
#endif
