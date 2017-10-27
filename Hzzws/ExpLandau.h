
#ifndef __HZZWS_EXPLANDAU_H__
#define __HZZWS_EXPLANDAU_H__

#include <string>
#include <fstream>
#include <map>
using namespace std;

#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooProduct.h>
#include <RooAbsReal.h>
#include <RooStats/HistFactory/RooBSplineBases.h>

#include "RooRealVar.h"
#include "RooBinning.h"

#include "Hzzws/SampleBase.h"


class RooWorkspace;
namespace RooStats {
  namespace HistFactory {
    class FlexibleInterpVar;
  }
}

class ExpLandau : public SampleBase {

  public:

    ExpLandau(const char* name, // used to construct PDF
        const char* input,        // input text file contains parameters
        const char* shape_sys,    // input text file contains shape variations
        bool _doSys
        );
    virtual ~ExpLandau();

    virtual bool setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys);
    virtual RooAbsPdf* getPDF();
    virtual bool addShapeSys(const TString& npName);


  private:

    RooProduct* variable(const string& parname);

    RooStats::HistFactory::FlexibleInterpVar* flexibleInterpVar(const string& fivName, 
        vector<string>& names, vector<double>& lowValues, vector<double>& highValues);

  private:
    RooWorkspace* workspace;
    map<string, map<string, string> > para_dic;
    bool doSys;
    vector<string>* shape_sys_names_;
};
#endif
