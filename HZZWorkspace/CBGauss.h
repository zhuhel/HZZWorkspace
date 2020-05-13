
#ifndef __HZZWS_CBGAUSS_H__
#define __HZZWS_CBGAUSS_H__

#include <string>
#include <fstream>
#include <map>
#include <vector>
using namespace std; // don't remove it otherwise RooBSpline crashes

#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooProduct.h>
#include <RooAbsReal.h>
#include <RooFitExtensions/RooBSplineBases.h>

#include "RooRealVar.h"

#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/Helper.h"


class RooWorkspace;
namespace RooStats {
  namespace HistFactory {
    class FlexibleInterpVar;
  }
}

using namespace std;

class CBGauss : public SampleBase
{

  public:

    CBGauss(const char* name, // used to construct PDF
        const char* input,        // input text file contains parameters
        const char* shape_sys,    // input text file contains shape variations
        int _doConv, //bool to do convolution (LWA) or not (NWA)
        bool _doSys,
        bool _add_int=false
        );
    virtual ~CBGauss();

    virtual RooAbsPdf* getPDF();

    virtual bool setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys);

    virtual bool addShapeSys(const TString& npName);

    // Add shape systematics
    void addSysShape(float m, const string& mean_input, const string& sigma_input);
    //void addSysShape(float mass,
    //        const map<string, map<string, string> >& sys_mean,
    //        const map<string, map<string, string> >& sys_sigma);
    void addSysShape(
            const map<float, TString>& input_mean,
            const map<float, TString>& input_sigma
            );

    void addSysSigma(float m, string file_name);
    void addSysSigma(float mass, const map<string, map<string, string> >& sys_sigma);
    void addSysSigma(const map<float, TString>& input);


    // set parameters used in CBGauss
    bool setParameters(const map<string, pair<double, double> >& input);

  private:

    void makeCBGParameterization();
    void loadTextInputFile();
    bool readTextInputFile(string ParameterName, pair<double,double>& pars);
    RooAbsReal* variable(const string& parname);

    RooStats::HistFactory::FlexibleInterpVar* flexibleInterpVar(const string& fivName, 
        vector<string>& names, vector<double>& lowValues, vector<double>& highValues);
    void BuildBases();

    RooAbsReal* getShapeSys(std::string name);

  protected:

    string inputParameterFile;

    map< string, pair<double, double> > textInputParameterValues;

  private:
    RooWorkspace* workspace;
    RooRealVar* mH;
    int doConv;
    bool doSys;
    bool addInt;

    int order_;
    vector<double>* masses_;
    vector<SysText*> shape_mean_sys_;
    vector<SysText*> shape_sigma_sys_;
    RooStats::HistFactory::RooBSplineBases* bases_;

    vector<string>* shape_sys_names_;
};
#endif
