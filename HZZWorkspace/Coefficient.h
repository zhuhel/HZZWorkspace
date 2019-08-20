#ifndef __HZZWS_COEFFICIENT_H__
#define __HZZWS_COEFFICIENT_H__

#include "HZZWorkspace/SysText.h"
#include "HZZWorkspace/Helper.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include <RooFitExtensions/RooBSplineBases.h>

class Coefficient {

  public:

    //constructor
    Coefficient(strmap& input);

    //destructor
    ~Coefficient();

    void setName(const std::string& in);
    bool setChannel(const char* channelName, bool with_sys);

    RooAbsReal* getCoefficient(std::string customname="");
    RooArgList* getCoefficientList();

    bool AddSys(const TString& npName);
    bool AddFactor(RooAbsReal*); //add your own factor to the coefficient

    const std::string& get_channel(){ return m_channel; }

    void SetCutoff (float in);

    // get the associated NPs for this coefficient
    const RooArgSet& GetNPs() { return np_;};
    const RooArgSet& GetGammas() {return gamma_;};

  private:
    RooArgList* m_arglist;
    strmap m_args;
    std::map<int, std::map<float, SysText*> > m_sysHandler;
    vector<TString> m_global_sys;

    //for normalization that depends on some poi (like mH)
    std::map<int, RooStats::HistFactory::RooBSplineBases*> m_bspline_bases;
    std::map<int, RooRealVar* > m_base_var;

    RooAbsReal* GetGenericFactor(const std::string& p, const std::string& a);
    RooAbsReal* GetGenericFactorUsingNormDic(const std::string& p, const std::string& a);
    RooAbsReal* GetGenericFactorUsingConfigDic(const std::string& p, const std::string& a);
    RooAbsReal* GetSystematicFactor(const std::string& p);
    RooAbsReal* GetMCStatFactor(const std::string& p);
    RooAbsReal* GetPOI(const std::string& p);
    RooAbsReal* GetGlobalFactor(const std::string& p);

    bool BuildCoefficient();

    std::string m_fullname;
    std::string m_nickname;
    std::string m_channel;

    RooAbsReal* m_builtCoef; //Coeffient does NOT own this pointer - user responsible for delete call
    std::pair<bool, float> m_customCutoff;

    strvec poi_names;
    
    RooArgSet np_;    // nuisance parameters associated to the coefficient
    RooArgSet gamma_; // MC stat NPs associated to the coefficient
};

#endif
