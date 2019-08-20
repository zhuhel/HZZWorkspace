
#ifndef __HZZWS_ANALYTICHMBKG_H__
#define __HZZWS_ANALYTICHMBKG_H__

#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/SysText.h"


class AnalyticHMBkg : public SampleBase {

  public:

    AnalyticHMBkg(const char* name, // used to construct PDF
        const char* input,        // input text file contains parameters
        const char* shape_sys);    // input text file contains shape variations
    virtual ~AnalyticHMBkg();

    virtual bool addShapeSys(const TString& npName);

    virtual RooAbsPdf* getPDF();

    virtual bool setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys);

  private:

    bool m_doSys;
    std::string m_parameterFile;
    std::string m_sysFile;
    std::map<string, double> m_par;
    std::map<string, SysText*> m_sys;


  protected:


  private:
    ClassDef(AnalyticHMBkg, 1) // Your description goes here...

};
#endif
