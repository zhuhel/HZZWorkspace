#ifndef __HZZWS_SIMPLEMORPH_H__
#define __HZZWS_SIMPLEMORPH_H__

#include <string>
#include <vector>

#include "HZZWorkspace/SampleHist.h"
#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/Coefficient.h"

#include "RooAbsPdf.h"
#include "RooArgList.h"


using namespace std;
class SimpleMorph : public SampleBase {

  public:

    SimpleMorph(const char* name, // used to construct PDF
        const char* configfile);
    virtual ~SimpleMorph();

    virtual bool setChannel(const RooArgSet&, const char* channelName, bool with_sys);

    virtual RooAbsPdf* getPDF();

    virtual bool addShapeSys(const TString& npName);

  private:
    RooAbsPdf* ConvertToParamKeys(RooAbsPdf* in, RooArgList& obs, const float mass, RooAbsReal *sf);

  protected:

    std::map<float, SampleBase*> m_shapemap;
    RooRealVar* m_basevar;
    TString m_morphType;
    RooStats::HistFactory::RooBSplineBases* m_bases;
    std::map<std::string, Coefficient*> m_paramScaleFactor;
    
    bool m_keepDisconnected;


};
#endif
