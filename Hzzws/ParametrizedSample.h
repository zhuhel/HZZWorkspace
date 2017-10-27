/* parametrize only SampleKeys, does not work well with SampleHist
 */
#ifndef __HZZWS_PARAMETRIZEDSAMPLE_H__
#define __HZZWS_PARAMETRIZEDSAMPLE_H__

#include <vector>
#include <string>

#include <Hzzws/SampleBase.h>
using namespace std;

#include <RooStats/HistFactory/RooBSplineBases.h>

class ParametrizedSample : public SampleBase
{
public:
    ParametrizedSample(const char* name, const char* para_name = "mH", float low = 118., float hi = 119.);
    virtual ~ParametrizedSample();

    // from sample base
    virtual bool setChannel(const RooArgSet& observable, const char* channelName, bool with_sys);
    virtual bool addShapeSys(const TString& npName);
    virtual RooAbsPdf* getPDF();
    //virtual RooAbsReal*  getCoeff();
    // virtual RooAbsPdf* get_mc_constraint();

    bool AddSample(SampleBase* signal);
    void SetParaRange(const string& range_name, double low_value, double hi_value);
    void SetParaOrder(int order){ order_ = order; }
private:
    void BuildBases();

private:
    string para_name_;
    int order_;

    RooRealVar* mH_;
    vector<SampleBase*>* signal_samples_;
    vector<double>* masses_;

    RooStats::HistFactory::RooBSplineBases* bases_;
};

#endif
