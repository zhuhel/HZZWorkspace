/* Shapes are models by RooParamKeysPdf, then parametrizaed by ParametrizedSample
 * It's only used for SM Higgs
 */
#ifndef __HZZWS_SAMPLEKEYS_H__
#define __HZZWS_SAMPLEKEYS_H__

#include <TCut.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>

#include <vector>
#include <string>

#include "Hzzws/SampleBase.h"
#include "Hzzws/SysText.h"
#include "Hzzws/Helper.h"


using namespace std;
class SampleKeys : public SampleBase {
public:
    explicit SampleKeys(const char* name,
            double mH, double mH_low, double mH_hi,
            const char* minitree_dir, //minitree
            const char* shape_sys);
    virtual ~SampleKeys();
    bool SetMHRange(double low, double hi) { 
        if(low > mass_ || hi < mass_) return false;
        log_err("range is not right [%.2f, %.2f]", low, hi);
        mh_low_ = low; mh_hi_ = hi; 
        return true;
    }
    void SetRho(double rho) { rho_ = rho; }
    virtual bool setChannel(const RooArgSet&, const char* channelName, bool with_sys);
    virtual bool addShapeSys(const TString& npName);
    virtual RooAbsPdf* getPDF();
    void SaveDataSet(RooWorkspace* ws);

private:
    int get_channel_index(const char* ch_name);
    void getShapeSys();
    void init();
private:
    bool is_himass_;
    double mh_low_;
    double mh_hi_;
    double rho_;
    string minitree_dir_;
    string tree_name_;    
    vector<string>* channel_cuts_;
    vector<RooDataSet*> ds_list_;

    int ch_cut_index;
    SysText* shape_sys_handler_;
    RooDataSet* ds;
    RooRealVar* mH_var_;
    RooRealVar* m4l_var_;
};

#endif
