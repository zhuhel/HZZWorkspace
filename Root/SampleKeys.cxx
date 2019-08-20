/*
 */
#include "HZZWorkspace/SampleKeys.h"
#include "HZZWorkspace/Helper.h"

#include <RooDataSet.h>
#include <RooFitExtensions/RooParamKeysPdf.h>
#include <TString.h>

#include <stdlib.h>
#include <fstream>
#include <stdexcept>

SampleKeys::SampleKeys(const char* name,
        double mH, double mH_low, double mH_hi, 
        const char* minitree_dir, 
        const char* shape_sys) :
    SampleBase(name),
    mh_low_(mH_low),
    mh_hi_(mH_hi),
    minitree_dir_(minitree_dir)
{
    mass_ = mH;
    is_himass_ = mH > 140;
    shape_sys_handler_ = new SysText(Form("%s/%s", Helper::getInputPath().c_str(), shape_sys));
    init();
}

void SampleKeys::init(){
    tree_name_ = "tree_incl_all";
    /* channel_cuts_ is analysis dependant 
     * should not be here..
     * Try to improve that later.
     * ***/
    channel_cuts_ = new vector<string>();
    channel_cuts_->push_back("event_type==0"); // 4mu
    if(is_himass_){
        channel_cuts_->push_back("event_type==2 || event_type==3"); // 2mu2e or 2e2mu
    } else {
        channel_cuts_->push_back("event_type==2"); // 2mu2e
        channel_cuts_->push_back("event_type==3"); // 2e2mu
    }
    channel_cuts_->push_back("event_type==1"); // 4e
    
    ch_cut_index = -1;
    rho_ = 1.0;
    ds = NULL;
    mH_var_ = nullptr; 
    m4l_var_ = nullptr;
}

SampleKeys::~SampleKeys(){
    if(channel_cuts_) delete channel_cuts_;
    if(shape_sys_handler_) delete shape_sys_handler_;
    if(ds_list_.size() > 0) {
        for(auto& ds : ds_list_){
            if(ds) delete ds;
        }
    }
    if(mH_var_) delete mH_var_;
    if(m4l_var_) delete m4l_var_;
}

bool SampleKeys::setChannel(const RooArgSet& obs, const char* ch_name, bool with_sys)
{
    SampleBase::setChannel(obs, ch_name, with_sys);
    if(with_sys && shape_sys_handler_) shape_sys_handler_->SetChannel(ch_name);
    ch_cut_index = get_channel_index(ch_name); 
    return true;
}

bool SampleKeys::addShapeSys(const TString& npName)
{
    if(shape_sys_handler_) return shape_sys_handler_->AddSys(npName);
    else return false;
}

int SampleKeys::get_channel_index(const char* ch_name){
    TString chName(ch_name);
    int res = -1;
    if(chName.Contains("4mu")){
        res = 0; 
    } else if (chName.Contains("2mu2e")) {
        res = 1;
    } else if(chName.Contains("2e2mu")) {
        res = is_himass_?1:2;
    } else if(chName.Contains("4e")) {
        res = is_himass_?2:3;
    } else {
        log_err("I don't know this channel %s", ch_name);
    }
    cout <<"Channel index: " << res << endl;
    return res;
}

RooAbsPdf* SampleKeys::getPDF()
{
    if(mH_var_ == nullptr) mH_var_ = new RooRealVar("mH", "mH", mass_, mh_low_, mh_hi_);
    if(m4l_var_ == nullptr) m4l_var_ = new RooRealVar("m4l_constrained", "m4l", mh_low_, mh_hi_);
    // RooRealVar m4l_var_("m4l_constrained", "m4l", mh_low_, mh_hi_);
    m4l_var_->Print();
    RooAbsReal* multiplier = nullptr; 
    if(shape_sys_handler_) {
        RooAbsReal* tmp = shape_sys_handler_->GetSys(Form("fiv_shape_%s", base_name_.Data()));
        if(tmp) multiplier = tmp;
        else multiplier = new RooRealVar("one", "one", 1.0);
    } else {
        multiplier = new RooRealVar("one", "one", 1.0);
    }
    TChain* chain =  Helper::loader(minitree_dir_, tree_name_);
    // Add variables used in Cut
    RooRealVar event_type("event_type", "event_type", -10, 10);
    RooRealVar prod_type("prod_type", "prod_type", -10, 10);
    RooRealVar weight("weight", "weight", 0, 10000);

    ds = new RooDataSet(Form("ds_%s",base_name_.Data()), "data set", 
            RooArgSet(*m4l_var_, event_type, prod_type, weight), 
            RooFit::Import(*chain), 
            RooFit::Cut(channel_cuts_->at(ch_cut_index).c_str())
            // RooFit::WeightVar("weight")
            );
    m4l_var_->SetName("m4l"); //FIXME cheap hardcoding
    ds->changeObservableName("m4l_constrained","m4l"); //FIXME cheap hardcoding
    ds_list_.push_back(ds);
    ds->Print();
    auto* keyspdf = new RooParamKeysPdf(Form("key_%.f_%s",mass_, base_name_.Data()), 
            "key pdf", *m4l_var_, *mH_var_, mass_, *multiplier, 
            *ds, RooParamKeysPdf::NoMirror,  rho_);
    keyspdf->Print();
    delete chain;
    return keyspdf;
}

void SampleKeys::SaveDataSet(RooWorkspace* ws) {
    if(ds) {
        ws->import(*ds);
    }
}
