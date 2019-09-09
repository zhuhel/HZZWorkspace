// =========================================================================
// Parameterized histogram shape PDF class
//    Description:  Derives from  SampleBase 
//
// ==========================================================================
#include "HZZWorkspace/SampleHistParam.h"
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <RooDataHist.h>
#include <TKey.h>
#include <TString.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include <RooFitExtensions/RooExpandedDataHist.h>
#include <RooFitExtensions/RooExpandedHistPdf.h>
#include "RooNDKeysPdf.h"
#include <RooFitExtensions/Roo1DMomentMorphFunction.h>

#include "HZZWorkspace/Helper.h"

SampleHistParam::SampleHistParam(const char* _name,
        const char* _input
        ) :
    SampleBase(_name)
{
    cout <<"creating SampleHistParam"<<endl;
    cout <<"name: " << _name << endl;
    cout <<"input: " << _input << endl;

    // expect _input has:
    // one root file for signal and one root file for background
    vector<string> inputNames;
    Helper::tokenizeString(_input, ';', inputNames);
    if (inputNames.size() != 2){
        log_err("SampleHistParam: Please provide root file for signal and background, separated by ,");
        exit(1);
    }

    sig_file_ = TFile::Open(Form("%s/%s", Helper::getInputPath().c_str(), inputNames.at(0).c_str()));
    bkg_file_ = TFile::Open(Form("%s/%s", Helper::getInputPath().c_str(), inputNames.at(1).c_str()));

    if(!sig_file_ || sig_file_->IsZombie() ||
            !bkg_file_ || bkg_file_->IsZombie() )
    {
        log_err("%s does not exist", sig_file_->GetName());
    }

    h_sig_ = nullptr;
    h_hb_ = nullptr;
    h_hH_ = nullptr;
    h_bonly_ = nullptr;
    mu_ = new RooRealVar("mu", "mu", 0, 1);
}

SampleHistParam::~SampleHistParam(){
    if(sig_file_) sig_file_->Close();
    if(bkg_file_) bkg_file_->Close();
}

bool SampleHistParam::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
    SampleBase::setChannel(_obs, _ch_name, with_sys);

    if (obs_list_.getSize() != 1){
        log_err("What you are doing? Only 1D is supported");
    }
    obsname = string(obs_list_.at(0)->GetName());

    return true;
}


RooAbsPdf* SampleHistParam::getPDF(){
    return makeNominalPDF();
}

RooHistParamPdf* SampleHistParam::makeNominalPDF()
{
    // get histograms for this channel
    if(! getHistograms() ){
        log_err("Cannot find enought inputs!");
        return nullptr;
    }
    string pdf_name(Form("%s_%s", base_name_.Data(), obsname.c_str()));
    RooRealVar* obs = (RooRealVar*) this->obs_list_.at(0);
    return new RooHistParamPdf(pdf_name.c_str(), pdf_name.c_str(),
            *obs, *mu_, *h_sig_, *h_hb_, *h_hH_, *h_bonly_);

}

bool SampleHistParam::getHistograms()
{
    // Hard code here to remove...
    TString reduced_channel_name(category_name_.c_str());
    // reduced_channel_name.ReplaceAll("ggF_", "");
    // reduced_channel_name.ReplaceAll("vv", "");
    //
    // string name_(Form("%s_%s", obsname.c_str(), category_name_.c_str()));
    string name_(Form("%s_%s", obsname.c_str(), reduced_channel_name.Data()));
    log_info("SampleHistParam searches: %s", name_.c_str());
    h_sig_   = (TH1*)sig_file_->Get(Form("%s_signal", name_.c_str()));
    h_hb_    = (TH1*)sig_file_->Get(Form("%s_HB", name_.c_str()));
    h_hH_    = (TH1*)sig_file_->Get(Form("%s_hH", name_.c_str()));
    h_bonly_ = (TH1*)bkg_file_->Get(Form("%s-Nominal-%s", obsname.c_str(), reduced_channel_name.Data()));

    return (h_sig_ && h_hb_ && h_hH_ && h_bonly_);
}
