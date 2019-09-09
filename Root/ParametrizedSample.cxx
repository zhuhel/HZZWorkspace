/*
 */
using namespace std;
#include <stdlib.h>

#include <RooFitExtensions/RooBSpline.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooProduct.h>
#include <RooRealSumPdf.h>

#include "HZZWorkspace/ParametrizedSample.h"

//-----------------------------------------------------------------------------
// PDF class to create RooParamKeysPDFs 
//-----------------------------------------------------------------------------

ParametrizedSample::ParametrizedSample(const char* name,
        const char* para_name, float low, float hi):
    SampleBase(name),
    para_name_(para_name)
{
    signal_samples_ = new vector<SampleBase*>();
    masses_ = new vector<double>();
    mH_ = new RooRealVar(para_name_.c_str(),para_name_.c_str(), (low+hi)/2.0, low, hi);
    bases_ = NULL;
    order_ = 3;
}

ParametrizedSample::~ParametrizedSample(){
    for(auto sample : *signal_samples_)
        if(sample) delete sample;

    if(signal_samples_) delete signal_samples_;
    if(masses_) delete masses_;
    if(mH_) delete mH_;
}

bool ParametrizedSample::AddSample(SampleBase* signal)
{
    if(!signal) return false;
    signal_samples_->push_back(signal);
    cout << "mass: " << signal->get_mass() << endl;
    masses_->push_back(signal->get_mass());
    return true;
}

bool ParametrizedSample::setChannel(const RooArgSet& observable,
        const char* channelName, bool with_sys)
{
    SampleBase::setChannel(observable, channelName, with_sys);

    if(signal_samples_->size() < 1) return false;
    bool result = true;
    for(const auto& sample : *signal_samples_){
        if(!sample->setChannel(observable, channelName, with_sys)){
            result = false;
        }
    }
    category_name_ = string(channelName);
    base_name_ = TString(Form("%s_%s", pdf_name_.c_str(), channelName));
    return result;
}

bool ParametrizedSample::addShapeSys(const TString& npName)
{
    if(signal_samples_->size() < 1) return false;
    bool result = true;
    for(const auto& sample : *signal_samples_){
        if(! sample->addShapeSys(npName)) result = false;
    }
    return result;
}

RooAbsPdf* ParametrizedSample::getPDF()
{
    if(!bases_) BuildBases();
    RooArgList* controlPoints = new RooArgList();
    for(const auto& sample : *signal_samples_){
        controlPoints->add(*(sample->getPDF()) );
    }
    const char* bs_pdf_name = Form("bs_%s", base_name_.Data());
    auto* bs_pdf = new RooStats::HistFactory::RooBSpline(bs_pdf_name, bs_pdf_name, *controlPoints, *bases_, RooArgSet());
    const char* category_pdf_name = Form("%s_Para", base_name_.Data());
    auto* one = new RooRealVar("one", "one", 1.0);
    return new RooRealSumPdf(category_pdf_name, category_pdf_name, RooArgList(*bs_pdf), RooArgList(*one));
}

void ParametrizedSample::SetParaRange(const string& range_name, double low_value, double hi_value){
    mH_->setRange(range_name.c_str(), low_value, hi_value);
}

void ParametrizedSample::BuildBases()
{
    if(masses_->size() < 1) { throw std::runtime_error("Call AddSample first please!"); }
    bases_ = new RooStats::HistFactory::RooBSplineBases(Form("bases_%s", nickname_.c_str()),
            Form("bases_%s", nickname_.c_str()), order_, *masses_, *mH_);
}
