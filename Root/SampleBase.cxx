// =========================================================================
//
//    Description:
//
// ==========================================================================
#include "HZZWorkspace/SampleBase.h"
#include <iostream>

#include <RooRealVar.h>
#include <RooUniform.h>

#include "HZZWorkspace/Helper.h"

SampleBase::SampleBase(const char* _name):
    pdf_name_(_name)
{
    nickname_ = pdf_name_.substr(pdf_name_.find_last_of('_')+1);
    mc_constraint = nullptr;
    use_adpt_bin_ = false;
    nom_pdf = nullptr;
    coef_=NULL;
}

SampleBase::~SampleBase(){
    if(nom_pdf != nullptr)  delete nom_pdf;
    if(mc_constraint != nullptr) delete mc_constraint;
}

bool SampleBase::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
    // if (mc_constraint != nullptr) delete mc_constraint; // if do so, program crash....
    mc_constraint = nullptr;
    nom_pdf = nullptr;

    // set category name
    category_name_ = string(_ch_name);
    log_info(" SampleBase: %s set to category: %s",pdf_name_.c_str() ,category_name_.c_str());

    // set observables
    obs_list_.removeAll();
    auto* iter =  _obs.createIterator();
    TIter next(iter);
    RooRealVar* var;
    while ( (var = (RooRealVar*) next()) ){
        obs_list_.add(*var);
    }

    //set normalization channel
    if (coef_){
      coef_->setName(pdf_name_);
      if (!coef_->setChannel(_ch_name, with_sys)) return false;
    }
    else
      log_warn("Trying to set channel of coefficient of %s, but no coefficient exists! I hope you know what you're doing!", pdf_name_.c_str());

    base_name_ = TString(Form("%s_%s", pdf_name_.c_str(), category_name_.c_str()));

    return true;
}

bool SampleBase::addNormSys(const TString& npName){
  return coef_->AddSys(npName);
}


RooAbsPdf* SampleBase::get_mc_constraint(){
    return mc_constraint;
}

void SampleBase::useAdaptiveBinning(){
    use_adpt_bin_ = true;
}

void SampleBase::addCoefficient(Coefficient* c){
  if (c)
    coef_ = c;
  else 
    log_err("Empty coefficient given! :(");
}

RooAbsReal* SampleBase::getCoefficient(const char* customname){
  return coef_->getCoefficient(customname);
}

RooAbsPdf* SampleBase::getPDF()
{
    string pdfname(Form("%s_COUNT", base_name_.Data()));
    return new RooUniform(pdfname.c_str(), pdfname.c_str(), obs_list_);
}
