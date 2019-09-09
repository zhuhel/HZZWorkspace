// =========================================================================
// Counting PDF
//    Description:
//  * Just a number
//  * Shape is a uniform distribution
//  * Derives from SampleBase 
// ==========================================================================
#include "HZZWorkspace/SampleCount.h"
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <RooDataHist.h>
#include <TKey.h>
#include <TString.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <RooProdPdf.h>
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include <RooFitExtensions/RooExpandedDataHist.h>
#include <RooFitExtensions/RooExpandedHistPdf.h>
#include "RooNDKeysPdf.h"
#include <RooFitExtensions/Roo1DMomentMorphFunction.h>

#include "HZZWorkspace/Helper.h"
#include "HZZWorkspace/BinningUtil.h"

SampleCount::SampleCount(const char* _name)
    : SampleBase(_name)
{
}

SampleCount::~SampleCount(){
    if(nom_pdf != nullptr)  delete nom_pdf;
    if(mc_constraint != nullptr) delete mc_constraint;
}

void SampleCount::setMCCThreshold(float thresh){
    if (thresh > 0) {
        use_mcc_ = true;
        thresh_ = thresh;
    } else {
        use_mcc_ = false;
    }
}

RooAbsPdf* SampleCount::get_mc_constraint(){
    log_info("in %s ",__func__);

    RooArgSet gammas=coef_->GetGammas();

    TIterator* iter = gammas.createIterator();
    TIter next(iter);
    // RooRealVar* gamma;
    RooRealVar* var;
    while ( (var = (RooRealVar*) next()) ){
        if(! nuisance_obs_list_.contains(*var) ) nuisance_obs_list_.add(*var);
        // gamma = var;
        TString np(var->GetName());
        //add poisson constraints
        if(np.Contains("MCSTAT") && np.Contains(category_name_)) {
          auto* pois = Helper::createMCStatConstraint(np.Data(), var, &global_obs_list_);
          return pois;
        }
    }

    return NULL;
}
