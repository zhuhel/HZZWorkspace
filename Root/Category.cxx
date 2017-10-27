#include "Hzzws/Category.h"

#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <RooRealVar.h>
#include <RooAddPdf.h>
//#include <RooRealSumPdf.h>
#include <RooGaussian.h>
#include <RooProdPdf.h>
#include <RooMCHistConstraint.h>
#include "RooAddition.h"
#include <RooHistPdf.h>
#include <RooDataHist.h>

#include "Hzzws/Helper.h"
//________________________________________________________________________
    Category::Category(const string& label)
    : m_label(label)
{
    sample_container_.clear();
    m_statThreshold = -1;
}

//________________________________________________________________________
Category::~Category() 
{
    for(auto sample : sample_container_){
        if(sample) delete sample;
    }
}

void Category::addSample(SampleBase* sample, SystematicsManager* sysMan){
    sample_container_.push_back(sample);
    bool with_sys = sysMan->totalNP() > 0;
    if (!sample->setChannel(this->obs, this ->m_label.c_str(), with_sys)){
      log_err("failed to set channel for sample %s, I am dropping it from category %s",sample->get_pdf_name().c_str(), m_label.c_str());
      return;
    }
    vector<TString>* nps = NULL;
    if(with_sys) {
        log_info("adding systematics to sample %s ",(sample->get_pdf_name()).c_str());
        nps = sysMan->add_sys(sample);
        if (nps != NULL && nps ->size() > 0)
        {
            for(auto& _np : *nps){
                nps_set.insert( _np );
            }
        }
        delete nps;
    }
    //Add the coefficient
    RooAbsReal* coefficient = sample->getCoefficient();
    if (!coefficient){
        log_err("Cannot find coefficient for %s in %s", sample->get_pdf_name().c_str(), m_label.c_str());
        return;
    }
    if (coefficient->getVal()==0){
        // log_info("Coefficient of %s in category %s has coefficient==0, it will be dropped from this category",sample->get_pdf_name().c_str(), m_label.c_str());
        log_info("Coefficient of %s in category %s has coefficient==0, CHECK please",sample->get_pdf_name().c_str(), m_label.c_str());
        // return;
    }
    if (coefficient->getVal()!=0){
        log_info("Coefficient of %s in category %s has coefficient is %f",sample->get_pdf_name().c_str(), m_label.c_str(), coefficient->getVal());
    }

    //Add the PDF
    RooAbsPdf* sample_pdf = sample->getPDF();
    if(!sample_pdf){
        log_err("Cannot find pdf for %s in %s", sample->get_pdf_name().c_str(), m_label.c_str());
        return;
    } 

    RooHistPdf* histPdf = dynamic_cast<RooHistPdf*> (sample_pdf);
    if(histPdf)
    {
        double ratio = float(sample->getStats())/histPdf->dataHist().numEntries();
        if(ratio < m_statThreshold)
        {
            log_info("CheckStats - stats less than cutoff threshold dropping %s",sample->get_pdf_name().c_str());
            log_info("CheckStats - cat: %s",m_label.c_str());
            
            return;
        }
    }
   log_info("adding coefficient:");
    coefficient->Print();
    coefList.add(*coefficient);

    sample_pdf->Print();
    pdfList.add(*sample_pdf);

    RooMCHistConstraint* mc_constrt =
        dynamic_cast<RooMCHistConstraint*>(sample->get_mc_constraint());
    // add MC constraint terms
    if (mc_constrt != nullptr){
        constraintList.add(*mc_constrt);
        TIter next_global(mc_constrt->getGlobalObservables().createIterator());
        TIter next_nuis(mc_constrt->getNuisanceParameters().createIterator());
        RooAbsReal* global;
        RooAbsReal* nuisance;
        while ((
                    global = (RooAbsReal*) next_global(),
                    nuisance = (RooAbsReal*) next_nuis()
               ))
        {
            if (nuisance->isConstant()) continue; // not include the ones below threshold.
            global_obs_list_.add(*global);
            nuisance_obs_list_.add(*nuisance);
        }
    }
}

void Category::setObservables(const RooArgSet& _obs)
{
    TIterator* iter = _obs.createIterator();
    TIter next(iter);
    RooRealVar* var;
    while ( (var = (RooRealVar*) next()) ){
        obs.add(*var);
    }
}


RooAbsPdf* Category::getPDF(){
    log_info("getPdf() creating a RooAddPdf for category %s",m_label.c_str());

    TString pdfname( Form( "modelunc_%s", m_label.c_str() ) );
    if (nps_set.empty()) pdfname = Form( "model_%s", m_label.c_str() );
    RooAbsPdf* uncon_pdf = new RooAddPdf(pdfname.Data(), pdfname.Data(), pdfList, coefList);
    uncon_pdf->Print();
    if (nps_set.empty()) return uncon_pdf;

    //add gaussian constraints
    for(auto& np : nps_set){
        auto* gauss = Helper::createConstraint(np.Data());
        constraintList.add(*gauss);
    }
    constraintList.add(*uncon_pdf);


    TString modelname( Form( "model_%s", m_label.c_str() ) );
    auto* prod_pdf = new RooProdPdf(modelname.Data(), modelname.Data(), constraintList);


    return prod_pdf;
}
