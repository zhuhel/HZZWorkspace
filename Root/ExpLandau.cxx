// =========================================================================
// 
//    Description:  
// 
// ==========================================================================
#include "Hzzws/ExpLandau.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <RooDataHist.h>
#include <TKey.h>
#include <TString.h>
#include <TSystem.h>

#include "RooGlobalFunc.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include <RooStats/HistFactory/RooBSpline.h>
#include "RooExpandedDataHist.h"
#include "RooExpandedHistPdf.h"
#include "RooNDKeysPdf.h"
#include "Roo1DMomentMorphFunction.h"
#include "RooPolyVar.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "RooArgSet.h"
#include "RooFormula.h"
#include "RooGenericPdf.h"

#include "Hzzws/Helper.h"

using namespace RooStats;
using namespace HistFactory;

ExpLandau::ExpLandau(const char* _name, 
        const char* _input,  
        const char* _shape_sys,
        bool _doSys ) : SampleBase(_name)
    , workspace(new RooWorkspace("ExpLandau"))
    , doSys(_doSys)
{

    Helper::readConfig(Form("%s/%s", Helper::getInputPath().c_str(), _input), '=', para_dic);
    shape_sys_names_ = new vector<string>();
}

ExpLandau::~ExpLandau()
{
    delete workspace;
    if(shape_sys_names_) delete shape_sys_names_;
}

bool ExpLandau::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
    SampleBase::setChannel(_obs,_ch_name,with_sys);

    shape_sys_names_->clear();
    return true;
}

RooAbsPdf* ExpLandau::getPDF()
{
    RooRealVar* obs = (RooRealVar*) obs_list_.at(0);
    RooProduct* p0 = variable("p0");
    RooProduct* p1 = variable("p1");
    RooProduct* p2 = variable("p2");
    RooProduct* p3 = variable("p3");
    RooProduct* n0 = variable("n0");
    RooProduct* n1 = variable("n1");
    RooProduct* n2 = variable("n2");
    RooProduct* n3 = variable("n3");
    RooProduct* n4 = variable("n4");
    RooProduct* n5 = variable("n5");
    
    RooArgList* arg_list = new RooArgList();
    arg_list->add(*p0);
    arg_list->add(*p1);
    arg_list->add(*p2);
    arg_list->add(*p3);
    arg_list->add(*n0);
    arg_list->add(*n1);
    arg_list->add(*n2);
    arg_list->add(*n3);
    arg_list->add(*n4);
    arg_list->add(*n5);
    arg_list->add(*obs);
    string func_high = Form("TMath::Exp(%s+%s*%s+%s*%s*%s+%s*TMath::Power(%s,2.7))",
            p0->GetName(), p1->GetName(), obs->GetName(), p2->GetName(),
            obs->GetName(), obs->GetName(), p3->GetName(), obs->GetName());
    string func_low = Form("%s*TMath::Landau(%s,%s,%s)+%s*%s+%s*%s*%s+%s*%s*%s*%s",
            n0->GetName(), obs->GetName(), n1->GetName(), n2->GetName(),
            n3->GetName(), obs->GetName(), 
            n4->GetName(), obs->GetName(), obs->GetName(),
            n5->GetName(), obs->GetName(), obs->GetName(), obs->GetName());

    RooFormula formula_low("formula_low", func_low.c_str(), 
            RooArgList(*obs, *n0, *n1, *n2, *n3, *n4, *n5));
    RooFormula formula_hi("formula_hi", func_high.c_str(), 
            RooArgList(*obs, *p0, *p1, *p2, *p3));
    obs->setVal(300);

    const RooArgSet* obs_set = new RooArgSet(*obs);
    double low_300 = formula_low.eval(obs_set);
    double hi_300 = formula_hi.eval(obs_set);

    double val_ratio = low_300==0?0:hi_300/low_300;
    cout << "ratio: in " << category_name_ << " : " << val_ratio << endl;
    auto* ratio = new RooRealVar(Form("%s_ratio", base_name_.Data()), "ratio", val_ratio);

    arg_list->add(*ratio);

    string formula = Form("%s*(%s>300)+(%s<=300)*%s*%s",func_high.c_str(), 
            obs->GetName(), obs->GetName(), ratio->GetName(), func_low.c_str());
    string output_name = Form("%s_generic", base_name_.Data());
    auto* pdf_cat = new RooGenericPdf(output_name.c_str(),
            "pdf", formula.c_str(), *arg_list);
    pdf_cat->Print();
    pdf_cat->Print("v");
    workspace->import(*pdf_cat);
    
    delete arg_list; 
    delete obs_set;
    delete ratio;
    delete pdf_cat;
    return workspace->pdf(output_name.c_str());
}

bool ExpLandau::addShapeSys(const TString& npName)
{
    std::cout<<"called addShapeSys "<<npName<<std::endl;
    /* In this specific channel the NP name will be:
     * pdf_name + category_name + parameter_name
     * ATLAS_Signal_ggH_ggF_4mu_13TeV_GA_s_p0
     */
    bool resultCB = false;

    if(npName.Contains(category_name_.c_str()))         //Cb parameterization systematics
    {
        for(const auto& sys_name : *shape_sys_names_){
            if(npName.Contains(sys_name.c_str())){
                resultCB = true;
                break;
            }
        }
    }

    return resultCB;
}


RooProduct* ExpLandau::variable(const string& parname)
{
    string val_str;
    try {
        val_str = para_dic.at(category_name_).at(parname);
    } catch (const out_of_range& oor) {
        log_err("cannot find %s in channel %s", parname.c_str(), category_name_.c_str());
        return NULL;
    }
    vector<string> tokens;
    Helper::tokenizeString(val_str, ',', tokens);
    if(tokens.size() < 1) return NULL;
    double value = (double) atof(tokens.at(0).c_str());
    double error = 0;
    if(tokens.size() > 1) {
        error = (double) atof(tokens.at(0).c_str());
    }
    double ratio = value == 0?0:error/value;

    string name(Form("%s_%s", base_name_.Data(),parname.c_str()));
    
    // ignore the systematics that have less than per-mill effect
    if (doSys && fabs(ratio) > 1e-3)
    {
        shape_sys_names_->push_back(parname);
        RooRealVar var( Form("%s_nom", name.c_str()), Form("nominal %s", parname.c_str()), value);

        vector<string> names;
        vector<double> lowValues;
        vector<double> highValues;

        names.push_back(name);
        if(ratio > 1) {
            lowValues.push_back(1E-6);
        } else {
            lowValues.push_back(1. - fabs(ratio));
        }
        if (ratio > 2){
            log_warn("The error of %s is over 100%%!: %.4f", name.c_str(), ratio);
        }
        highValues.push_back(1. + fabs(ratio));

        FlexibleInterpVar* fiv = flexibleInterpVar(name, names, lowValues, highValues);

        RooProduct prod(name.c_str(), name.c_str(), RooArgList(var,*fiv));

        workspace->import(prod);
    }
    else {
        RooRealVar var( name.c_str(), name.c_str(), value);
        workspace->import(var);
    }

    return (RooProduct*)workspace->obj(name.c_str());
}

FlexibleInterpVar* ExpLandau::flexibleInterpVar(const string& fivName, vector<string>& names, 
        vector<double>& lowValues, vector<double>& highValues)
{
    RooArgList variables;

    for (int inp=0; inp<(int)names.size(); inp++) {
        RooRealVar* np = Helper::createNuisanceVar(names[inp].c_str());
        variables.add(*np);
    }

    FlexibleInterpVar* fiv = new FlexibleInterpVar(("fiv_"+fivName).c_str(), ("fiv_"+fivName).c_str(),
            variables, 1., lowValues, highValues);
    fiv->setAllInterpCodes(4);

    return fiv;
}
