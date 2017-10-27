#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <TString.h>
#include <TSystem.h>
#include <RooAddPdf.h>
#include "RooAddition.h"
#include "RooFormulaVar.h"

#include "Hzzws/Helper.h"
#include "Hzzws/CBGaussSum.h"
#include "Hzzws/SampleFactory.h"
#include "Hzzws/Coefficient.h"
#include "Hzzws/CoefficientFactory.h"

CBGaussSum::CBGaussSum(
        const char* _name,
        const char* _input
        ) : SampleBase(_name)
{
    // create CBGauss
    map<string, map<string, string> > all_dic;
    Helper::readConfig((Helper::getInputPath()+_input).c_str(), '=', all_dic);
    Helper::printDic<string>(all_dic);
    for(const auto& cat_info : all_dic) {

      //PDF 
      std::cout<<"reading CBGauss pdf info for subchan "<<cat_info.first<<" under field 'para'"<<std::endl;
      const string& para = (cat_info.second).at("para");
      strvec args;
      Helper::tokenizeString(para, ',', args);
      auto cb_gauss = SampleFactory::CreateSample("CBGauss", args);

      //Fracs
      std::cout<<"reading CBGauss fraction info for subchan "<<cat_info.first<<" under field 'frac' (will be renormalized)"<<std::endl;
      const string& arg = (cat_info.second).at("frac");
      auto frac = CoefficientFactory::CreateCoefficient(arg);
      
      cb_gauss->addCoefficient(frac);
      
      //store them
      signalContainer_[cat_info.first] = cb_gauss;
    }
}

CBGaussSum::~CBGaussSum()
{
  for(auto signal : signalContainer_) {
    delete signal.second;
  }
}

bool CBGaussSum::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
  bool ret = SampleBase::setChannel(_obs,_ch_name,with_sys);

  std::cout<<"CBGaussSum setting category for sub CBGauss pdfs"<<std::endl;
  for (auto& s: signalContainer_){
    if (!ret) break;
    std::cout<<"setting channel to "<<s.first.c_str()<<std::endl;
    if (!s.second->setChannel(_obs,s.first.c_str(),with_sys))
      ret=false;
  }
  std::cout<<"setChannel returning: "<<(ret?"true":"false")<<std::endl;
  return ret;
}


bool CBGaussSum::addShapeSys(const TString& npName)
{
  bool res = true;
  for (auto& s: signalContainer_){
    if(!s.second->addShapeSys(npName)) {
      res = false;
      break;
    }
  }
  return res;
}

RooAbsPdf* CBGaussSum::getPDF()
{
  //Find sum of fractions to normalize it
    RooArgList allFracList_ = RooArgList();
    std::map<std::string, RooAbsReal*> coefmap;
    for (auto& s: signalContainer_) {
      std::string nick = Helper::extractDecay(s.first);
      auto f = s.second->getCoefficient(Form("%s_rel%s",base_name_.Data(),nick.c_str()));
      coefmap[s.first]=f;
      allFracList_.add(*f);
    }
    auto allFrac = new RooAddition(Form("%s_sumRel",base_name_.Data()),"",allFracList_);


    RooArgList pdfList_ = RooArgList();
    RooArgList fracList_ = RooArgList();
    for (auto& s: signalContainer_){
      std::string nick = Helper::extractDecay(s.first);
      pdfList_.add(*s.second->getPDF());
      auto f = new RooFormulaVar(Form("%s_frac%s",base_name_.Data(),nick.c_str()),"@0/@1",RooArgList(*coefmap[s.first],*allFrac));
      fracList_.add(*f);
    }


    string pdf_name(Form("%s_cbga_sum", base_name_.Data()));
    auto final_pdf = new RooAddPdf(pdf_name.c_str(), pdf_name.c_str(), pdfList_, fracList_);  
    return final_pdf;
}
