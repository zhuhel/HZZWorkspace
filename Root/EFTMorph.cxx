// =========================================================================
// 
//    Description:  
// 
// ==========================================================================
#include "HZZWorkspace/CoefficientFactory.h"
#include "HZZWorkspace/EFTMorph.h"

#include "TMath.h"

#include <algorithm> // std::transform
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <sys/stat.h>

#include "RooAddition.h"

#include "RooStringVar.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooRealSumPdf.h"


EFTMorph::EFTMorph(const char* name,
        const char* configfile,bool _shape_BSM_only) : SampleBase(name),
    m_eftfunc(NULL),
    m_morphcoefs(NULL),
    m_morphfuncs(NULL),
    onlyShapeBSMsensitive(_shape_BSM_only)
{

  //Read inputs from config file
  Helper::readConfig((Helper::getInputPath()+configfile).c_str(), '=', morph_dic);
  Helper::printDic<string>(morph_dic);

  //Create internal coefficients here
  if (morph_dic.find("coefficients")!=morph_dic.end()){
    for (auto& c : morph_dic["coefficients"]){
      log_info("Creating custom internal coefficient for sample %s",(c.first).c_str());
      auto coef = CoefficientFactory::CreateCoefficient(c.second);
      m_sampleCoefMap[c.first]=coef;
    }
  }

}

EFTMorph::~EFTMorph(){
  if (m_eftfunc) delete m_eftfunc;
  for (auto& c: m_sampleCoefMap) delete c.second;
  if (m_morphcoefs) delete m_morphcoefs;
  if (m_morphfuncs) delete m_morphfuncs;
}


bool EFTMorph::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
  SampleBase::setChannel(_obs, _ch_name, with_sys);

  if (m_eftfunc) delete m_eftfunc;

  //set channel for sample coefficients
  for (auto& c : m_sampleCoefMap){
    c.second->setName((c.first).c_str());
    c.second->setChannel(_ch_name, with_sys);
  }
  RooArgSet* decMorphPara = new RooArgSet("DecayCouplings");
  RooArgSet* prodMorphPara = new RooArgSet("ProductionCouplings");
  RooArgList* samples = new RooArgList("Samples");

  //check morphing inputs are present
  if (morph_dic.find(_ch_name)==morph_dic.end()) { log_err("morphing config file does not contain channel %s",_ch_name); return false; }
  for (const auto s : {"file","folder","productioncouplings","decaycouplings","samples"}) {
    // basisname is optional, treated later below
    if (morph_dic[_ch_name].find(s)==morph_dic[_ch_name].end()) { log_err("morphing config file for channel %s is missing field %s",_ch_name, s); return false; }
    log_info("for channel %s found %s as %s",_ch_name,s,morph_dic[_ch_name][s].c_str());
  }

  //Load information from config as inputs for morphing
  if (morph_dic[_ch_name].find("basisname") == morph_dic[_ch_name].end()) {
      // fallback case: if basisname is not provided, assume it is
      // HiggsCharacterisation for backward compatibility
      log_info("for channel %s: 'basisname' property is not specified, assume it is HiggsCharacterisation", _ch_name);
      decMorphPara->add(createHCMorphParaSet(morph_dic[_ch_name]["decaycouplings"]));
      prodMorphPara->add(createHCMorphParaSet(morph_dic[_ch_name]["productioncouplings"]));
  }
  else { // switch according to chosen basis
      auto basisName = morph_dic[_ch_name]["basisname"];
      // Go for case insensitive property.
      std::transform(basisName.begin(), basisName.end(), basisName.begin(), ::tolower);
      if (basisName == "hc" or basisName == "higgscharacterisation") {
          log_info("for channel %s recognized 'basisname' to be HiggsCharacterisation", _ch_name);
          decMorphPara->add(createHCMorphParaSet(morph_dic[_ch_name]["decaycouplings"]));
          prodMorphPara->add(createHCMorphParaSet(morph_dic[_ch_name]["productioncouplings"]));
      }
      else if (basisName == "smeft" or basisName == "warsaw" or basisName == "higgs") {
          log_info("for channel %s recognized 'basisname' to be SMEFT or Warsaw or Higgs", _ch_name);
          decMorphPara->add(createGeneralMorphParaSet(morph_dic[_ch_name]["decaycouplings"]));
          prodMorphPara->add(createGeneralMorphParaSet(morph_dic[_ch_name]["productioncouplings"]));
      }
      else {
          log_err("Basis name '%s' not recognized for morphing.", morph_dic[_ch_name]["basisname"].c_str());
          return false;
      }
  } // end switch basis name

  //Add samples
  std::vector<std::string> samplesvec;
  Helper::tokenizeString(morph_dic[_ch_name]["samples"],',',samplesvec);
  for (auto& s: samplesvec) samples->add(*(new RooStringVar(s.c_str(),s.c_str(),s.c_str())));

  m_eftfunc = 
    new RooLagrangianMorphFunc(Form("%s_mf",base_name_.Data()),
        Form("%s_mf",base_name_.Data()),
        (Helper::getInputPath()+morph_dic[_ch_name]["file"]).c_str(), 
        (morph_dic[_ch_name]["folder"] + _obs.first()->GetName()).c_str() ,
        *prodMorphPara,
        *decMorphPara,
        *samples);
  // fix obsname automatically created by RooLagrangianMorphFunc to name provided in top level config file
  m_eftfunc->getVariables()->find(m_eftfunc->getObservable()->GetName())->SetName(obs_list_.first()->GetName());

//   information for cross section morphing
  m_eftfunc->writeMatrixToFile(m_eftfunc->getInvertedMatrix(),Form("coeffsInvMatrix_%s.txt",m_eftfunc->GetName()));
  m_eftfunc->writeFormulas(Form("formulas_%s.txt",m_eftfunc->GetName()));
  m_eftfunc->writePhysics(Form("phys_%s.txt",m_eftfunc->GetName()));
  return true;
}

RooAbsPdf* EFTMorph::getPDF(){
  
  //manually build RooRealSumPdf with my extra factors in it
  RooAbsPdf* finalpdf = new RooRealSumPdf(base_name_,Form("RooRealSumPdf of morphing terms for %s in %s",nickname_.c_str(),category_name_.c_str()), *m_morphfuncs, *m_morphcoefs);
  return finalpdf;
}


RooArgSet EFTMorph::createHCMorphParaSet(std::string parlist){

  
//   RooRealVar* lambda = dynamic_cast<RooRealVar*>(couplingsDatabase.find("Lambda"));
//   
//   if(!lambda){
//       lambda= new RooRealVar("Lambda","Lambda",1000.);
//       lambda->setConstant();
//       couplingsDatabase.add(*lambda);
//       Helper::addPoiName(lambda->GetName());
//   }

  RooRealVar* cosa = dynamic_cast<RooRealVar*>(couplingsDatabase.find("cosa"));
  if(!cosa){
      cosa= new RooRealVar("cosa","cosa",1./(TMath::Sqrt(2.)),0,1);
      couplingsDatabase.add(*cosa);
      Helper::addPoiName(cosa->GetName());
  } 

  RooArgSet* morphPara = new RooArgSet();

  std::vector<std::string> varlist;
  Helper::tokenizeString(parlist,',',varlist);

  if(varlist.size() == 0){
    log_err("WARNING createHCMorphParaSet() did not find any couplings in input list %s",parlist.c_str());
    log_err("Example: {kAzz} {kAzz,kHzz} {kHzz}");
  }

  const TString bsmpara = "kHzz,kAzz,kHdz,kHww,kAww,kHdwR,kHdwI,kHda";

  for(auto var : varlist){

    const TString sk = TString(var);
    const TString cname = (TString(var)).ReplaceAll("k","_g");// name of final coupling

    const bool isCPodd = sk.Contains("kA");
//     const bool isBSM = bsmpara.Contains(sk);

    TString formulak = isCPodd ? (TString("sqrt(1-(cosa*cosa))*")+sk) : (TString("cosa*")+sk);
//     formulak = isBSM ? (formulak + TString("/Lambda")) : formulak;

    RooRealVar* k = dynamic_cast<RooRealVar*>(couplingsDatabase.find(sk.Data()));
    if(!k){
      k= new RooRealVar(sk.Data(),sk.Data(),1.,-100,100);//FIXME I default the value to 1.0!
      couplingsDatabase.add(*k);
    } 
    Helper::addPoiName(k->GetName());

//     RooArgList list = isBSM  ? RooArgList(*cosa,*k,*lambda) : RooArgList(*cosa,*k);
    RooArgList list = RooArgList(*cosa,*k);

    RooFormulaVar* _g = dynamic_cast<RooFormulaVar*>(formulaDatabase.find(cname.Data()));
      if(!_g){
        _g= new RooFormulaVar(cname.Data(),formulak.Data(),list);
        formulaDatabase.add(*_g);
      } 
    log_info("formula %s %s",_g->GetName(),_g->GetTitle());
    _g->Print("t");
    morphPara->add(*_g);
  }
  return (*morphPara);
}


RooArgSet EFTMorph::createGeneralMorphParaSet(std::string parlist)
{
    std::vector<std::string> varlist;
    Helper::tokenizeString(parlist, ',', varlist);
    if(varlist.size() == 0){
        log_err("WARNING createGeneralMorphParaSet() did not find any couplings in input list %s", parlist.c_str());
        log_err("Example: cHW,cHB,cHWB");
    }

    RooArgSet* morphPara = new RooArgSet();
    for (auto var : varlist) {
        const TString sk = TString(var);
        RooRealVar* k = dynamic_cast<RooRealVar*>(couplingsDatabase.find(sk.Data()));
        if (!k) {
            k = new RooRealVar(sk.Data(), sk.Data(), 1., -100, 100); // FIXME I default the value to 1.0!
            couplingsDatabase.add(*k);
        }
        Helper::addPoiName(k->GetName());
        morphPara->add(*k);
    }
    return *morphPara;
}


RooAbsReal* EFTMorph::getOverallNormalization(){

  RooRealSumPdf* temppdf = m_eftfunc->getPdf();

  if (m_morphcoefs) delete m_morphcoefs;
  if (m_morphfuncs) delete m_morphfuncs;
  RooArgList* tempmorphcoefs = new RooArgList(temppdf->coefList());
  m_morphcoefs = new RooArgList();
  m_morphfuncs = new RooArgList(temppdf->funcList());

  if (tempmorphcoefs->getSize()!=m_morphfuncs->getSize()){
    log_err("morphcoefs different size than morphfuncs!"); return NULL;
  }

  RooArgList* allexp = new RooArgList(); //list to store sum events from each term

  for (int i(0);i<tempmorphcoefs->getSize();++i){
    RooArgList* prodlist = new RooArgList(); //list to store factors to make events for each term


    //get raw coefficient from morphing
    RooAbsReal* coefi = (RooAbsReal*)&((*tempmorphcoefs)[i]); 
    
    //add custom coefficient (if one exists) to add systematics to that coefficient
    RooAbsReal* sysTerm=NULL;
    TString name = (*m_morphfuncs)[i].GetName();
    std::string match="";
    for (auto & s: m_eftfunc->getSamples())
      if (name.Contains(s.c_str())) {match=s;break;}
    if (match!="" && m_sampleCoefMap.find(match)!=m_sampleCoefMap.end()){
      sysTerm =  m_sampleCoefMap[match]->getCoefficient(Form("internalcoef_%s_%s",base_name_.Data(),match.c_str()));
      sysTerm->Print();
      coefi = new RooProduct(Form("wgt_%s_%s",base_name_.Data(),match.c_str()),"",RooArgList(*coefi,*sysTerm));
    }

    m_morphcoefs->add(*coefi);

    //add that coefficient to the prod list
    prodlist->add(*coefi);

    //get morphing functions and add them to the prod list
    RooProduct* func = (RooProduct*)(&((*m_morphfuncs)[i]));
    RooArgList funccomps(func->components());
    for (int j(0);j<funccomps.getSize();++j){
      if (funccomps[j].dependsOn(obs_list_)){ //factors which depends on the observables need to be integrated out (the RooHistFuncs)
        RooAbsReal* in = ((RooAbsReal*)(&(funccomps[j])))->createIntegral(obs_list_);
        prodlist->add(*in);
      }
      else                                  //factors which do not depends on the observables can be taken as-is (the w_ formulae)
        prodlist->add(funccomps[j]);
    }
    log_info("component %i normalization built from factors:",i);
    prodlist->Print("v");

    RooProduct* p = new RooProduct(Form("coefprod_%s",(*m_morphfuncs)[i].GetName()),"",*prodlist);
    allexp->add(*p);
  }
  log_info("all morphing components have normalization from terms:");
  allexp->Print("v");

  //add up all these integrals: this is the expected # of events
  RooAddition* morphcoefsum = new RooAddition(Form("morph_coef_sum_%s",base_name_.Data()),"",*allexp);
  log_info("total expected from morphing model: %.2f",morphcoefsum->getVal());
  
  return morphcoefsum;
}

RooAbsReal* EFTMorph::getCoefficient(const char* customname){

  RooAbsReal* expectedEv = getOverallNormalization();
  if (!expectedEv) log_err("could not construct expected events term!");
  if (!onlyShapeBSMsensitive)coef_->AddFactor(expectedEv);

  return SampleBase::getCoefficient(customname);
}
  

bool EFTMorph::addShapeSys(const TString& npName){
  bool ret=false;
  //add systematics to sample coefficients
  for (auto& c: m_sampleCoefMap){
    if (c.second->AddSys(npName)) ret=true;
  }
  return ret;
}























