// =========================================================================
// 
//    Description:  
// 
// ==========================================================================
#include "HZZWorkspace/SampleFactory.h"
#include "HZZWorkspace/SimpleMorph.h"
#include "HZZWorkspace/CoefficientFactory.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <sys/stat.h>

#include <RooFitExtensions/RooBSpline.h>
#include <RooFitExtensions/RooParamKeysPdf.h>
#include "RooMomentMorph.h"
#include "TVectorD.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooRealSumPdf.h"


SimpleMorph::SimpleMorph(const char* name,
        const char* configfile) : SampleBase(name),
  m_bases(NULL),
    m_keepDisconnected(false)
{
    //Read inputs from config file
    strdic cfg_dic;
    Helper::readConfig((Helper::getInputPath()+configfile).c_str(), '=', cfg_dic);
    Helper::printDic<string>(cfg_dic);

    m_morphType = cfg_dic["config"]["type"].c_str(); //"spline", or "momentmorph", or "paramspline"
    m_morphType.ToLower();

    m_keepDisconnected = (cfg_dic["config"]["keepDisconnected"]=="true");

    //configure the morphing sample and their bases
    for (auto& s: cfg_dic["bases"]){
        TString sbase = s.first.c_str();
        std::string ssample = s.second;
        float base = sbase.Atof();

        strvec samplecfg;
        Helper::tokenizeString(ssample,':',samplecfg);
        std::string sampletype = samplecfg[0];
        strvec sampleargs;
        Helper::tokenizeString(samplecfg[1],',',sampleargs);

        auto newsample = SampleFactory::CreateSample(sampletype,sampleargs);
        if (newsample) m_shapemap[base]= newsample;
    }

    //configure the morphing variable
    std::string varname = cfg_dic["config"]["morphvar"];
    float min = m_shapemap.begin()->first;
    float max = m_shapemap.rbegin()->first;



    //record any extra syst to be added
    std::string paramSyst="";
    auto paramSystDic = cfg_dic["paramSyst"];
    if (paramSystDic.size() && !m_morphType.Contains("param") && m_keepDisconnected==false)
        log_warn("received ParamSyst arguments but not using ParamSpline morph type... these will be ignored unless you enabled keepDisconnected");
    for (auto& p : paramSystDic){
            m_paramScaleFactor[p.first] = CoefficientFactory::CreateCoefficient(p.second);
            m_paramScaleFactor[p.first]->setName(Form("multiplyBy_SF%s%s",p.first.c_str(),nickname_.c_str()));
            if (p.first=="mH") m_paramScaleFactor[p.first]->SetCutoff(Helper::getSysCutoff("scale"));
            if (p.first=="Res") m_paramScaleFactor[p.first]->SetCutoff(Helper::getSysCutoff("res"));
    }

    m_basevar = new RooRealVar(varname.c_str(),varname.c_str(),min,max);
    Helper::addPoiName(varname);

}

SimpleMorph::~SimpleMorph(){}

bool SimpleMorph::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
  SampleBase::setChannel(_obs, _ch_name, with_sys);
  bool ret=true;
  for (auto& s : m_shapemap) if (!s.second->setChannel(_obs,_ch_name,with_sys)) ret=false;
  for (auto& p : m_paramScaleFactor) if (!p.second->setChannel(_ch_name,with_sys)) ret=false;
  return ret;
}

RooAbsPdf* SimpleMorph::getPDF(){

    if (m_morphType.Contains("momentmorph")){
        std::vector<double> basevals;
        for (auto& s:m_shapemap) basevals.push_back(s.first);
        TVectorD basevalsV(basevals.size(),&(basevals[0]));

        RooArgList controlPoints;
        for (auto& s:m_shapemap) controlPoints.add(*(s.second->getPDF()));

        auto p = new RooMomentMorph(Form("%s_Morph",base_name_.Data()),"", *m_basevar, obs_list_, controlPoints, basevalsV);
        p->useHorizontalMorphing(true);
        return p;
    }
    else if (m_morphType.Contains("spline")){
        std::cout<<"inside "<<__FILE__<<" "<<__func__<<std::endl;
        if (!m_bases){
            std::vector<double> basevals;
            for (auto& s:m_shapemap) basevals.push_back(s.first);
            m_bases = new RooStats::HistFactory::RooBSplineBases(Form("bases_%s", nickname_.c_str()),"", 3, basevals , *m_basevar);
        }

        std::map<std::string, RooAbsReal*> realCoefMap;
        for (auto& p: m_paramScaleFactor)
            realCoefMap[p.first] = p.second->getCoefficient(Form("%s_systFactor_%s",p.first.c_str(),category_name_.c_str()));

        RooAbsReal* multiplicativeScaleFactor = NULL;
        for (auto& a : realCoefMap){
            if (a.first==m_basevar->GetName()) multiplicativeScaleFactor = a.second;
            else if (m_keepDisconnected) Helper::getDisconnectedArgs().add(*a.second);
            else delete a.second;
        }

        RooArgList controlPoints;
        for (auto& s:m_shapemap) {
            RooAbsPdf* pdf = s.second->getPDF();
            if (m_morphType.Contains("param")) { 
                RooAbsPdf* pdftmp = ConvertToParamKeys(pdf, obs_list_, s.first, multiplicativeScaleFactor); 
                delete pdf; pdf = pdftmp; 
            }
            controlPoints.add(*pdf);
        }

        auto* bs_pdf = new RooStats::HistFactory::RooBSpline(Form("bs_%s",base_name_.Data()), "", controlPoints, *m_bases, RooArgSet());
        return new RooRealSumPdf(Form("%s_Morph",base_name_.Data()),"",RooArgList(*bs_pdf),RooArgList(RooFit::RooConst(1.0)));
    }
    else {
        std::cout<<__FILE__<<" received unrecognized morph type (received: "<<m_morphType<<"). Recognized options are:\n \tMomentMorph,\n \tSpline,\n \tParamSpline (converts pdfs to ParamKeys before doing spline) "<<std::endl;
        return NULL;
    }
}

bool SimpleMorph::addShapeSys(const TString& npName){
  bool ret=false;
  for (auto& s:m_shapemap) if (s.second->addShapeSys(npName)) ret=true;
  for (auto& s:m_paramScaleFactor) if (s.second->AddSys(npName)) ret=true;
  return ret;
}

RooAbsPdf* SimpleMorph::ConvertToParamKeys(RooAbsPdf* in, RooArgList& obs, const float mass, RooAbsReal* sf){
  log_info("converting %s into a RooParamKeysPdf", in->GetName());

  if (obs.getSize()>1) {
    log_err(" no support for ParamSpline with greater than 1D observables!");
  }

  RooRealVar& x = *((RooRealVar*)(&obs[0]));
  RooRealVar wgt("wgt","",0,1000);

  RooDataSet* ds = new RooDataSet("dstmp","",RooArgSet(x,wgt),"wgt");

  TH1* hist = in->createHistogram("htemp",x,RooFit::Binning(400));
  for (int b(0);b<=hist->GetNbinsX();++b){
    x.setVal(hist->GetBinCenter(b));
    ds->add(RooArgSet(x,wgt), hist->GetBinContent(b));
  }
  delete hist;


  auto pdf = new RooParamKeysPdf(Form("%s_paramkeys",in->GetName()),"",x, *m_basevar, mass, (sf?*sf:RooFit::RooConst(1.)), *ds, RooParamKeysPdf::MirrorBoth, 0.1);
  
  /*
  auto frame = x.frame();
  ds->plotOn(frame);
  m_basevar->setVal(mass);
  pdf->plotOn(frame);
  TCanvas* c = new TCanvas("c","c",800,800);
  frame->Draw();
  c->Print(Form("plots/Validplot_ParamKeys_%s.png",in->GetName()));
  */

  delete ds;

  return pdf;

}





























