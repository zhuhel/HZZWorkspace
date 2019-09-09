#include "HZZWorkspace/AnalyticHMBkg.h"
#include "HZZWorkspace/Helper.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooProduct.h"

#include <sstream>

using namespace RooFit;

//-----------------------------------------------------------------------------
// Class to produce qqZZ/ggZZ background PDF for high mass analysis
//-----------------------------------------------------------------------------


ClassImp(AnalyticHMBkg)

AnalyticHMBkg::AnalyticHMBkg(const char* _name,
    const char* _inputParams,
    const char* _inputSystematics):
  SampleBase(_name)
{
  m_parameterFile = Helper::getInputPath()+"/"+_inputParams;
  m_sysFile = Helper::getInputPath()+"/"+_inputSystematics;
}

AnalyticHMBkg::~AnalyticHMBkg(){}

//This will just load parameters into the map
bool AnalyticHMBkg::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys){

  SampleBase::setChannel(_obs,_ch_name,with_sys);
  m_doSys=with_sys;

  std::map<std::string, std::map<string, string> > inputPars;
  Helper::readConfig(m_parameterFile.c_str() , '=', inputPars);

  if (inputPars.find(_ch_name)==inputPars.end()){
    log_err("failed to find parameters for Analytic Background in file %s in channel %s",m_parameterFile.c_str(),_ch_name);
    return false;
  }

  //Read in the parameter values
  for (auto& m : inputPars[_ch_name])
    m_par[m.first] = atof(m.second.c_str());

  log_info("for channel %s found parameters:",_ch_name);
  for (auto& m : m_par)
    log_info("%s = %f",m.first.c_str(), m.second);

  //Set up the sys text
  if (m_doSys){
    for (auto& m : m_par){
      m_sys[m.first] = new SysText(m_sysFile.c_str());
      //look for systematics for each parameter under [<cat> <par>]
      m_sys[m.first]->SetChannel(Form("%s %s",_ch_name,m.first.c_str()));
    }
  }

  return true;
}

bool AnalyticHMBkg::addShapeSys(const TString& npName){
  bool result=false;
  for (auto& s: m_sys)
    if (s.second->AddSys(npName)) result=true;
  return result;
}

//create the pdf
RooAbsPdf* AnalyticHMBkg::getPDF(){

    if (obs_list_.getSize()!=1){
        log_err("AnalyticalHMBkg only supports 1D pdfs! has %d observables",obs_list_.getSize());
        return NULL;
    }
    RooAbsReal* x = (RooAbsReal*)(&obs_list_[0]);

    std::map<std::string, RooAbsReal*> m_var;
    for (auto& p:m_par){
        m_var[p.first] = new RooConstVar(
                Form("%s_%s_%s",pdf_name_.c_str(),category_name_.c_str(),p.first.c_str()),
                Form("%s_%s_%s",pdf_name_.c_str(),category_name_.c_str(),p.first.c_str()),
                p.second);
    }

    //Add systematics here by replacing objects in the map
    if (m_doSys){
        for (auto & v: m_var){
            auto sys = m_sys[v.first]->GetSys(
                    Form("%s_%s_fiv_%s",pdf_name_.c_str(),category_name_.c_str(),v.first.c_str())
                    );
            if (sys){
                v.second = new RooProduct(
                        Form("%s_%s_%s_withsys",pdf_name_.c_str(),category_name_.c_str(),v.first.c_str()),"@0*@1",RooArgList(*v.second,*sys)
                        );
            }
        }
    }

    std::cout<<"Final AnalyticHM_Bkg params including systematics"<<std::endl;
    for (auto & m: m_var)
        m.second->Print();

    //Defining the pdf
    //
    auto x0 = m_var["x0"];
    std::cout<<"var x0:"; x0->Print();

    //f1(x)
    const char* f1form = "TMath::Exp(@1 + @2*@0 + @3*TMath::Power(@0,2))";
    auto f1 =new RooFormulaVar(Form("%s_%s_f1",pdf_name_.c_str(),category_name_.c_str()), f1form,     RooArgList(*x ,*m_var["a1"],*m_var["a2"],*m_var["a3"]));
    auto f1x0 =new RooFormulaVar(Form("%s_%s_f1x0",pdf_name_.c_str(),category_name_.c_str()), f1form, RooArgList(*x0,*m_var["a1"],*m_var["a2"],*m_var["a3"]));

    //f2(x)
    const char * f2form = "(0.5 + 0.5*TMath::Erf((@0-@1)/@2)) * (1./(1.0+TMath::Exp((@0-@1)/@3)))";
    auto f2 = new RooFormulaVar(Form("%s_%s_f2",pdf_name_.c_str(),category_name_.c_str()),f2form,     RooArgList(*x ,*m_var["b1"],*m_var["b2"],*m_var["b3"]));
    auto f2x0 = new RooFormulaVar(Form("%s_%s_f2x0",pdf_name_.c_str(),category_name_.c_str()),f2form, RooArgList(*x0,*m_var["b1"],*m_var["b2"],*m_var["b3"]));

    //f3(x) = exp( polynominal(c1, c2, c3, c4), 'c4' should be 'c3.7' for ggZZ
    // find how many parameters of c* in the array.
    vector<string> poly_order_names;
    for (map<string, RooAbsReal*>::iterator it = m_var.begin();
            it != m_var.end(); ++it){
        string para_name(it->first);
        if (para_name[0] == 'c'){
            poly_order_names.push_back(para_name);
            log_info("Found %s", para_name.c_str());
        }
    }
    log_info("Found %ld parameters for poly",  poly_order_names.size());
    stringstream oss;
    float eps = 1E-6;
    RooArgList f3_x0_list(*x0);
    RooArgList f3_x_list(*x);
    for(size_t i = 0; i < poly_order_names.size(); i += 1)
    {
        float order = atof(poly_order_names.at(i).substr(1).c_str()) - 1.;
        if( fabs(order) < eps && i!=0){
            log_err("please put zero order in the beginning!");
            log_err("Correct it.");
            break;
        }

        if( fabs(order) < eps) {
            oss << "@" << i+1;
        } else if (fabs(order - 1) < eps){
            oss << " + @" << i+1 << " *@0 ";
        } else if (fabs(order - 2) < eps){
            oss << " + @" << i+1 << " *@0*@0 ";
        } else{
            oss << " + @" << i+1 << " *TMath::Power(@0, "<< order << ")";
        }
        f3_x0_list.add(*m_var[poly_order_names.at(i)]);
        f3_x_list.add(*m_var[poly_order_names.at(i)]);
    }

    // const char* f3form = "TMath::Exp(@1+@2*@0+@3*TMath::Power(@0,2.0)+@4*TMath::Power(@0,2.7))";
    const char* f3form = Form("TMath::Exp(%s)", oss.str().c_str());
    auto f3 = new RooFormulaVar(Form("%s_%s_f3",pdf_name_.c_str(),category_name_.c_str()),f3form,     f3_x_list);
    auto f3x0 = new RooFormulaVar(Form("%s_%s_f3x0",pdf_name_.c_str(),category_name_.c_str()), f3form , f3_x0_list);

    //Cnorm
    auto cnorm = new RooFormulaVar(Form("%s_%s_cnorm",pdf_name_.c_str(),category_name_.c_str()),"@2/(@0+@1)",RooArgList(*f1x0,*f2x0,*f3x0));
    //step(x<x0)
    auto lowstep = new RooFormulaVar(Form("%s_%s_ltx0",pdf_name_.c_str(),category_name_.c_str()),"1.0*(@0<@1)",RooArgList(*x,*x0));
    //step(x>x0)
    auto highstep = new RooFormulaVar(Form("%s_%s_gtx0",pdf_name_.c_str(),category_name_.c_str()),"1.0*(@0>=@1)",RooArgList(*x,*x0));

    //final formula
    auto fx = new RooFormulaVar(Form("%s_%s_form",pdf_name_.c_str(),category_name_.c_str()),"(@0+@1)*@2*@3 + @4*@5",RooArgList(*f1,*f2,*lowstep,*cnorm,*f3,*highstep));


    return new RooGenericPdf(Form("%s_%s_ana_shape",pdf_name_.c_str(),category_name_.c_str()),"@0",RooArgList(*fx));
}
