/*
 */
#include <stdlib.h>
#include <sstream>
#include <stdexcept>

#include <RooRealVar.h>
#include "RooStats/HistFactory/FlexibleInterpVar.h"

#include "HZZWorkspace/SysText.h"
#include "HZZWorkspace/Helper.h"

//-----------------------------------------------------------------------------
// Create a systematic uncertainty object based on a text file
// * Apply systematics cutoff ("all" case)
// * Does not process shape systematics -- norm only
// * Shape sysematics handled through SampleBase and its derived classes 
//-----------------------------------------------------------------------------

SysText::SysText(const char* text_path):
    file_name_(text_path)
{
    sys_all_.clear();
    sys_global_.clear();
    Helper::readConfig(text_path, '=', sys_all_);

    // please comment it out after debugging!
    // Helper::printDic<string>(sys_all_);

    ch_name_ = "";
    SetCutoff(Helper::getSysCutoff("all"));
    SetMCStatCutoff(Helper::getSysCutoff("mcstat"));
    Clear();
    m_cutoff = 0.;
    m_cutoff_mcstat = 0.;
}

SysText::SysText(const map<string, map<string, string> >& sys_all)
{
    sys_all_.clear();
    sys_global_.clear();
    for(const auto sys_item : sys_all) {
        sys_all_[sys_item.first] = sys_item.second;
    }
    ch_name_ = "";
    SetCutoff(Helper::getSysCutoff("all"));
    SetMCStatCutoff(Helper::getSysCutoff("mcstat"));
    Clear();
}

SysText::SysText(){
    sys_all_.clear();
    ch_name_ = "";
    SetCutoff(Helper::getSysCutoff("all"));
    SetMCStatCutoff(Helper::getSysCutoff("mcstat"));
    Clear();
}

SysText::~SysText(){
    TIter np_iter(np_.createIterator());
    RooRealVar* np;
    while((np=(RooRealVar*)np_iter())){
        delete np;
    }
}

void SysText::Clear() {
    sys_channel_.clear();
    low_.clear();
    high_.clear();
    np_.removeAll();
}

bool SysText::SetChannel(const char* ch_name){
    Clear();
    ch_name_ = string(ch_name);
    bool has_channel = false;
    log_info("sys cutoff set to %.6f",m_cutoff);
    for(const auto& sec : sys_all_) {
        if(sec.first.compare(ch_name) == 0){
            has_channel=true;
            for(const auto& npName : sec.second) {
                istringstream iss(npName.second);
                float low_value, high_value;
                iss >> low_value >> high_value ;

                if (std::isnan(low_value) || std::isnan(high_value)) continue;

                float abs_low, abs_high, abs_var;
                abs_low = fabs(low_value-1.0);
                abs_high = fabs(high_value-1.0);
                abs_var = fabs(high_value-low_value);

                // mc stat cutoff
                if( (npName.first).find("MCSTAT") != string::npos ) {
                  if (!(abs_low>m_cutoff_mcstat || abs_high>m_cutoff_mcstat || 0.5*abs_var>m_cutoff_mcstat)) {
                      log_info("mc_stat below cutoff (%.5f)! dropping %s in chan %s with effect %.5f %.5f %.5f",m_cutoff_mcstat,(npName.first).c_str(),ch_name,abs_low,abs_high,abs_var);
                      continue;
                  }
                  log_info("mc_stat above cutoff (%.5f)! keeping %s in chan %s with effect %.5f %.5f %.5f",m_cutoff_mcstat,(npName.first).c_str(),ch_name,abs_low,abs_high,abs_var);
                } else {
                // systematic cutoff
                  if (!(abs_low>m_cutoff || abs_high>m_cutoff || 0.5*abs_var>m_cutoff)) {
                      log_info("sys below cutoff (%.5f)! dropping %s in chan %s with effect %.5f %.5f %.5f",m_cutoff,(npName.first).c_str(),ch_name,abs_low,abs_high,abs_var);
                      continue;
                  }
                  log_info("sys above cutoff (%.5f)! keeping %s in chan %s with effect %.5f %.5f %.5f",m_cutoff,(npName.first).c_str(),ch_name,abs_low,abs_high,abs_var);
                }

                vector<float>  sys;
                sys.push_back(low_value);
                sys.push_back(high_value);

                sys_channel_[TString(npName.first)] = sys;
            }
        }
    }
    if(!has_channel) log_err("%s does not have channel: %s", file_name_.c_str(), ch_name);
    return true;
}

bool SysText::AddSys(const TString& npName)
{
    vector<float>* norm_varies;
    try{
        norm_varies =&( this->sys_channel_.at(npName));
    } catch (const out_of_range& oor) {
        // try global systematics
        try{
            norm_varies = &( this->sys_global_.at(npName));
        } catch (const out_of_range& oor) {
            return false;
        }
    }
    addSys(npName, norm_varies->at(0), norm_varies->at(1));
    return true;
}

bool SysText::AddSys(const TString& ch_name, const TString& npName)
{
    if(! ch_name.EqualTo(ch_name_.c_str()) ) SetChannel(ch_name.Data());
    return AddSys(npName);
}

void SysText::addSys(const TString& npName, double low, double up)
{
    RooRealVar* npVar;
    string npNameIn=Form("%s", npName.Data());
    TString npVarName;
    if(!npName.Contains("MCSTAT")) {
      npVar = Helper::createNuisanceVar(npName.Data());
      np_.add(*npVar);
      low_.push_back(low);
      high_.push_back(up);
      log_info("added systematic: %s", npName.Data());
    } else { // MC stat
      // scale factor to be applied to the yield: 1 + stat_errr
      // nominal:1
      // stat_err: poisson constraint
      if(npNameIn.find("gamma_stat") == string::npos)
        npVarName="gamma_stat_"+npName;
      Double_t sigma = fabs(up-low)/2.;
      npVar = new RooRealVar(npVarName.Data(), npVarName.Data(), 1.0, 0, 1+5*sigma); // use 5*sigam by default
      gamma_.add(*npVar);
      log_info("added MC stat: %s", npName.Data());
    }
}

RooAbsReal* SysText::GetSys(const char* name) {
    string nameIn=Form("%s", name);

    if (low_.size() < 1){
        log_info("%s: Tried GetSys but I don't have any systematics... returning NULL", file_name_.c_str());
        return NULL;
    }
    if(nameIn.find("gamma_stat") == string::npos) {
      auto* fiv = new RooStats::HistFactory::FlexibleInterpVar(name, name, np_, 1., low_, high_);
      // 4: piece-wise log outside boundaries, polynomial interpolation inside
      fiv->setAllInterpCodes(4);
      return fiv;
    } else { // MC stat
      if(gamma_.getSize() < 1) {
        log_info("%s: Tried GetSys for MC stat but I don't have any ... returning NULL", file_name_.c_str());
        return NULL;
      }
      RooProduct*  McStat = new RooProduct(name, name, RooArgSet(gamma_));
      return McStat;
    }
}

RooAbsReal* SysText::GetSys(){
    return GetSys(Form("fiv_%s", ch_name_.c_str()));
}

bool SysText::ReadConfig(const char* text_file){
    file_name_ = text_file;
    sys_all_.clear();
    Helper::readConfig(text_file, '=', sys_all_);
    return true;
}

bool SysText::AddConfig(const char* text_file){
  map<string, map<string, string> > sys_temp; // systematics: all categories (different value)
  Helper::readConfig(text_file, '=', sys_temp);
  for (auto sec : sys_temp)
    for (auto var : sec.second)
      sys_all_[sec.first][var.first]=var.second;
  return true;
}

void SysText::AddGlobalSys(const char* npName, float low, float up)
{
    vector<float> sys;
    sys.push_back(low);
    sys.push_back(up);
    sys_global_[TString(npName)] = sys;
}

//Happy debugging
void SysText::Print() const{
    log_info("Printing SysText: %s in %s",file_name_.c_str(),ch_name_.c_str());
    for (int i(0);i<np_.getSize();++i)
        log_info("%s: lo=%.2f, hi=%.2f",np_[i].GetName(),low_[i],high_[i]);
}


void SysText::SetCutoff(float in){
    log_info("received call to SysText::SetCutoff... overriding %.6f to %.6f",m_cutoff,in);
    m_cutoff=in;
}

void SysText::SetMCStatCutoff(float in){
    log_info("received call to SysText::SetMCStatCutoff... overriding %.6f to %.6f",m_cutoff_mcstat,in);
    m_cutoff_mcstat=in;
}
