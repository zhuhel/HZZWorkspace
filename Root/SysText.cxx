/*
 */
#include <stdlib.h>
#include <sstream>
#include <stdexcept>

#include <RooRealVar.h>
#include "RooStats/HistFactory/FlexibleInterpVar.h"

#include "Hzzws/SysText.h"
#include "Hzzws/Helper.h"


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
    Clear();
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
    Clear();
}

SysText::SysText(){
    sys_all_.clear();
    ch_name_ = "";
    SetCutoff(Helper::getSysCutoff("all"));
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

                // systematic cutoff
                if (!(fabs(low_value-1.0)>m_cutoff || fabs(high_value-1.0)>m_cutoff || 0.5*fabs(high_value-low_value)>m_cutoff)) {
                    log_info("sys below cutoff (%.5f)! dropping %s in chan %s with effect %.5f %.5f %.5f",m_cutoff,(npName.first).c_str(),ch_name,fabs(low_value-1.0),fabs(high_value-1.0),fabs(high_value-low_value));
                    continue;
                }
                log_info("sys above cutoff (%.5f)! keeping %s in chan %s with effect %.5f %.5f %.5f",m_cutoff,(npName.first).c_str(),ch_name,fabs(low_value-1.0),fabs(high_value-1.0),fabs(high_value-low_value));


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
    RooRealVar* npVar = Helper::createNuisanceVar(npName.Data());
    np_.add(*npVar);
    low_.push_back(low);
    high_.push_back(up);
    log_info("added systematic: %s", npName.Data());
}

RooAbsReal* SysText::GetSys(const char* name) {
    if (low_.size() < 1){ 
        log_info("%s: Tried GetSys but I don't have any systematics... returning NULL", file_name_.c_str());
        return NULL;
    }
    auto* fiv = new RooStats::HistFactory::FlexibleInterpVar(name, name, np_, 1., low_, high_);
    // 4: piece-wise log outside boundaries, polynomial interpolation inside
    fiv->setAllInterpCodes(4); 
    return fiv;
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

