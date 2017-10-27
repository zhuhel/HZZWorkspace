#ifndef __HZZWS_CATEGORY_H__
#define __HZZWS_CATEGORY_H__

#include <string>
#include <vector>
#include <iostream>
#include <set>

#include <RooArgSet.h>
#include <RooArgList.h>

#include "Hzzws/SampleBase.h"
#include "Hzzws/SystematicsManager.h"
using namespace std;

class Category {

 public:

  // Constructor using pre-defined labels.
  explicit Category(const string& label);
  virtual ~Category();

  const string& label() const { return m_label; }
  /////////////////////////////////////
  // Add sample to this category, tell sample to work in this category
  // Add systematics in SystematicsManager
  /////////////////////////////////////
  void addSample(SampleBase* sample, SystematicsManager* sysMan);
  void setObservables(const RooArgSet& _obs);

  // get the final PDF for this category
  RooAbsPdf* getPDF();
  const RooArgList& getGlobal() { return global_obs_list_; }
  const RooArgList& getNuisance() { return nuisance_obs_list_; } 

  // To clean PDFs based on what stat they have
  void setStatThreshold(double statThreshold){m_statThreshold = statThreshold;};

 private:

  string m_label;
  // bool is_2D ; // TODO
  RooArgSet obs ;
  set<TString> nps_set;
  vector<SampleBase*> sample_container_; // used to clean up samples

  RooArgList pdfList;
  RooArgList coefList;
  RooArgList constraintList; 
  RooArgList global_obs_list_;
  RooArgList nuisance_obs_list_;


  double m_statThreshold;
};

#endif 
