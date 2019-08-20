//
//    Description:  Combining each category
//
#ifndef __HZZWS_COMBINER_H__
#define __HZZWS_COMBINER_H__
#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/Coefficient.h"
#include "HZZWorkspace/Category.h"
#include "HZZWorkspace/SystematicsManager.h"

#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <RooArgSet.h>
#include <RooSimultaneous.h>
#include <RooWorkspace.h>

#include <vector>
#include <map>
using namespace std;
class Combiner{

    public:
        explicit Combiner(const char* _name, const char* _configName);

        void combine();
        virtual ~Combiner();

    private:
        // attributes
        TString ws_name_; // name of workspace
        string file_path_;
        vector<string> all_categories_;
        RooWorkspace* workspace;
        string config_name_;

        map<string, map<string, string> > all_dic; // maps for whole configuration file
        RooArgSet nuisanceSet_; // nuisance parameters
        RooArgSet globalobsSet_; // global observables
        RooArgSet obs_; // observables for minitree
        RooArgSet obs_ws_; // observables for workspace
        strmap rename_map_;

        TChain* data_chain;
        TChain* mc_chain;
        map<string, Coefficient*> allCoefficients; // map for co-efficients
        SystematicsManager* sysMan;
	
	std::string weight_var_name;//name of weight variable in minitree, default is 'weight'
	

        // operations
        string findCategoryConfig(const string& cat_name, const string& name);
        Coefficient* getCoefficient(string& name);
        void getObservables(const string& obs_str, RooArgSet& obs_ws, RooArgSet& obs_minitree, bool& adaptive);
        SampleBase* createSample(const string& name, const string& sampleargs);
        void addCutVariables(RooArgSet& ch_obs_minitree, const string& cut);
        void addDataChan(map<string, RooDataSet*>& map, TChain* chain, RooArgSet& obs, const std::string& catname, bool weighted);

        // read the overall configurations
        void readConfig(const char* _name);
        void configWorkspace(RooWorkspace* ws);
};
#endif
