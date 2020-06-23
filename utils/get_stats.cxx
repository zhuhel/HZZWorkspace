// ////////////////////////////////////////////////////////////
// This is the main interface for statistical analysis of workspaces, after you have built your workspace
// It provides flexibilities of setting cerntain parameters to constant values.
// It's very useful when one intends to evaluate the results without systematics
// For documentation see https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HiggsZZRunIIWorkspaces#Doing_statistics_with_your_works
//
// Intended for high mass limit setting 
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TSystem.h>
#include <RooWorkspace.h>
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooSimultaneous.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooMinimizer.h"

#include "HZZWorkspace/RooStatsHelper.h"
#include "HZZWorkspace/Helper.h"
#include "HZZWorkspace/runAsymptoticsCLsCorrectBands.h"



using namespace std;

int main(int argc, char** argv)
{
    if ((argc > 1 && string(argv[1]) == "help") ||
            argc < 5)
    {
        cout << argv[0] << " combined.root ws_name mu_name data_name mc_name strategy var=value,var=value option obs,exp other_poi" << endl;
        cout <<"in 'var=value', if noSys used, all systematics will be set to constant" << endl;
        cout << "option: pvalue,limit" << endl;
        cout << "other_poi: 0:1 meaning other poi set to 0 and fixed" << endl;
        return 0;
    }

    RooStatsHelper::setDefaultMinimize();
    const string input_name(argv[1]);
    const string wsName(argv[2]);
    const string muName(argv[3]);
    const string dataName(argv[4]);
    const string mcName(argv[5]);
    int opt_id = 6;
    int strategy =  1;
    if(argc > opt_id) {strategy = atoi(argv[opt_id]); } opt_id ++;
    string fix_variables = "";
    if(argc > opt_id) { fix_variables = argv[opt_id]; }  opt_id ++;
    string options = ""; // limit or pvalue
    if(argc > opt_id) { options = argv[opt_id]; }  opt_id ++;
    string data_opt = ""; // observed or expected
    if(argc > opt_id) { data_opt = argv[opt_id]; }  opt_id ++;
    float other_poi_value = 1;
    bool is_other_poi_const = true;
    string fix_other_str("");
    if(argc > opt_id) {
        vector<string> tokens;
        fix_other_str = argv[opt_id];
        Helper::tokenizeString(fix_other_str, ':', tokens);
        if(tokens.size() > 0){
            other_poi_value = atof(tokens.at(0).c_str());
        }
        if(tokens.size() > 1){
            is_other_poi_const = (bool) atoi(tokens.at(1).c_str());
        }
    }
    opt_id ++;
    if (fix_variables == "none") { fix_variables = "";}

    // summary of options
    cout<<" Input: " << input_name << endl;
    cout<<" wsName: " << wsName << endl;
    cout<<" muName: [" << muName <<"]"<< endl;
    cout<<" dataName: " << dataName << endl;
    cout<<" mcName: " << mcName << endl;
    cout<<" strategy: " << strategy << endl;
    cout<<" Fix variables: " << fix_variables << endl;
    cout<<" options: " << options << endl;
    cout<<" data option: " << data_opt << endl;
    cout<<" other_poi: " << fix_other_str << endl;
    cout<<" other_poi: " << other_poi_value << ", isFixed: " << is_other_poi_const << endl;
    // Load extra package for LWA workspace.
    // gSystem->Load("/afs/cern.ch/user/x/xju/public/RelativisticBWInt_cxx.so");
    // gSystem->Load("/afs/cern.ch/user/x/xju/public/h4l_classes/RooHistParamPdf_cxx.so");

    // read the workspace.
    auto input_file = TFile::Open(input_name.c_str(), "read");
    auto ws = (RooWorkspace*) input_file->Get(wsName.c_str());
    auto mc = (RooStats::ModelConfig*) ws->obj(mcName.c_str());
    auto obs_data = (RooDataSet*) ws->data(dataName.c_str());
    auto poi = (RooRealVar*) ws->var(muName.c_str());
    // const RooArgSet* observables = mc->GetObservables();
    // auto nuisances = mc->GetNuisanceParameters();
    // string asimov_data_name = "asimovNull";
    string asimov_data_name = "asimovData";
    ws->Print("v");

    // setup other POIs before set constant
    auto pois = const_cast<RooArgSet*>(mc->GetParametersOfInterest());
    RooStatsHelper::setOtherPOIs(pois, muName+fix_variables, other_poi_value, is_other_poi_const);

    // TIter itr_nuis(mc->GetNuisanceParameters()->createIterator());
    // RooRealVar* tmp_obj;
    // while ( (tmp_obj = (RooRealVar*) itr_nuis()) ){
    //     tmp_obj->setConstant(true);
    // }

    // set nuisance and global to be zero
    
    RooRealVar* tmp_obj;
    if ( ! ws->loadSnapshot("nominalNuis") ){
      TIter itr_nuis(mc->GetNuisanceParameters()->createIterator());
      while ( (tmp_obj = (RooRealVar*) itr_nuis()) ){
         TString tmp_str(tmp_obj->GetName());
         tmp_str.ToLower();
         if( tmp_str.Contains("bin") ){
            tmp_obj->setVal(1.0);
         }else {
            tmp_obj->setVal(0.0);
            tmp_obj->setError(1.0);
         }
      }
    }
    if (! ws->loadSnapshot("nominalGlobs") ){
      TIter itr_glob(mc->GetGlobalObservables()->createIterator());
      while ( (tmp_obj = (RooRealVar*) itr_glob()) ){
         TString tmp_str(tmp_obj->GetName());
         tmp_str.ToLower();
         if(tmp_str.Contains("bin")) continue;
         tmp_obj->setVal(0.0);
      }
    }
    

    // setup constants...
    RooStatsHelper::fixVariables(ws, fix_variables, mc);
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(strategy);

    // auto simPdf =dynamic_cast<RooSimultaneous*>(mc->GetPdf());
    // Check Nuisance parameters
    // RooStatsHelper::CheckNuisPdfConstraint(nuisances, simPdf->getAllConstraints(*observables, *const_cast<RooArgSet*>(nuisances), false));

    stringstream out_ss;

    double obs_p0 = -1;
    double exp_p0 = -1;
    
    RooRealVar* mu_ZZ = ws->var("muZZ");
    if (mu_ZZ){
        mu_ZZ->setVal(1.);
        mu_ZZ->setConstant(false);
    }

    RooRealVar* mu_ZZ1 = ws->var("muZZ_VBF");
    if (mu_ZZ1){
        mu_ZZ1->setVal(1.);
        mu_ZZ1->setConstant(false);
    }

    RooRealVar* mu_ZZ2 = ws->var("muZZ_ggF");
    if (mu_ZZ2){
        mu_ZZ2->setVal(1.);
        mu_ZZ2->setConstant(false);
    }

    RooRealVar* mu_ZZ3 = ws->var("muZZ_rest");
    if (mu_ZZ3){
        mu_ZZ3->setVal(1.);
        mu_ZZ3->setConstant(false);
    }
   
    RooRealVar* mu_llll_ZZ = ws->var("mu_llll_ZZ");
    if (mu_llll_ZZ){
        mu_llll_ZZ->setVal(1.);
        mu_llll_ZZ->setConstant(false);
    }

    RooRealVar* mu_llll_ZZ_VBF = ws->var("mu_llll_ZZ_VBF");
    if (mu_llll_ZZ_VBF){
        mu_llll_ZZ_VBF->setVal(1.);
        mu_llll_ZZ_VBF->setConstant(false);
    }

    RooRealVar* mu_llll_ZZ_ggF = ws->var("mu_llll_ZZ_ggF");
    if (mu_llll_ZZ_ggF){
        mu_llll_ZZ_ggF->setVal(1.);
        mu_llll_ZZ_ggF->setConstant(false);
    }

    RooRealVar* mu_llll_ZZ_rest = ws->var("mu_llll_ZZ_rest");
    if (mu_llll_ZZ_rest){
        mu_llll_ZZ_rest->setVal(1.);
        mu_llll_ZZ_rest->setConstant(false);
    }

    RooRealVar* mu_llvv_ZZ = ws->var("mu_llvv_ZZ");
    if (mu_llvv_ZZ){
        mu_llvv_ZZ->setVal(1.);
        mu_llvv_ZZ->setConstant(false);
    }

    RooRealVar* XS_ggF = ws->var("XS_ggF");
    if(XS_ggF) XS_ggF->setMin(0.);
    RooRealVar* XS_VBF = ws->var("XS_VBF");
    if(XS_VBF) XS_VBF->setMin(0.);
    RooRealVar* mu_ggF = ws->var("mu_ggF");
    if(mu_ggF) mu_ggF->setMin(0.);
    RooRealVar* mu_VBF = ws->var("mu_VBF");
    if(mu_VBF) mu_VBF->setMin(0.);

    if (options == "" || options.find("pvalue") != string::npos) 
    {
        std::cout<<"running pvalue code"<<std::endl;
        if (data_opt == "" || data_opt.find("obs") != string::npos) {
            obs_p0 = RooStatsHelper::getPvalue(ws, mc, obs_data, poi->GetName());
        }
        if (data_opt == "" || data_opt.find("exp") != string::npos) {
            bool do_profile = false; // hard coded
            auto asimov_data = (RooDataSet*) ws->data(asimov_data_name.c_str());
            if(!asimov_data) {
                // set other POIs to other_poi_value and fix
                RooStatsHelper::setOtherPOIs(pois, muName+fix_variables, other_poi_value, true);
                asimov_data = RooStatsHelper::makeAsimovData(ws, 1.0, 0.0,
                        poi->GetName(), mcName.c_str(), dataName.c_str(), do_profile);
                RooStatsHelper::setOtherPOIs(pois, muName+fix_variables, other_poi_value, is_other_poi_const);
            }
            if(asimov_data)
                exp_p0 = RooStatsHelper::getPvalue(ws, mc, asimov_data, poi->GetName());
        }
        cout << "expected p0: " << exp_p0 << endl;
        cout << "obs p0: " << obs_p0 << endl;
    }
    out_ss << fix_variables << " " << obs_p0 << " " << exp_p0 << " " ;

    if (options.find("limit") != string::npos){
        std::cout<<"running limit code"<<std::endl;
        RooRealVar* brw_kappa = ws->var("RBW_kappa");
        if (brw_kappa){

            // when it is LWA workspace, need recursively calculate the limit
            stringstream out_ss_temp;
            float old_exp_limt = 1000.;
            brw_kappa->setVal(0.5);
            brw_kappa->setConstant();

            TString exp_limit;
            Limit::run_limit(ws, mc, obs_data, poi, asimov_data_name.c_str(), &out_ss_temp);
            for (int ii=0; ii<3; ii++) out_ss_temp>>exp_limit;
            float curr_exp_limit = exp_limit.Atof();

            while( fabs((curr_exp_limit-old_exp_limt)/curr_exp_limit) < 0.01 )
            {
                old_exp_limt = curr_exp_limit;
                brw_kappa->setVal(curr_exp_limit/ws->function("SM_xs")->getVal());
                out_ss_temp.str("");
                Limit::run_limit(ws, mc, obs_data, poi, asimov_data_name.c_str(), &out_ss_temp);
                for (int ii=0; ii<3; ii++) out_ss_temp>>exp_limit;
                curr_exp_limit = exp_limit.Atof();
            }
            out_ss << out_ss_temp.str();

        }else{
            Limit::run_limit(ws, mc, obs_data, poi, asimov_data_name.c_str(), &out_ss);
        }
    } else {
        out_ss << endl;
    }

    fstream file_out("stats_results.txt",  fstream::app);

    file_out << out_ss.str();
    file_out.close();
    return 0;
}
