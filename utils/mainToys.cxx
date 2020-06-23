/*
 * =====================================================================================
 *
 *    Description:  Master class for obtain toys limits
 *
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <string>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include "HZZWorkspace/RooStatsHelper.h"

#include "HZZWorkspace/Helper.h"
#include "TMath.h"

void print_help(char* argv)
{
    cout << argv << " action combined.root mass poi_value" << endl;
    cout << "action:  " <<"toys, expected, observed" << endl;
    cout << "options work only if a `=` is added!" << endl;
    cout << "options: " <<"-w= workspace name [combined]" << endl;
    cout << "         " <<"-d= dataset name [obsData]" << endl;
    cout << "         " <<"-p= POI name [SigXsecOverSM]" << endl;
    cout << "         " <<"-n= number of toys[20]" << endl;
    cout << "         " <<"-s= seed" << endl;
    // cout << "         " <<"-v= debug [0/1]" << endl;
}
using namespace std;
int main(int argc, char** argv)
{
    enum Action {
        TOYS = 0,
        EXPECTED,
        OBSERVED,
        NOTHING
    };

    if ((argc > 1 && string(argv[1]) == "help") ||
            argc < 5) {
        print_help(argv[0]);
        return 1;
    }
    const string action(argv[1]);
    Action my_action = NOTHING;
    if(action == "toys"){
        my_action = TOYS;
    } else if (action == "expected"){
        my_action = EXPECTED;
    } else if (action == "observed"){
        my_action = OBSERVED; 
    } else {
        print_help(argv[0]);
        return 1;
    }

    const string input_name(argv[2]);
    double mH_val = atof(argv[3]);
    const string str_poi_val(argv[4]);
    
    // start to loop over options

    string wsName("combined");
    string dataName("obsData");
    string muName("SigXsecOverSM");
    int ntoys = 20;
    int seed_init = 1;
    // bool debug = false;
    
    for(int i=4; i < argc; i++){
        const char* key = strtok(argv[i],"=") ;
        const char* val = strtok(0," ") ;
        if (strcmp(key,"-w")==0) wsName = string(val);
        if (strcmp(key,"-d")==0) dataName = string(val);
        if (strcmp(key,"-p")==0) muName = string(val);
        if (strcmp(key,"-n")==0) ntoys = atoi(val);
        if (strcmp(key,"-s")==0) seed_init = atoi(val);
        // if (strcmp(key,"-v")==0) debug = (bool) atoi(val);
    }

    const string out_dir("./");
    
    // setup output
    string out_name = "test.root";
    double poi = -1;
    if (my_action == TOYS){
        poi = atof(str_poi_val.c_str());
        out_name = string(Form("toys_mH%.0f_%s%.4f_seed%d.root", 
                    mH_val, muName.c_str(), poi, seed_init)); //don't change!
    } else if(my_action == EXPECTED) {
        out_name = string(Form("expected_mH%d_seed%d.root", (int)mH_val, seed_init));
    } else {
        out_name = string(Form("observed_mH%d.root", (int)mH_val));
    }

    cout<<"input = "<< input_name <<std::endl;
    cout<<"poiName = "<< muName <<std::endl;
    cout<<"poi str = "<< str_poi_val <<std::endl;
    cout<<"ntoys = "<< ntoys <<std::endl;
    cout<<"input seed = "<< seed_init <<std::endl;
    // seed_init = seed_init+int(mH_val) + int(100*poi_val);
    cout<<"mH: " << mH_val<< endl;
    cout<<"outName " << out_name << endl;
    

    auto file_in = TFile::Open(input_name.c_str());
    auto workspace = (RooWorkspace*) file_in->Get(wsName.c_str());

    const string mcName("ModelConfig");
    auto mc = dynamic_cast<RooStats::ModelConfig*>(workspace->obj(mcName.c_str()));
    auto obs_data = dynamic_cast<RooDataSet*>(workspace->obj(dataName.c_str()));
    auto pois = const_cast<RooArgSet*>(mc->GetParametersOfInterest());
    RooRealVar* poi_val;
    TIter poi_itr(pois->createIterator());

    if(workspace->var("mH")){
        workspace->var("mH")->setVal(mH_val);
        workspace->var("mH")->setConstant(true);
    }
    // set POIs to be positve
    poi_itr.Reset();
    while( (poi_val = (RooRealVar*)poi_itr()) ){
        poi_val->setRange(0, 10);
    }

    // save nominal snapshot nominalNP,nominalGO
    workspace->saveSnapshot("nominalGO",*mc->GetGlobalObservables());
    workspace->saveSnapshot("nominalNP",*mc->GetNuisanceParameters());

    RooStatsHelper::setDefaultMinimize();

    TStopwatch timer = TStopwatch();
    timer.Start();
    if (obs_data && action != "observed")
    {
        // profile nuisance parameters to observed data
        // that depends on the "action"
        // not for "observed"
        if(my_action == EXPECTED) {
            // profile to background-only 
            poi_itr.Reset();
            while( (poi_val = (RooRealVar*)poi_itr()) ){
                if( string(poi_val->GetName()).find("ZZ") != string::npos) continue;
                poi_val->setVal(0);
                poi_val->setConstant(true);
            }
        }
        if(my_action == TOYS) {
            workspace->var(muName.c_str())->setVal(poi);
            workspace->var(muName.c_str())->setConstant(true);
        }
        unique_ptr<RooNLLVar> nll(RooStatsHelper::createNLL(obs_data, mc));
        RooStatsHelper::minimize(nll.get());
    } else {
        std::cout<<"Skipping conditional step!"<<std::endl;
    }
    workspace->saveSnapshot("condGO",*mc->GetGlobalObservables());
    workspace->saveSnapshot("condNP",*mc->GetNuisanceParameters());
    workspace->saveSnapshot("condPOI",*mc->GetParametersOfInterest());
    pois->Print("v");

    timer.Stop();
    double minutes = timer.RealTime()/60.;
    cout << "[Timer] " << minutes << " min for conditional fit" << endl;
    timer.Reset(); 

    // prepare tree for output
    auto fout = TFile::Open(Form("%s/%s",out_dir.c_str(),out_name.c_str()), "recreate");
    TTree* physics = new TTree("physics", "physics");
    map<string, double> res;
    res["nll_SB_free"] = -1;
    res["poi_SB_free"] = -1;
    res["status_SB_free"] = -1;
    res["nll_SB_fixed"] = -1;
    res["poi_SB_fixed"] = -1;
    res["status_SB_fixed"] = -1;
    res["nll_B_fixed"] = -1;
    res["poi_B_fixed"] = -1;
    res["status_B_fixed"] = -1;
    res["mu"] =0;
    res["mH"] = 0;
    res["seed"] = 0;
    res["qzero"] = 0; // likelihood ratio for background-only asimov data.
    res["bonly_muhat"] = 0; // best-fitted mu for background-only asimov data.
    res["bonly_sbFree_status"] = 0;
    res["bonly_sbFix_status"] = 0;

    for( auto& dic : res){
        physics->Branch(dic.first.c_str(), 
                &(dic.second), Form("%s/D",dic.first.c_str()));
    }

    RooStatsHelper::setDefaultMinimize();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    if(my_action == TOYS) {
        res["mH"] = mH_val;
        res["mu"] = poi;

        timer.Start();
        for(int itoy = 0; itoy < ntoys; itoy++)
        {
            std::cout<<"\n\nmH="<<mH_val<<" mu="<< poi <<" on toy "<<itoy<<"/"<<ntoys<<std::endl;
            int curr_seed = itoy + seed_init;
            res["seed"] = curr_seed;
            RooStatsHelper::generateToy(workspace, muName.c_str(), poi,
                    curr_seed, res);
            physics->Fill();
        }
    } else {
        vector<double> poi_list;
        vector<string> items;
        Helper::tokenizeString(str_poi_val.c_str(), ',', items);
        cout <<"poi input: " << str_poi_val << endl;
        for(auto& str_in : items){
            cout <<"adding: " << str_in << endl;
            poi_list.push_back(atof(str_in.c_str()));
        }
        res["mH"]= mH_val;
        if(my_action == EXPECTED){
            for(int itoy = 0; itoy < ntoys; itoy++){
                int curr_seed = itoy + seed_init;
                // generate background-only pseudo-data
                poi_itr.Reset();
                while( (poi_val = (RooRealVar*)poi_itr())) {
                    if( string(poi_val->GetName()).find("ZZ") != string::npos) continue;
                    poi_val->setVal(0);
                    poi_val->setConstant(false);
                }
                RooDataSet* toyData = dynamic_cast<RooDataSet*>(RooStatsHelper::generatePseudoData(workspace, muName.c_str(), curr_seed));
                // fit them with different mu values
                for (int i = 0; i < (int)poi_list.size(); i++)
                {
                    double poival = poi_list.at(i);
                    res["mu"] = poival;
                    // std::cout<<"\n\n mu="<<poival<<std::endl;
                    RooStatsHelper::fitData(workspace, mcName.c_str(), toyData, muName.c_str(), poival, res);
                    physics->Fill();
                }
                delete toyData;
            }
        } else { // OBSERVED
            // create asimov-data
            poi_itr.Reset();
            while( (poi_val = (RooRealVar*)poi_itr())) {
                if( string(poi_val->GetName()).find("ZZ") != string::npos) continue;
                poi_val->setVal(0);
                poi_val->setConstant(true);
            }
            bool do_profile = false;
            unique_ptr<RooDataSet> asimov_data( RooStatsHelper::makeAsimovData(workspace, 0.0, 0.0, 
                        muName.c_str(), mcName.c_str(), dataName.c_str(), do_profile) );

            auto poi_var = (RooRealVar*) workspace->var(muName.c_str());
            for (int i = 0; i < (int)poi_list.size(); i++)
            {
                double poival = poi_list.at(i);
                res["mu"] = poival;
                // if(debug) std::cout<<"\n\n mu="<<poival<<std::endl;
                RooStatsHelper::fitData(workspace, mcName.c_str(), dataName.c_str(), muName.c_str(), poival, res);

                // Fit asimov-data as well
                poi_var->setVal(poival);
                poi_var->setConstant(1);
                unique_ptr<RooNLLVar> nll_asimov(RooStatsHelper::createNLL(asimov_data.get(), mc));
                res["bonly_sbFix_status"] = (RooStatsHelper::minimize(nll_asimov.get(), workspace))->status();
                double nll_bonly_val = nll_asimov->getVal();

                poi_var->setVal(0);
                poi_var->setConstant(0);
                res["bonly_sbFree_status"] = (RooStatsHelper::minimize(nll_asimov.get(), workspace))->status();
                res["bonly_muhat"] = poi_var->getVal();

                if(poi_var->getVal() < 0){ //do tilde
                    poi_var->setVal(0);
                    poi_var->setConstant(1);
                    RooStatsHelper::minimize(nll_asimov.get(), workspace);
                }
                res["qzero"] = 2.*(nll_bonly_val - nll_asimov->getVal());

                physics->Fill();
            }
        }

    }

    timer.Stop();
    minutes = timer.RealTime()/60.;
    cout << "[Timer] " << minutes << " min for studying " << ntoys << " toys." << endl;

    fout->cd();
    physics->Write();
    fout->Close();
    file_in->Close();
}
