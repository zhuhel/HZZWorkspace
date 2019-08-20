/*
 * =====================================================================================
 *
 *       Filename:  testToyStudies.cxx
 *
 *    Description:  toy stduies for low mass 
 *
 *        Version:  1.0
 *        Created:  11/03/2015 11:41:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *   Organization:  
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


using namespace std;
int main(int argc, char** argv)
{

    if ((argc > 1 && string(argv[1]) == "help") ||
            argc < 6) {
        cout << argv[0] << " combined.root poi_name poi_value mH ntoys seed" << endl;
        return 1;
    }

    const string input_name(argv[1]);
    const string poi_name(argv[2]);
    double poi_val = atof(argv[3]);
    const string wsName("combWS");
    const string mcName("ModelConfig");
    const string dataName("combData");
    // const string fix_variables(argv[4]);
    double mH_val = atof(argv[4]);
    int ntoys = atoi(argv[5]);
    int seed_init=0;
    if (argc>6) {
      seed_init = atoi(argv[6]);
    }
    //const string out_dir("/afs/cern.ch/user/g/gcree/work/public/Analysis/HZZWorkspaces/Trunk/plots/toys/"); //change this as you want
    const string out_dir("./");

    const string out_name(Form("toys_mH%.0f_%s%.4f_seed%d.root", mH_val, poi_name.c_str(), poi_val, seed_init)); //don't change!

    cout<<"input = "<< input_name <<std::endl;
    cout<<"poiName = "<< poi_name <<std::endl;
    cout<<"ntoys = "<< ntoys <<std::endl;
    cout<<"input seed = "<< seed_init <<std::endl;
    seed_init = seed_init+int(mH_val) + int(100*poi_val);
    cout<<"used seed = "<<seed_init<<std::endl;
    // cout<<" Fix variables: " << fix_variables << endl;
    cout<<"mH: " << mH_val<< endl;
    cout<<"outName " << out_name << endl;
    

    auto file_in = TFile::Open(input_name.c_str());
    auto workspace = (RooWorkspace*) file_in->Get(wsName.c_str());
    auto mc = dynamic_cast<RooStats::ModelConfig*>(workspace->obj(mcName.c_str()));
    auto obs_data = dynamic_cast<RooDataSet*>(workspace->obj(dataName.c_str()));
    // RooStatsHelper::fixVariables(workspace, fix_variables, mc);
    if(workspace->var("mH")){
        workspace->var("mH")->setVal(mH_val);
        workspace->var("mH")->setConstant(true);
    }

    // save nominal snapshot nominalNP,nominalGO
    workspace->saveSnapshot("nominalGO",*mc->GetGlobalObservables());
    workspace->saveSnapshot("nominalNP",*mc->GetNuisanceParameters());

    RooStatsHelper::setDefaultMinimize();

    // condition NPs on observed data
    TStopwatch timer = TStopwatch();
    timer.Start();
    if (obs_data){
        // profile to the S+B fit
        unique_ptr<RooNLLVar> nll_SBfixed(RooStatsHelper::createNLL(obs_data, mc));
        RooStatsHelper::minimize(nll_SBfixed.get());
    } else {
      std::cout<<"could not find dataset:"<<dataName<<". Skipping conditional step!"<<std::endl;
    }
    workspace->saveSnapshot("condGO",*mc->GetGlobalObservables());
    workspace->saveSnapshot("condNP",*mc->GetNuisanceParameters());
    workspace->saveSnapshot("condPOI",*mc->GetParametersOfInterest());

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
    for( auto& dic : res){
        physics->Branch(dic.first.c_str(), 
                &(dic.second), Form("%s/D",dic.first.c_str()));
    }

    RooStatsHelper::setDefaultMinimize();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    res["mH"] = mH_val;
    res["mu"] = poi_val;
    
    timer.Start();
    for(int itoy = 0; itoy < ntoys; itoy++)
    {
        std::cout<<"\n\nmH="<<mH_val<<" mu="<< poi_val <<" on toy "<<itoy<<"/"<<ntoys<<std::endl;
        int curr_seed = itoy + seed_init;
        res["seed"] = curr_seed;
        RooStatsHelper::generateToy(workspace, poi_name.c_str(), poi_val,
                curr_seed, res);
        physics->Fill();
    }
    timer.Stop();
    minutes = timer.RealTime()/60.;
    cout << "[Timer] " << minutes << " min for studying " << ntoys << " toys." << endl;

    fout->cd();
    physics->Write();
    fout->Close();
    file_in->Close();
}
