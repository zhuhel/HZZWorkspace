/*
 * =====================================================================================
 *
 *    Description:  toy studies for global pvalue
 *
 *        Version:  1.0
 *        Created:  09/04/2017 12:51:49 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *
 * =====================================================================================
 Input:
 - Workspace ROOT file name
 - POI name
 - number of toys
 - optional:
    - initial value (seed) for randomizer
    - Variables in a list to be fixed to supplied values [var=value,var=value...]

Intended to evaluate global p-value for an excess (High mass search)
------------------------------------------------------------------------------
 */
#include <stdlib.h>
#include <string>
#include <map>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include "HZZWorkspace/RooStatsHelper.h"
#include "RooMinimizer.h"

struct FitResult_t{
    int status; // fitted status
    double mH;  // which mass
    double sigma; // what signficance
    double poi;
    double nll;
};

bool descend_on_sigma(FitResult_t& res1, FitResult_t& res2){
    return res1.sigma > res2.sigma;
}

using namespace std;
int main(int argc, char** argv)
{

    if ((argc > 1 && string(argv[1]) == "help") ||
            argc < 4) {
        cout << argv[0] << " combined.root poi_name ntoys [seed=0] [var=value,var=value]" << endl;
        return 1;
    }

    const string input_name(argv[1]);
    const string poi_name(argv[2]);
    int ntoys = atoi(argv[3]);
    int seed_init=0;
    if (argc > 4) {
      seed_init = atoi(argv[4]);
    }
    string fix_variables = "";
    if (argc > 5) {
        fix_variables = argv[5];
    }

    const string wsName("combined");
    const string mcName("ModelConfig");
    const string dataName("obsData");

    const string out_dir("./");

    const string out_name(Form("toys_%s_seed%d.root", poi_name.c_str(), seed_init)); //don't change!

    cout<<"input = "<< input_name <<std::endl;
    cout<<"poiName = "<< poi_name <<std::endl;
    cout<<"ntoys = "<< ntoys <<std::endl;
    cout<<"input seed = "<< seed_init <<std::endl;
    cout<<" Fix variables: " << fix_variables << endl;
    cout<<"outName " << out_name << endl;

    auto file_in = TFile::Open(input_name.c_str());
    auto workspace = (RooWorkspace*) file_in->Get(wsName.c_str());
    auto mc = dynamic_cast<RooStats::ModelConfig*>(workspace->obj(mcName.c_str()));
    auto obs_data = dynamic_cast<RooDataSet*>(workspace->obj(dataName.c_str()));
    auto poi_var = dynamic_cast<RooRealVar*>(workspace->var(poi_name.c_str()));
    auto mH_var = dynamic_cast<RooRealVar*>(workspace->var("mH"));
    poi_var->setVal(0.);

    // set some variables constants as necessary
    RooStatsHelper::fixVariables(workspace, fix_variables, mc);

    // save nominal snapshot nominalNP,nominalGO
    workspace->saveSnapshot("nominalGO",*mc->GetGlobalObservables());
    workspace->saveSnapshot("nominalNP",*mc->GetNuisanceParameters());

    RooStatsHelper::setDefaultMinimize();
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // condition NPs on observed data
    TStopwatch timer = TStopwatch();
    timer.Start();
    if (false && obs_data){
        // perform a background-only fit and save the best fitted nuisance parameters
        // they will be used in generating toys.
        // If observed data not available, just use nominal values
        poi_var->setVal(0.);
        poi_var->setConstant(true);
        RooNLLVar* nll_ = RooStatsHelper::createNLL(obs_data, mc);
        RooStatsHelper::minimize(nll_);
        delete nll_;
    }
    workspace->saveSnapshot("condGO",*mc->GetGlobalObservables());
    workspace->saveSnapshot("condNP",*mc->GetNuisanceParameters());

    timer.Stop();
    double minutes = timer.RealTime()/60.;
    cout << "[Timer] " << minutes << " min for conditional fit" << endl;
    timer.Reset();

    // prepare tree for output
    auto fout = TFile::Open(Form("%s/%s",out_dir.c_str(),out_name.c_str()), "recreate");
    TTree* physics = new TTree("physics", "physics");
    map<string, double> res;
    res["nll_B_fixed"] = -1;
    res["poi_B_fixed"] = -1;
    res["status_B_fixed"] = -1;

    res["poi_s1"] = -1;
    res["mH_s1"] = -1;
    res["sigma_s1"] = -1;
    res["status_s1"] = -1;
    res["nll_s1"] = -1;

    res["poi_s2"] = -1;
    res["mH_s2"] = -1;
    res["sigma_s2"] = -1;
    res["status_s2"] = -1;
    res["nll_s2"] = -1;
    for( auto& dic : res){
        physics->Branch(dic.first.c_str(),
                &(dic.second), Form("%s/D",dic.first.c_str()));
    }

    // generate mass points for scan according to mass resolution.
    vector<double> masses;
    for(double mH_ = 200; mH_ < 300; mH_ += 2) masses.push_back(mH_);
    for(double mH_ = 300; mH_ < 500; mH_ += 4) masses.push_back(mH_);
    for(double mH_ = 500; mH_ < 700; mH_ += 6) masses.push_back(mH_);
    for(double mH_ = 700; mH_ < 900; mH_ += 8) masses.push_back(mH_);
    for(double mH_ = 900; mH_ < 1000; mH_ += 10) masses.push_back(mH_);
    for(double mH_ = 1000; mH_ < 1200; mH_ += 13) masses.push_back(mH_);
    cout <<"total mass points: "<< masses.size() << endl;
    mH_var->setRange(120, 2000);
    // RooWorkspace* new_ws = new RooWorkspace("toys_data", "toys_data");


    timer.Start();
    for(int itoy = 0; itoy < ntoys; itoy++)
    {
        std::cout<<"\n on toy "<<itoy<<"/"<<ntoys<<std::endl;
        // generate background toys
        poi_var->setVal(0.);
        poi_var->setConstant(true);
        mH_var->setVal(400); // otherwise fit does not converge.
        mH_var->setConstant(true); // otherwise fit does not converge.
        cout << "MH is at: " << mH_var->getVal() << endl;
        RooDataSet* pseduo_data = dynamic_cast<RooDataSet*>(
                RooStatsHelper::generatePseudoData(workspace, poi_name.c_str(), itoy+seed_init)
                );

        // new_ws->import(*pseduo_data);

        // perform a Background-only fit
        poi_var->setVal(0.);
        poi_var->setConstant(true);
        RooNLLVar* nll_bonly = RooStatsHelper::createNLL(pseduo_data, mc);
        int status_bonly = (RooStatsHelper::minimize(nll_bonly) == NULL)?0:1;
        res["nll_B_fixed"] = nll_bonly->getVal();
        res["status_B_fixed"] = status_bonly * 1.0;
        res["poi_B_fixed"] = poi_var->getVal();

        // calculate p0 for each mH.
        poi_var->setConstant(false);
        vector<FitResult_t> fitResultsVec;
        for(double mH_ : masses){
            mH_var->setVal(mH_);
            mH_var->setConstant(true);

            cout <<" mH: before" << mH_var->getVal() << " " << endl;
            RooNLLVar* nll_ = RooStatsHelper::createNLL(pseduo_data, mc);
            int status_ = (RooStatsHelper::minimize(nll_) == NULL)?0:1;
            double sigma_ = RooStatsHelper::calculateSignificance(nll_->getVal(), res["nll_B_fixed"]);
            cout <<" mH: " << mH_var->getVal() << " " << sigma_ << endl;

            // Save the results.
            FitResult_t fit_results;
            fit_results.status = status_;
            fit_results.mH = mH_;
            fit_results.sigma = sigma_;
            fit_results.poi = poi_var->getVal();
            fit_results.nll = nll_->getVal();
            fitResultsVec.push_back(fit_results);

            delete nll_;
        }
        delete nll_bonly;
        delete pseduo_data;

        // sort results based on the significance
        sort(fitResultsVec.begin(), fitResultsVec.end(), descend_on_sigma);

        // keep the first two mass points
        res["poi_s1"] = fitResultsVec.at(0).poi;
        res["mH_s1"] = fitResultsVec.at(0).mH;
        res["sigma_s1"] = fitResultsVec.at(0).sigma;
        res["status_s1"] = fitResultsVec.at(0).status;
        res["nll_s1"] = fitResultsVec.at(0).nll;

        res["poi_s2"] = fitResultsVec.at(1).poi;
        res["mH_s2"] = fitResultsVec.at(1).mH;
        res["sigma_s2"] = fitResultsVec.at(1).sigma;
        res["status_s2"] = fitResultsVec.at(1).status;
        res["nll_s2"] = fitResultsVec.at(1).nll;

        physics->Fill();
    }
    timer.Stop();
    minutes = timer.RealTime()/60.;
    cout << "[Timer] " << minutes << " min for studying " << ntoys << " toys." << endl;

    fout->cd();
    physics->Write();
    // new_ws->Write();
    fout->Close();
    file_in->Close();
}
