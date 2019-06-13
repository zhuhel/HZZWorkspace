#include <stdlib.h>
#include <iostream>
#include <string>
#include <utility>

#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <TFile.h>
#include <TTree.h>


#include "Hzzws/RooStatsHelper.h"
#include "Hzzws/Helper.h"

using namespace std;
int main(int argc, char** argv)
{

    RooStatsHelper::setDefaultMinimize();
    string out_name("scan_mu.root");
    string wsName = "combined";
    string mcName = "ModelConfig";
    string dataName = "obsData";
    string muName = "mu";
    int err_id = 0;
    if((argc > 1 && string(argv[1]).find("help") != string::npos) ||
            argc < 7)
    {
        cout << argv[0] << " combined.root out.root ws_name mu_name data_name mH:100:120:130 mH:125.,mu:1.0" << endl;
        cout << "ws_name: name of RooWorkspace" << endl;
        cout << "mu_name: POI name" << endl;
        cout << "mH:100:120:130, scan mH in range of (120, 130) with 100 bins " << endl;
        cout << "data_name: observed data or asimovData_1_0." << endl;
        cout << "\tAsimovData follow convention: asimovData_1_0, 1 is mu value, 0 is mu value profiled to (-1 means no profile)" << endl;
        cout << "mH:125,mu:1.0, set variable to this value and constant" << endl;
        exit(++err_id);
    }

    string input_name(argv[1]);
    int opt_id = 2;
    if (argc > opt_id) out_name = string(argv[opt_id]);
    opt_id++;
    if (argc > opt_id) wsName = string(argv[opt_id]);
    opt_id++;
    if (argc > opt_id) muName = string(argv[opt_id]);
    opt_id++;
    if (argc > opt_id) dataName = string(argv[opt_id]);
    opt_id++;
    
    // set binning and range for POI
    int nbins = 100;
    double low = 0, hi = 3;
    string scan_var_name = "";
    if (argc > opt_id) {
        string options(argv[opt_id]);
        vector<string> tokens;
        Helper::tokenizeString(options, ':', tokens);
        cout <<"scan variable: " << options << endl;
        if (tokens.size() != 4) {
            cout << "scan variable setting is wrong: " << options << endl;
            cout << "e.g.: mH:100:120:130" << endl;
            exit(++err_id);
        } else {
            scan_var_name = tokens.at(0).c_str();
            nbins = atoi(tokens.at(1).c_str());
            low = (double) atof(tokens.at(2).c_str());
            hi = (double) atof(tokens.at(3).c_str());
        }
    }
    opt_id ++;

    // set range for other parameters
    map<string, double> map_var_range;
    if(argc > opt_id) {
        string options(argv[opt_id]);
        cout << "constant: " << options << endl;
        vector<string> tokens;
        Helper::tokenizeString(options, ',', tokens);
        for(auto iter = tokens.begin(); iter != tokens.end(); iter++){
            string token(*iter);
            vector<string> values;
            Helper::tokenizeString(token, ':', values);
            if (values.size() != 2) {
                cout << "constant is wrong: " << token << endl;
                exit(++err_id);
            }
            map_var_range[values.at(0)] = static_cast<double>(atof(values.at(1).c_str()));
        }
    }
    opt_id ++;

    auto* file_in = TFile::Open(input_name.c_str(), "read");
    auto* workspace = (RooWorkspace*) file_in->Get(wsName.c_str());

    if (map_var_range.size() > 0) {
        for(auto it = map_var_range.begin(); it != map_var_range.end(); it++){
            string var_name(it->first);
            auto par = (RooRealVar*) workspace->var(var_name.c_str());
            auto value = it->second;
            if(!par){
                log_warn("%s does not exist!", var_name.c_str());
            } else {
                log_info("%s set to [%.2f]", var_name.c_str(), value);
                par->setVal(value);
                par->setConstant(1);
            }
        }
    }

    // load snapshot for combined workspace
    if ( ! workspace->loadSnapshot("nominalNuis") ){
        log_warn("nominal nuisance not there");
    }
    if ( ! workspace->loadSnapshot("nominalGlobs") ){
        log_warn("nominal global values not there");
    }

    if (dataName.find("asimov") != string::npos){
        log_info("processing asimov data: %s", dataName.c_str());
        auto data = (RooDataSet*) workspace->data(dataName.c_str());
        if (!data) {
            log_info("%s does not exist in current workspace, make one", dataName.c_str());
            vector<string> tokens;
            Helper::tokenizeString(dataName, '_', tokens);
            double mu_val = 1.0, profileMu = 0.0;
            bool do_profile = false;
            size_t n_opts = tokens.size();
            if(n_opts > 1) mu_val = (double) atof(tokens.at(1).c_str());
            if(n_opts > 2) profileMu = (double) atof(tokens.at(2).c_str());
            if(profileMu < 0) do_profile = false;
            data = RooStatsHelper::makeAsimovData(workspace, mu_val, profileMu, muName.c_str(),
                    mcName.c_str(), dataName.c_str(), do_profile);
            dataName = data->GetName();
        }
    }
    auto poi =  workspace->var(muName.c_str());
    poi->setConstant(false);
    TTree* physics = new TTree("physics", "physics");
    RooStatsHelper::ScanPOI(workspace, dataName,
            scan_var_name.c_str(), nbins, low, hi, physics);

    auto* file_out = TFile::Open(out_name.c_str(), "RECREATE");
    file_out ->cd();
    physics->Write();
    file_out->Close();

    delete physics;
    file_in->Close();
    return 0;
}
