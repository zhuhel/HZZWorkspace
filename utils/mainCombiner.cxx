#include "Hzzws/Combiner.h"
#include <iostream>

///////////////////////////////////
// This is the main function of the workspace framework, used to create a workspace based on your config file.
// For documentation on how to build your config file, please see
// https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsZZRunIIWorkspaces#The_config_file
//
// Usage:
// > mainCombiner config.ini wsname.root


using namespace std;

int main(int argc, char* argv[]){
    string configname = "test.ini";
    if (argc > 1 && string(argv[1]) == "help"){
        cout << argv[0] << " config_file combined.root " << endl;
        return 0;
    }
    string out_name = "combined.root";
    if (argc > 1){
        configname = string(argv[1]);
    }
    if(argc > 2) out_name = string(argv[2]);
    cout << configname << endl;
    Combiner* combiner = new Combiner(out_name.c_str(), configname.c_str());
    combiner->combine();
    delete combiner;
    return 0;
}
