#include "HZZWorkspace/SystematicsManager.h"
#include <fstream>
using namespace std;

//-----------------------------------------------------------------------------
// Managing class to handle systematics
// * Assigns which systematics are used based on NPList given in [main] in
//   the workspace configuration file  
//-----------------------------------------------------------------------------

SystematicsManager::SystematicsManager(){
    all_nps = new std::vector<TString>();
}

SystematicsManager::~SystematicsManager(){
    delete all_nps;
}

SystematicsManager::SystematicsManager(const char* fileName){
    all_nps = new std::vector<TString>();
    log_info("Systematic list: %s",fileName);
    this->readNPs(fileName);
}

void SystematicsManager::readNPs(const char* fileName){
    if(all_nps && all_nps->size() > 1){
        log_info("SystematicsManager HAVE a set of NPs!");
        return;
    }
    ifstream ifile(fileName, ifstream::in);
    while (ifile.good() && !ifile.eof()) {
        char line[512];
        //         ifile.getline(line, 512);
        log_info("reading... %s",line);
        // 	log_info("%i",ifile.good());
        // 	log_info("%i",(line[0] != '#'));
        // 	log_info("%i",!string(line).empty());
        if (ifile.getline(line, 512) && line[0] != '#' && !string(line).empty() && line[0] != '[' ) {
            TString np_name(line);
            np_name.ReplaceAll(" ","");
            all_nps ->push_back(np_name);
            log_info("np name %s",np_name.Data());
        }
    }
    log_info("SystematicsManager reads in %lu NPs",all_nps->size());

}

vector<TString>* SystematicsManager::add_sys(SampleBase* sample){
    if (all_nps->size() < 1) return NULL;
    auto* nps_vec = new vector<TString>();
    for(unsigned int i=0; i < all_nps->size(); i++){
        TString& np = all_nps->at(i);
        bool has_norm = sample->addNormSys(  np );
        bool has_shape = sample->addShapeSys( np );
        if(has_shape || has_norm){
            nps_vec ->push_back( np );
        }
    }
    return nps_vec;
}
