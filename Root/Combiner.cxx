#include "HZZWorkspace/Combiner.h"

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <boost/algorithm/string.hpp>

#include <RooRealVar.h>
#include <RooAddPdf.h>
#include <RooWorkspace.h>
#include <RooCategory.h>
#include "RooStats/ModelConfig.h"

#include "HZZWorkspace/Helper.h"
#include "HZZWorkspace/SampleFactory.h"
#include "HZZWorkspace/CoefficientFactory.h"
#include "HZZWorkspace/CBGauss.h"
#include "HZZWorkspace/SampleKeys.h"
#include "HZZWorkspace/RooStatsHelper.h"

//-----------------------------------------------------------------------------
// Source file for operational class Combiner
// *** THIS IS THE MAIN CODE ***
//
// "Combiner" class -
// This class reads a supplied configuration file and coordinates the creation of
// a RooWorkspace containing the data, MC, and likelihood described in the config.
// https://gitlab.cern.ch/HZZ/HZZSoftware/HZZWorkspace/wikis/Config/Introduction
//-----------------------------------------------------------------------------

// Combiner class constructor
Combiner::Combiner(const char* _name, const char* _configName):
    ws_name_(_name),
    config_name_(_configName)
{
    // Print the output workspace name and the configuration name
    cout<<"workspace name: "<< ws_name_ << endl;
    cout<<"configuration name: " << config_name_ << endl;
    // Initialize select class members
    file_path_ = "./";
    data_chain = nullptr;
    mc_chain = nullptr;
    workspace = nullptr;
    weight_var_name = ("weight");
}

// Combiner class destructor
Combiner::~Combiner(){
    if(workspace) delete workspace;
    if(data_chain) delete data_chain;
    if(mc_chain) delete mc_chain;
    if(sysMan) delete sysMan;
    for(auto& coef : allCoefficients) delete coef.second;
}

// Process the configuration file to populate the Combiner class members
void Combiner::readConfig(const char* configName)
{
    // Create a dictionary
    Helper::readConfig(configName, '=', all_dic);

    Helper::printDic<string>(all_dic);

    // Output workspace always has the internal name "combined"
    // IDEA: replace with output filename, stripping .root
    workspace = new RooWorkspace("combined");

    auto& main_dic = all_dic.at("main");

    ///////////////////////////////////
    // load input dir
    // ////////////////////////////////
    try {
        file_path_ = main_dic.at("fileDir")+"/";
    } catch (const out_of_range& oor) {
        cout << "'fileDir' not specific, look in current directory" << endl;
    }
    std::cout<<"looking for files in "<<file_path_<<std::endl;
    Helper::getInputPath(file_path_);
    Helper::getWorkspace(workspace);

    ////////////////////////////////////
    // load systematic cutoff
    ///////////////////////////////////
    try {
        std::string newcut = main_dic.at("SysCutoff");
        strvec cutoffs;
        Helper::tokenizeString(newcut,',',cutoffs);
        bool gotDefault=false;
        for (auto& s: cutoffs){
            strvec typeval;
            Helper::tokenizeString(s,':',typeval);
            if (typeval.size()==1){
                Helper::getSysCutoff("all",TString(typeval[0].c_str()).Atof());
                gotDefault=true;
            }
            else if (typeval.size()==2){
                log_info("Setting SysCutoff %s to %.6f",typeval[0].c_str(),TString(typeval[1].c_str()).Atof());
                Helper::getSysCutoff(typeval[0],TString(typeval[1].c_str()).Atof());
            }
            else {
                log_err("I don't understand 3 arguments per sys cutoff...");
            }
            if (!gotDefault){
                Helper::getSysCutoff("all",TString(typeval[typeval.size()-1].c_str()).Atof());
                gotDefault=true;
            }
        }
        Helper::getSysCutoff("shape",Helper::getSysCutoff("all"));//sets "shape" to "all", if it wasn't set already
    } catch (const out_of_range& oor){
        cout<<"No systematic cutoff provided - keeping any value!"<<std::endl;
    }

    Helper::getDisconnectedArgs();

    ///////////////////////////////////
    //load coefficients
    ///////////////////////////////////

    std::cout<<"reading in coefficients"<<std::endl;
    if (all_dic.find("coefficients")==all_dic.end()){
        log_err("No coefficients provided! aborting");
        exit(-1);
    }
    for (auto& sample : all_dic.at("coefficients")) {
        std::cout<<"on coef of sample "<<sample.first<<std::endl;

        //format
        // {sampleName} = {coefficient formula} : {{args}}

        auto newcoef = CoefficientFactory::CreateCoefficient(sample.second);

        if (!newcoef){
            log_err("Failed to add coefficient for sample %s! Critical failure: aborting Combiner!", sample.first.c_str()); exit(-1);
        }

        allCoefficients[sample.first] = newcoef;
    }

    cout << "Added " << allCoefficients.size() << " Coefficients" << endl;

    ///////////////////////////////////
    //load systematics
    ///////////////////////////////////
    try {
        string NP_list = main_dic.at("NPlist") ;
        sysMan = new SystematicsManager(Form("%s/%s", file_path_.c_str(), NP_list.c_str()));
    } catch (const out_of_range& oor) {
        sysMan = new SystematicsManager();
        cerr << "NPlist is not defined, meaning no systematics used" << endl;
    }

    ///////////////////////////////////
    // load the data
    ///////////////////////////////////
    try {
        string input_data = main_dic.at("data");
        cout << "you are adding DATA!!" << endl;
        data_chain = Helper::loader(input_data.c_str(), "tree_incl_all");
    } catch (const out_of_range& oor) {
        cout << "no data added! Analysis is blinded." << endl;
    }


    ///////////////////////////////////
    // load the MC! (for "full asimov")
    ///////////////////////////////////
    try {
        string input_mc = main_dic.at("mc");
	strvec mc_weight_config;
	Helper::tokenizeString(input_mc,',',mc_weight_config);
	// first argument must be input mc file
	input_mc = mc_weight_config.size() > 0 ? mc_weight_config[0] : "";
	// loop over additional arguments
	for(unsigned int iarg = 1; iarg < mc_weight_config.size() ; ++iarg){
	  string full_arg (mc_weight_config[iarg]);
	  strvec arg;
	  Helper::tokenizeString(full_arg,':',arg);
	  // in case the user wants to use a different variable as weight than the default (='weight')
	  if(arg.size()==2 && arg[0].compare("weight") == 0){
	    weight_var_name = arg[1];
	  }
	  else if(!full_arg.empty()){
	    log_warn(" did not recognize argument given to 'mc' configuration: '%s' . Allowed arguments are: 'weight : <weightvar>' ",full_arg.c_str());
	  }
	}
        mc_chain = Helper::loader(input_mc.c_str(), "tree_incl_all");
	if(mc_chain){
	  // check if input tree contains weight variable
	  if(!mc_chain->GetListOfBranches()->FindObject(weight_var_name.c_str())){
	    log_err("MC weight var %s not contained in input tree -> no MC will be added",weight_var_name.c_str());
	    mc_chain = 0;
	  }
	  else{
	    rename_map_[weight_var_name] = "weightVar"; //add the weight variable to the workspace
	    workspace->factory("weightVar[-1000,1000]");
	    log_info("you are adding weighted MC!! As weight the minitree variable with name '%s' is used", weight_var_name.c_str());
	  }
	}
    } catch (const out_of_range& oor){
        log_info("no MC was added");
    }

    ///////////////////////////////////
    //add categories
    ///////////////////////////////////
    istringstream iss_cat( main_dic["categories"] );
    string category_name;
    RooCategory channelCat("channelCat", "channelCat");
    map<string, RooAbsPdf*> pdfMap;
    map<string, RooDataSet*> dataMap;
    map<string, RooDataSet*> mcMap;
    int catIndex = 0;
    char delim = ',';


    // Threshold for dropping PDF
    double pdfStatThreshold = -1;
    try {
        string thres = main_dic.at("pdfStatThres");
        cout << "you are going to remove PDF based on the avaliable stats" << endl;
        cout << "limit is "<< thres << endl;
        pdfStatThreshold    = (double) atof(thres.c_str());

    } catch (const out_of_range& oor) {
            pdfStatThreshold = -1;
    }


    while( getline( iss_cat, category_name, delim ))
    {
        boost::algorithm::trim(category_name);
        cout << "====================================" << endl;
        cout <<"On category: "<< category_name << endl;

        string mcsets = findCategoryConfig(category_name, "mcsets");
        if(mcsets == "") {
            log_err("No corresponding MC set! %s will be dropped!", category_name.c_str());
            continue;
        }

        vector<string> mcsets_names;
        Helper::tokenizeString( mcsets, ',', mcsets_names) ;
        Category* category = new Category(category_name);
        category->setStatThreshold(pdfStatThreshold);

        ///////////////////////////////////
        // add observables
        ///////////////////////////////////
        RooArgSet ch_obs_minitree; // observables in mini-tree
        RooArgSet ch_obs_ws;       // observables in workspace
        // Create a string to store the observables
        string obs_str;
        if(all_dic.find("observables") != all_dic.end()) {
            // Look for a global observable definition
            obs_str = findCategoryConfig("observables", category_name);
        } else { // --> if not found, look for per-category definitions
            // Look for a section of the configuration file called "observables"
            // and retrieve the observables defined for each category
            // NB: global observables are not supported at this time.
            obs_str = findCategoryConfig( category_name, "observables");
        }
        // IDEA: enable mixing of global and per-category observables
        bool use_adaptive_binning = false;
        getObservables(obs_str, ch_obs_ws, ch_obs_minitree, use_adaptive_binning);
        category->setObservables(ch_obs_ws);
        cout << ch_obs_ws.getSize() << " observable in " << category_name << endl;
        obs_.add(ch_obs_minitree,true); // add observable in this channel to overall obs set.
        obs_ws_.add(ch_obs_ws,true); // add observable in this channel to overall obs set.
        //
        // Add Samples in listed in mcsets
        // Samples are created on-the-fly
        //
        for(auto& mcset_name : mcsets_names){
            string sampleargs = findCategoryConfig(category_name, mcset_name);
            SampleBase* sample = createSample(mcset_name, sampleargs);
            if(sample) {
                category->addSample(sample, sysMan);
                if(use_adaptive_binning){
                    sample->useAdaptiveBinning();
                }
            } else {
                log_err("Failed to add mcset %s! IGNORE THIS SAMPLE!",mcset_name.c_str());
                continue;
            }
        }

        cout <<"End add Samples: "<< category_name << endl;

        string catName(Form("%sCat",category_name.c_str()));
        channelCat.defineType(catName.c_str(), catIndex++);

        //////////////////////////////////////////////////////////////
        // Category will sum individual sample's pdf and add constraint terms
        //////////////////////////////////////////////////////////////
        RooAbsPdf* final_pdf = category->getPDF();
        string final_pdf_name(final_pdf->GetName());
        // final_pdf->getVal(); // To increast the fitting speed??
        std::cout<<"final_pdf @: "<<final_pdf<<std::endl;
        final_pdf->Print();

        workspace->import(*final_pdf, RooFit::RecycleConflictNodes());

        //////////////////////////////////////////////////////////////
        // Add nuisance parameters and global observables for MC statistic
        //////////////////////////////////////////////////////////////
        TIter next_global(category->getGlobal().createIterator());
        TIter next_nuis(category->getNuisance().createIterator());
        RooAbsReal* global;
        RooAbsReal* nuisance;

        while ((
                    global = (RooAbsReal*) next_global(),
                    nuisance = (RooAbsReal*) next_nuis()
               ))
        {
            nuisanceSet_.addClone(*nuisance);
            globalobsSet_.addClone(*global);
        }


        pdfMap[catName] = workspace->pdf(final_pdf_name.c_str());

        ////////////////////////////////////////////
        //  Add Data and MC (per channel)
        ////////////////////////////////////////////
        addDataChan(dataMap, data_chain, ch_obs_minitree, category_name, false);
        addDataChan(mcMap, mc_chain, ch_obs_minitree, category_name, true);

        cout <<"End category: "<< category_name << endl;
        cout <<"Summary: "<< endl;
        final_pdf->Print();
        cout << "====================================" << endl;
        delete category;
        delete final_pdf;
    }

    cout<<"End loop over categories:"<<endl;

    ////////////////////////////////////////////
    //  Add Data //
    ////////////////////////////////////////////
    TString orignames(""), newnames("");
    for(auto items : rename_map_) {
        if (orignames!="") orignames+=",";
        if (newnames!="") newnames+=",";
        orignames+=items.first;
        newnames+=items.second;
    }

    if(data_chain){
        std::cout<<"adding data in ALL categories"<<std::endl;
        obs_.add(channelCat);
        RooDataSet* obsData = new RooDataSet("obsData", "observed data", obs_, RooFit::Index(channelCat), RooFit::Import(dataMap));
        workspace->import(*obsData, RooFit::RenameVariable(orignames.Data(),newnames.Data()));
        //for (auto& d : dataMap) workspace->import(*d.second, RooFit::RenameVariable(orignames.Data(),newnames.Data()));
        obsData->Print();
    }

    if (mc_chain){
        std::cout<<"adding mc in ALL categories"<<std::endl;
        obs_.add(channelCat,true); //silence because it might already be there if you added data!
        static RooRealVar* wgt = new RooRealVar(weight_var_name.c_str(),weight_var_name.c_str(),-1000,1000);
        obs_.add(*wgt);
        RooDataSet* MC = new RooDataSet("MC","Weighted MC",obs_,RooFit::Index(channelCat),RooFit::Import(mcMap), RooFit::WeightVar(weight_var_name.c_str()));
        workspace->import(*MC, RooFit::RenameVariable(orignames.Data(),newnames.Data()));
        //for (auto& d : mcMap) workspace->import(*d.second, RooFit::RenameVariable(orignames.Data(),newnames.Data()));
        MC->Print();
    }

    ////////////////////////////////////////////
    //  Pdf to workspace
    ////////////////////////////////////////////

    auto* simPdf = new RooSimultaneous("simPdf", "simPdf", pdfMap, channelCat);
    workspace->import(*simPdf, RooFit::RecycleConflictNodes(), RooFit::Silence(), RooFit::RenameVariable(orignames.Data(), newnames.Data()));

    this->configWorkspace(workspace);

    ////////////////////////////////////////////
    //  default POI configuration
    ////////////////////////////////////////////
    if(all_dic.find("DefaultPOIInfo") != all_dic.end())
    {
        for(const auto& defaultPoi: all_dic.at("DefaultPOIInfo"))
        {
            vector<string> poiInfo;
            Helper::tokenizeString(defaultPoi.second, ',', poiInfo);

            if(poiInfo.size() != 3)
            {
                log_err("Poi default does not contain 3 arguemnts. Format is <poiName> = <centralVal>, <lower limit>, <upper limit>: aborting Combiner!"); exit(-1);
            }

            float centralVal    = (double) atof(poiInfo.at(0).c_str());
            float lwrLim        = (double) atof(poiInfo.at(1).c_str());
            float uprLim        = (double) atof(poiInfo.at(2).c_str());

            if(lwrLim > uprLim)
            {
                log_err("Lower limit for poi greater than the upper limit: aborting Combiner!"); exit(-1);
            }
            TString poiName(defaultPoi.first);

            if(workspace->var(poiName))
            {
                cout<<"Setting default value for: "<<poiName<< endl;
                cout<<"central val: "<<centralVal<<" upr lim: "<<uprLim<<" lwr lim: "<<lwrLim<< endl;
                workspace->var(poiName)->setVal(centralVal);
                workspace->var(poiName)->setMin(lwrLim);
                workspace->var(poiName)->setMax(uprLim);
            }
            else
            {
                cout<<"POI not found... Skipping: "<<poiName<< endl;
            }
        }
    }

    if (workspace->obj("ModelConfig")){
      std::cout<<"\n PRINTING POI"<<std::endl;
      ((RooStats::ModelConfig*)workspace->obj("ModelConfig"))->GetParametersOfInterest()->Print("v");
      std::cout<<std::endl;
    }

    //////////////////////////////////////////////
    // Add Asimov
    // //////////////////////////////////////////

    workspace->saveSnapshot("ParamsBeforeAsimovGeneration",workspace->allVars());

    auto mc = (RooStats::ModelConfig*)workspace->obj("ModelConfig");

    for (auto & opt : all_dic){
        if (opt.first.find("asimov")==std::string::npos) continue;

        TString asimovName=opt.first.c_str();
        asimovName.ReplaceAll(" ","");
        asimovName.ReplaceAll("asimov:","");

        std::cout<<"going to make asimov dataset called: "<<asimovName<<std::endl;

        for (auto & par : opt.second){
            RooRealVar* var = workspace->var(par.first.c_str());
            if (var)
                var->setVal( atof(par.second.c_str()) );
            else {
                log_err("Trying to build asimov but received an invalid parameter: %s",par.first.c_str());
                exit (-1);
            }
        }

        auto adata = RooStatsHelper::makeUnconditionalAsimov(workspace, mc , asimovName.Data());
        if (!workspace->data(asimovName.Data())) workspace->import(*adata); // currently done inside makeUnconditionalAsimov() so we don't really have to do this
    }

    workspace->loadSnapshot("ParamsBeforeAsimovGeneration");

    //////////////////////////////////////////
    // Import Disconnected Objects marked for saving
    /////////////////////////////////////////
    workspace->import(Helper::getDisconnectedArgs(), RooFit::RecycleConflictNodes(),  RooFit::Silence());

    //////////////////////////////////////////
    // Done
    /////////////////////////////////////////
    // workspace->addClassDeclImportDir(getenv("HZZWSCODEDIR"));
    // workspace->addClassImplImportDir(getenv("HZZWSCODEDIR"));
    workspace->importClassCode("*", false); //import code for any custom classes not in standard ROOT build
    workspace ->writeToFile(ws_name_);
    std::cout<<"\n\nfinal workspace:\n\n"<<std::endl;
    workspace->Print();
}

Coefficient* Combiner::getCoefficient(string& name)
{
    Coefficient* coef = NULL;
    try{
        coef = allCoefficients.at(name);
    }catch(const out_of_range& oor){
        cerr << "ERROR (Combiner::getCoefficient): Coefficient " << name << " not defined! " << endl;
    }
    return coef;
}

string Combiner::findCategoryConfig(const string& sec_name, const string& key_name)
{
    string token = "";
    try { // Try locating the section name as is own section
        // Parse the section according to the keys (categories)
        token = all_dic.at(sec_name).at(key_name);
    } catch (const out_of_range& orr) {
        try { // Try locating section name in the [main] section (global namespace)
            token = all_dic.at("main").at(key_name);
            //token = all_dic.at( "main" ).at( sec_name );
        } catch (const out_of_range& orr) {
            // Loop over the section names and try to match sec_name
            int ntimes = 0;
            for (auto it : all_dic) {
                if(it.first.find(sec_name) != string::npos){
                    if(it.second.find(key_name) != it.second.end()) {
                        token = it.second.at(key_name);
                    }
                    ntimes ++;
                }
            }
            if(ntimes > 1){
                log_err("%s are found more than twice!", sec_name.c_str());
            }
        }
    }
    return token;
}

void Combiner::configWorkspace(RooWorkspace* ws)
{
    cout << "Configuring the workspace" << endl;

    ws->Print();
    ////////////////////
    // define set
    ////////////////////
    RooArgSet* poiSet = (RooArgSet*)ws->allVars().selectByName(Helper::addPoiName().c_str());

    if (poiSet->getSize()==0)
        log_err("I cannot find POI!");

    //////////////////////////////
    // define nusiance parameters and global observables
    //////////////////////////////

    for(auto& np : *(sysMan->getNP()))
    {
        string nuisanceName(Form("alpha_%s", np.Data()));
        string globalName  (Form("nom_%s", np.Data()));
        if (! ws->var(nuisanceName.c_str()) ) {
            log_err("no nusiance parameter: %s in workspace", nuisanceName.c_str());
            continue;
        }
        if (! ws->var(globalName.c_str())) {
            log_err("no global observable: %s in workspace", globalName.c_str());
            continue;
        }
        nuisanceSet_.addClone( *ws->var(nuisanceName.c_str()) );
        globalobsSet_.addClone( *ws->var(globalName.c_str()) );
    }
    ws->defineSet("obs", obs_ws_);
    ws->defineSet("nuisance", nuisanceSet_);
    ws->defineSet("globalobs", globalobsSet_);
    ws->defineSet("poi", *poiSet);

    ////////////////////
    // make model config
    ////////////////////
    RooStats::ModelConfig modelConfig("ModelConfig","H->4l example");
    modelConfig.SetWorkspace           ( *ws );
    modelConfig.SetPdf(*ws->pdf("simPdf"));
    modelConfig.SetParametersOfInterest( *ws->set("poi") );
    modelConfig.SetObservables         ( *ws->set("obs") );
    modelConfig.SetNuisanceParameters  ( *ws->set("nuisance") );
    modelConfig.SetGlobalObservables   ( *ws->set("globalobs") );
    ws->import(modelConfig);
    ws->saveSnapshot("nominalGlobs", *ws->set("globalobs"));
    ws->saveSnapshot("nominalNuis", *ws->set("nuisance"));

    cout << "end of configuring workspace" << endl;
    return;
}


void Combiner::combine(){
    readConfig(config_name_.c_str());
    std::cout<<"end of Combiner::combine"<<std::endl;
}

void Combiner::getObservables(
        const string& obs_str,
        RooArgSet& obs_ws, RooArgSet& obs_minitree,
        bool& adaptive)
{
    // @format
    // obs1_name_in_minitree:obs1_name_in_ws nbins low hi; obs2_name_in_minitree:obs2_name_in_ws nbins low hi

    if (obs_str=="COUNT"){//keyword - the only exception to the above
       std::cout<<"detected COUNT: only doing counting in this category, introduce dummy observable"<<std::endl;
       auto var_ws = new RooRealVar("dummyObs","dummyObs",0,1);
       var_ws->setBins(1);
       var_ws->Print();
       obs_ws.add(*var_ws);
    }
    vector<string> obs_list;
    Helper::tokenizeString(obs_str, ';', obs_list);
    if (obs_list.size() < 1) {
        log_err("observable input cannot be recognized: %s", obs_str.c_str());
        return ;
    }
    for (auto obs_input : obs_list){
        vector<string> obs_parameters;
        Helper::tokenizeString(obs_input, ',', obs_parameters);
        const string& obs_names = obs_parameters.at(0);
        size_t pos_semi_comma = obs_names.find(':');
        string minitree_name = "";
        string ws_name = "";
        if (pos_semi_comma == string::npos) {
            minitree_name = ws_name = obs_names;
        } else {
            minitree_name = obs_names.substr(0, pos_semi_comma);
            ws_name = obs_names.substr(pos_semi_comma+1, obs_names.size());
        }
        cout <<" MINITREE branch name: " << minitree_name << endl;
        cout <<" Observable name in WS: " << ws_name << endl;
        if(rename_map_.find(minitree_name) == rename_map_.end()) {
            rename_map_[minitree_name] = ws_name;
        }
        int nbins = 100;
        double low_val=9999, hi_val=9999;
        // catch variable binning

        std::vector<double> variableBinsList = getBinningFromObsListStr(obs_input);
        if (variableBinsList.size() ==0) {
            if(obs_parameters.size() == 3){
                low_val = (double) atof(obs_parameters.at(1).c_str());
                hi_val = (double) atof(obs_parameters.at(2).c_str());
                adaptive = true;
            } else if(obs_parameters.size() == 4) {
                nbins = atoi(obs_parameters.at(1).c_str());
                low_val = (double) atof(obs_parameters.at(2).c_str());
                hi_val = (double) atof(obs_parameters.at(3).c_str());
            } else {
                log_err("I don't understand: %s", obs_input.c_str());
                continue;
            }
            cout<<"Low: "<< low_val << ", High: " << hi_val << " bins:" << nbins << endl;
            float value = (low_val+hi_val)/2.;
            auto var_ws = new RooRealVar( ws_name.c_str(), ws_name.c_str(),
                    value, low_val, hi_val
                    );
            var_ws->Print();
            if(!adaptive) var_ws->setBins(nbins);
            obs_ws.add(*var_ws);
            auto var_minitree = new RooRealVar( minitree_name.c_str(), minitree_name.c_str(),
                        value, low_val, hi_val
                    );
            if(!adaptive) var_minitree->setBins(nbins);
            obs_minitree.add(*var_minitree);
            var_minitree->Print();
        }
        else {
            low_val = variableBinsList.front();
            hi_val = variableBinsList.back();
            float value = (low_val+hi_val)/2.;
            std::cout << "Detected var binning" << std::endl;
            auto var_ws = new RooRealVar( ws_name.c_str(), ws_name.c_str(),
                    value, low_val, hi_val
                    );
            RooBinning theBinning(variableBinsList.size()-1, &(variableBinsList.at(0)));
            var_ws->setBinning(theBinning);
            var_ws->Print();
            obs_ws.add(*var_ws);
            auto var_minitree = new RooRealVar( minitree_name.c_str(), minitree_name.c_str(),
                        value, low_val, hi_val
                    );
            var_minitree->setBinning(theBinning);
            var_minitree->Print();
            obs_minitree.add(*var_minitree);
        }
    }
}

std::vector<double>
Combiner::getBinningFromObsListStr(const std::string &ListStr) {

    std::vector<double> binEdges;
    auto found_opening_brace = ListStr.find("{"); 
    auto found_closing_brace = ListStr.find("}"); 
    if(found_opening_brace  == std::string::npos ||
    found_closing_brace == std::string::npos){
        return binEdges;
    } 
    // if we are here, this is a var binning
  // Build the list of categories
    std::string myStr = ListStr.substr(found_opening_brace+1, found_closing_brace - found_opening_brace - 1);
  size_t current, previous = 0;
  do {
    current = myStr.find(",", previous);
    // std::cout << " Trying to stod "
              // << myStr.substr(previous, current - previous) << std::endl;
    binEdges.push_back(std::stod(myStr.substr(previous, current - previous)));
    previous = current + 1;

  } while (current != std::string::npos);
  // std::cout << " BIN EDGES ";
  // for (auto p : binEdges) {
  //   std::cout << p << ", ";
  // }
  // std::cout << std::endl;
  return binEdges;
}

SampleBase* Combiner::createSample(const string& name, const string& sample_input)
{
    // {sampletype} : {{sampleargs}}
    strvec sampleargs;
    Helper::tokenizeString(sample_input, ':', sampleargs);
    if (sampleargs.size()!=2){
      log_err("Trying to create sample but requires exactly two colon-seperated fields \"sampletype : arg0, arg1, arg2 ...\"");
      log_err("I received \"%s\"", sample_input.c_str());
      return NULL;
    }
    std::string sampletype = sampleargs[0];
    std::string copyargs = sampleargs[1];
    Helper::tokenizeString(copyargs, ',', sampleargs);

    auto newsample = SampleFactory::CreateSample(sampletype, sampleargs);

    if (!newsample){
        log_err("Failed to add sample %s!,aborting Combiner!", name.c_str());
        return NULL;
    }

    if(allCoefficients.find(name) != allCoefficients.end()) {
        newsample->addCoefficient(allCoefficients[name]);
    } else {
        log_err("Tried to add coefficient for sample %s but I couldn't find it!",name.c_str());
        return NULL;
    }
    return newsample;
}
void Combiner::addCutVariables(RooArgSet& obs, const string& cut)
{
    if (cut == "") return;
    strvec name_list;
    Helper::getListOfNames(cut, name_list);
    for(auto name : name_list){
        if (!obs.find(name.c_str())){
            std::cout<<__func__<<" adding variable: "<<name<<std::endl;
            obs.add(*(new RooRealVar(name.c_str(), name.c_str(), -1E6, 1E6)));
        }
    }
}

void Combiner::addDataChan( map<string, RooDataSet*>& map, TChain* chain, RooArgSet& obs, const std::string& catname, bool weighted){
  if (!chain) return;

  std::cout<<"adding dataset in category "<<catname<<std::endl;
  string cut = findCategoryConfig("cuts", catname);
  cout << "Cut on the dataset: |" << cut << "|" << endl;
  TTree* cutChain = chain->CopyTree(cut.c_str());
  cout << "Observables in minitree: " << obs.getSize() << endl;
  string data_name(Form("%s_%s",(weighted?"mc":"data"),catname.c_str()));
  RooCmdArg wgt_arg = weighted ? RooFit::WeightVar(weight_var_name.c_str()) : RooCmdArg::none() ;
  RooArgSet thisObs(obs);
  if (weighted){
    if(!chain->GetListOfBranches()->FindObject(weight_var_name.c_str())){
      log_err("MC weight var %s not contained in input tree -> no MC added",weight_var_name.c_str());
      return;
    }
    static RooRealVar* wgt = new RooRealVar(weight_var_name.c_str(),weight_var_name.c_str(),-1000,1000);
    thisObs.add(*wgt);
  }
  std::cout<<"reading data type "<<(weighted?"MC":"data")<<std::endl;
  std::cout<<"observables in dataset:"<<std::endl;
  thisObs.Print();

  RooDataSet* data = new RooDataSet(data_name.c_str(), "data set", thisObs, RooFit::Import(*cutChain), wgt_arg);
  map[catname+"Cat"] = data;

  log_info("Read %s in Category %s:",(weighted?"MC":"Data"),catname.c_str());
  data->Print();

  delete cutChain;
}
