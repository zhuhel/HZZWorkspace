//    Description:  A set of helper functions
//
#ifndef _HZZWS_HELPER_H
#define _HZZWS_HELPER_H

#include "HZZWorkspace/Common.h"

#include <stdio.h>
#include <RooRealVar.h>
#include <RooProduct.h>
#include <RooConstVar.h>
#include <RooGaussian.h>
#include <RooPoisson.h>

#include <TChain.h>
#include <TH1.h>
#include "TStopwatch.h"

using namespace std;
namespace Helper{


    /////////////////////////
    // help to read text file
    // //////////////////////
    bool readConfig(const char* filename, const char delim, strdic& all_dic);
    void readNormTable(const char* filename, dbldic& norm_dic, double multiplier = 1.0);
    void readScaleFile(const char* file_name, map<string, double>& all_dic);

    void tokenizeString(const string& str, char delim, vector<string>& tokens);
    void tokenizeString(const char* str, char delim, vector<string>& tokens);
    void tokenizeString(const char* str, char delim, vector<TString>& tokens);
    tstrvec fileList(const char* pattern, std::map<float, TString>* m=NULL);
    template<typename T>
    void printDic( const map<string, map<string, T> >& all_dic )
    {
        for(auto& kv : all_dic){
            cout << "section: |" << kv.first << "|" << endl;
            for(auto& sec : kv.second){
                cout<< "\t |" << sec.first <<"| = |" << sec.second << "|" << endl;
            }
        }
    }
    void readAcceptancePoly(std::vector<double>& params, const char* prod, const char* chan, const char* sys="Nominal");


    // to have a uniformed name convention for nuisance parameters and global name
    RooRealVar* createNuisanceVar(const char* npName);
    RooRealVar* createGlobalVar(const char* npName);
    RooAbsPdf* createConstraint(const char* npName);
    RooAbsPdf* createMCStatConstraint(const char* npName, RooRealVar* np, RooArgList* globalSet=NULL);

    TChain* loader(const string& inFile_name, const string& chain_name);
    TChain* loader(const vector<string>& inFile_name, const string& chain_name);
    bool IsSafeTH1(TH1* h1);
    bool IsGoodTH1(TH1* h1);
    bool TH1FHasEmptyBin(TH1F* h);
    void printStopwatch(TStopwatch& timer);
    const std::string& getInputPath(std::string i=std::string("./"));
    const std::string& addPoiName(std::string i=std::string(""));
    RooWorkspace* getWorkspace(RooWorkspace* i=NULL);
    float getSysCutoff(std::string type, float setValue=-1);
    RooArgSet& getDisconnectedArgs();

    // convert a raw histogram to a format useful for building a RooHistPdf.
    // binWidthCorr tells the function whether or not to normalise each bin to unit width
    // normalise tells the function to subsequently normalise the entire thing to unit integral.
    TH1* prepareHistoInputForPdf(TH1* input, bool binWidthCorr=true, bool normalise=false);

    void getListOfNames(const string& cut, strvec& name_list, strmap& name_map = DEFAULT_STRMAP);

    void readObservable(const string& str, vector<string>& obs_str, string& branch_name);

    std::string removeSpaces(const std::string& s);
    std::string extractDecay(const std::string& s);

    bool isMathSyntax(const string& str);
    // -1, not a boolean type, 0: false, 1: true
    int isBoolean(const string& str);

    // sort string vector by length
    bool sort_str_len(const std::string A, const std::string B);
}
#endif
