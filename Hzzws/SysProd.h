#ifndef __SYSPROD_H__
#define __SYSPROD_H__
#include <string>
#include <map>
#include <vector>

#include "TString.h"
#include "TH1.h"
#include "RooArgSet.h"

using namespace std;

class SysProd{
   public:
        SysProd(const string& configFile);
        SysProd();
        ~SysProd();

        void readConfig(const string& configFile);
        bool checkConfig();
        void process();
        void getObs(string cat, string &oname, RooArgSet &treeobs);
        TH1* getHist(string  hName, string sample, string cat, vector<string> fname, string& obsname, string cuts, string weightLabel, string smooth);
        void getRho(string sample, string smooth, double &rho);
        void getRhoVec(string sample, string  smooth, double &rho_x, double &rho_y);
        void fillEmptyBins(string cat, TH1 *hist);
        void fillNormFile(const string& normFileName, vector<string> catList, vector<string> outSystHistName, vector<string> outNormSystHistName, vector<double> normalizationValues, vector<double> normalizationValuesW);
        bool checkNP(vector<string> systDirs);
        bool checkNormNP(vector<string> outNormSystWeightName, string fName, string tName);

   private:
        double lumi;
        const string configFile;
        std::string m_outDir;
        map<string, map<string, string> > p_dic;
        map<string, map<string, string> > NP_dic;

};

#endif

