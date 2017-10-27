#ifndef RUNASYMPTOTICSCLSCORRECTBANDS_H
#define RUNASYMPTOTICSCLSCORRECTBANDS_H

#include "RooNLLVar.h"
#include "RooAbsReal.h"
#include <string>
#include <sstream>

using namespace std;
using namespace RooFit;
using namespace RooStats;


namespace Limit{
//main
void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName,
		       const char* modelConfigName,
		       const char* dataName,
		       const char* asimovDataName,
		       string folder,
		       double CL, const char* muName = "mu", const string& fixother = " ");
void run_limit(RooWorkspace* ws, ModelConfig* mc, 
        RooDataSet* data, RooRealVar* firstPOI, 
        const char* asimovDataName = "asimovData_0",
        stringstream* out_ss = NULL); 
        
 //for backwards compatibility
void runAsymptoticsCLs(const char* infile,
		       const char* workspaceName = "combWS",
		       const char* modelConfigName = "ModelConfig",
		       const char* dataName = "combData",
		       const char* asimovDataName = "asimovData_0",
		       const char* conditionalSnapshot = "conditionalGlobs_0",
		       const char* nominalSnapshot = "nominalGlobs",
		       string folder = "test",
		       double CL = 0.95, const char* muName = "mu", const string& fixother =" ");

double getLimit(RooNLLVar* nll, double initial_guess = 0);
double getSigma(RooNLLVar* nll, double mu, double muhat, double& qmu);
double getQmu(RooNLLVar* nll, double mu);
void saveSnapshot(RooNLLVar* nll, double mu);
void loadSnapshot(RooNLLVar* nll, double mu);
void doPredictiveFit(RooNLLVar* nll, double mu1, double m2, double mu);
RooNLLVar* createNLL(RooDataSet* _data);
double getNLL(RooNLLVar* nll);
double findCrossing(double sigma_obs, double sigma, double muhat);
void setMu(double mu);
double getQmu95_brute(double sigma, double mu);
double getQmu95(double sigma, double mu);
double calcCLs(double qmu_tilde, double sigma, double mu);
double calcPmu(double qmu_tilde, double sigma, double mu);
double calcPb(double qmu_tilde, double sigma, double mu);
double calcDerCLs(double qmu, double sigma, double mu);
int minimize(RooNLLVar* nll);
int minimize(RooAbsReal* nll);

RooDataSet* makeAsimovData(bool doConditional, RooNLLVar* conditioning_nll, double mu_val,
        const char* muName, string* mu_str = NULL, string* mu_prof_str = NULL, double mu_val_profile = -999, bool doFit = true);
}
#endif

