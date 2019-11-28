// =========================================================================
//
//    Description:
//
// ==========================================================================
#include "HZZWorkspace/Gauss.h"
#include "HZZWorkspace/RelativisticBW.h"
#include "HZZWorkspace/RelativisticBWInt.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <RooDataHist.h>
#include <TKey.h>
#include <TString.h>
#include <TSystem.h>

#include "RooGlobalFunc.h"
#include "RooArgList.h"
#include "RooRealVar.h"
#include "RooFFTConvPdf.h"
#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include <RooFitExtensions/RooBSpline.h>
#include "RooPolyVar.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "HZZWorkspace/Helper.h"

using namespace RooStats;
using namespace HistFactory;
//-----------------------------------------------------------------------------
// Crystal-Ball + Gaussian (CBG) PDF for high mass analysis
// * NWA
// * LWA, no interference
// * LWA, with interference
// Breit-Wigner convoluted with CBG for LWA
//
// mean bounds currently hardcoded!!!
//-----------------------------------------------------------------------------

Gauss::Gauss(const char* _name,
        const char* _input,
        const char* _shape_sys,
        bool _doConv, // Large width
        bool _doSys,
        bool _add_int // Add interference given large width scenario
        ) : SampleBase(_name),
    inputParameterFile(_input),
    workspace(new RooWorkspace(Form("Gauss_%s", _name))),
    mean(new RooRealVar("mean","mean",-5.,5.)),
    doConv(_doConv),
    doSys(_doSys),
    addInt(_add_int),
    order_(3),
    masses_(new vector<double>()),
    bases_(NULL)
{
    if(doSys){
        shape_sys_names_ = new vector<string>();
        map<float, TString> mean_input;
        map<float, TString> sigma_input;
        Helper::fileList(Form("%s/mean_%s", Helper::getInputPath().c_str(), _shape_sys), &mean_input);
        Helper::fileList(Form("%s/sigma_%s", Helper::getInputPath().c_str(), _shape_sys), &sigma_input);
        this->addSysShape(mean_input, sigma_input);
    } else {
        shape_sys_names_ = NULL;
    }
    Helper::addPoiName(mean->GetName());
}

Gauss::~Gauss()
{
    if(mean) delete mean;
    //delete workspace;
    delete masses_;

    if(doSys) {
        for (size_t i=0; i<shape_mean_sys_.size(); ++i)
            delete shape_mean_sys_[i];

        for (size_t i=0; i<shape_sigma_sys_.size(); ++i)
            delete shape_sigma_sys_[i];
    }

    if(shape_sys_names_) delete shape_sys_names_;
}

bool Gauss::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
{
    SampleBase::setChannel(_obs,_ch_name,with_sys);
    textInputParameterValues.clear(); // clean up and get prepared

    if (doSys){
        std::cout<<"setting channel for shape_mean_sys_"<<std::endl;
        for (auto& s: shape_mean_sys_)
            s->SetChannel(category_name_.c_str());

        std::cout<<"setting channel for shape_sigma_sys_"<<std::endl;
        for (auto& s: shape_sigma_sys_)
            s->SetChannel(category_name_.c_str());

        shape_sys_names_->clear();
    }
    // makeCBGParameterization();
    return true;
}

bool Gauss::setParameters(const map<string, pair<double, double> >& input_parameters)
{
    if( !textInputParameterValues.empty() ){
        log_warn("Gauss: parameters will be replaced");
        textInputParameterValues.clear();
    }

    for(const auto& paras: input_parameters){
        textInputParameterValues[paras.first] = paras.second;
    }
    return true;
}

void Gauss::makeCBGParameterization()
{

    if (textInputParameterValues.empty()) {
        loadTextInputFile();
    }
    //The s parameter
    // RooAbsReal* a0ga = variable("sigma");
    RooAbsReal* sGA = variable("sigma");
    workspace->import(*sGA);
}

bool Gauss::addShapeSys(const TString& npName)
{
    std::cout<<"called addShapeSys "<< npName <<std::endl;
    /* In this specific channel the NP name will be:
     * pdf_name + category_name + parameter_name
     * ATLAS_Signal_ggH_ggF_4mu_13TeV_GA_s_p0
     */
    if(!doSys) return false;
    bool resultCB = false;

    if(npName.Contains(category_name_.c_str()))         //Cb parameterization systematics
    {
        for(const auto& sys_name : *shape_sys_names_){
            if(npName.Contains(sys_name.c_str())){
                resultCB = true;
                break;
            }
        }
    }

    bool resultMean = true;                                  //Mean shape systematics
    if (!resultCB){
        for(const auto& nSysHandler : shape_mean_sys_){
            if(! nSysHandler->AddSys(npName)) resultMean = false;
        }
    }

    bool resultSigma = true;                                  //Sigma shape systematics
    if (!resultCB){
        for(const auto& nSysHandler : shape_sigma_sys_){
            if(! nSysHandler->AddSys(npName)){
                resultSigma = false;
            }
        }
    }

    return resultSigma||resultMean||resultCB;
}

RooAbsPdf* Gauss::getPDF()
{
    makeCBGParameterization();
    RooRealVar& localObs=(RooRealVar&)this->obs_list_[0];

    //Get polynomial functions/variables/whatever from workspace
    RooAbsReal* sGA =     (RooAbsReal*)workspace->obj((base_name_+"_sigma").Data());

    if (!sGA ){
        std::cerr<<"ERROR! in "<<__func__<<" missing GA parameters"<<std::endl;
        return NULL;
    }

    if (doSys){

    }

    //build pdfs
    RooAbsPdf* finalPdf =new RooGaussian( ( base_name_ + "_ga" ).Data(), "Gaussian", localObs, *mean, *sGA );
    // RooAbsPdf* finalPdf=new RooAddPdf( (base_name_  + "_cbga" ).Data(), "Crystal Ball + Gaussian", *tmpcb , *tmpga , *fCB );


    workspace->import(*finalPdf,RooFit::RecycleConflictNodes());
    workspace->pdf(finalPdf->GetName())->getVal();//triggers caching
    workspace->pdf(finalPdf->GetName())->Print();
    return workspace->pdf(finalPdf->GetName());

}


void Gauss::loadTextInputFile()
{

    string param_file_name = Helper::getInputPath() + "/" + inputParameterFile+"_"+category_name_+".txt";

    //Read the whole file
    string fullFileName = param_file_name;

    ifstream inputFile( gSystem->ExpandPathName(fullFileName.c_str()) );
    unsigned int parsedLines = 0;


    if (inputFile.is_open()) {

        while ( inputFile.good() ) {
            //Read a line
            string newLine;
            getline( inputFile, newLine );

            std::vector<string> substr;
            Helper::tokenizeString(newLine,' ',substr);
            string blank("");
            substr.erase(std::remove(substr.begin(),substr.end(),blank), substr.end());
            if (substr.empty()) continue;

            //Parse the line
            if (substr.size()==3){
                parsedLines++;
                string parameterName( substr[0]);
                double parameterValue = atof(substr[1].c_str());
                double parameterError = atof(substr[2].c_str());

                //Uniqueness test
                if ( textInputParameterValues.find( parameterName ) == textInputParameterValues.end() ) {
                    //Store name : value pair
                    textInputParameterValues[ parameterName ] = make_pair(parameterValue,parameterError);
                    std::cout<<"storing "<<parameterName<<" = "<<parameterValue<<" +- "<<parameterError<<std::endl;
                }
                else {
                    cerr << "Doubly-defined parameter " << parameterName << " in file " << fullFileName << endl;
                    exit(1);
                }
            }
            else {
                //Ignore irrelevant line
                std::cout<<"failed to read exactly 3 items"<<std::endl;
                std::cout<<"actually read: "<<substr.size()<<" items"<<std::endl;
                for (auto&s : substr) std::cout<<"token: "<<s<<std::endl;
            }
        }
    }
    else {
        cout<<"problem with file "<<fullFileName<<endl;
        exit(3);
    }

    cout << "Parsed " << parsedLines << " lines from input text file "
        << gSystem->ExpandPathName(fullFileName.c_str()) << endl;
    if (parsedLines==0) {
        exit(3);
    }
}


bool Gauss::readTextInputFile(string ParameterName, std::pair<double,double>& pars)
{

    //std::cout<<"in "<<__func__<<" looking for parameter: "<<ParameterName<<" among options:"<<std::endl;
    //for (auto & t: textInputParameterValues) std::cout<<t.first<<" = "<<t.second.first<<","<<t.second.second<<std::endl;

    //Throw an error if the parameter name is not found
    if ( textInputParameterValues.find( ParameterName ) == textInputParameterValues.end() ) {
        cerr << "Parameter name " << ParameterName << " not found in cached input" << endl;
        return false;
    }
    else {
        //Return the parameter value and report it is used
        pars = textInputParameterValues[ ParameterName ];
        cout << "Using " << ParameterName << " = " << pars.first<<" +/- "<<pars.second << endl;
    }
    return true;
}

RooAbsReal* Gauss::variable(const string& parname)
{
    pair<double,double> pars;

    bool found = readTextInputFile(parname.c_str(), pars);
    if (!found) return NULL;

    double ratio = pars.first == 0? 0:fabs(pars.second/pars.first);

    string name(Form("%s_%s", base_name_.Data(),parname.c_str()));

  
    RooRealVar var( name.c_str(), name.c_str(), pars.first);
    workspace->import(var);

    return dynamic_cast<RooAbsReal*>(workspace->obj(name.c_str()));
}

FlexibleInterpVar* Gauss::flexibleInterpVar(const string& fivName, vector<string>& names,
        vector<double>& lowValues, vector<double>& highValues)
{
    RooArgList variables;

    for (size_t inp=0; inp<names.size(); inp++) {
        RooRealVar* np = Helper::createNuisanceVar(names[inp].c_str());
        variables.add(*np);
    }

    FlexibleInterpVar* fiv = new FlexibleInterpVar(("fiv_"+fivName).c_str(), ("fiv_"+fivName).c_str(),
            variables, 1., lowValues, highValues);
    fiv->setAllInterpCodes(4);

    return fiv;
}


void Gauss::BuildBases()
{
    if(masses_->size() < 1) { throw std::runtime_error("Call AddSample first please!"); }
    bases_ = new RooStats::HistFactory::RooBSplineBases(Form("bases_%s", nickname_.c_str()),
            Form("bases_%s", nickname_.c_str()), order_, *masses_, *mean);
}


void Gauss::addSysShape(float m,
        const string& mean_input, const string& sigma_input)
{

    cout << "Gauss: add sys mass = "<< m << " using: "<< mean_input << " " << sigma_input << endl;
    masses_->push_back(m);
    auto s2 = new SysText(mean_input.c_str());
    shape_mean_sys_.push_back(s2);
    auto s3 = new SysText(sigma_input.c_str());
    shape_sigma_sys_.push_back(s3);
}

void Gauss::addSysShape(
            const map<float, TString>& input_mean,
            const map<float, TString>& input_sigma
            )
{
    for(const auto& items: input_mean) {
        std::string meanfilename = input_mean.at(items.first).Data();
        std::string sigmafilename = input_sigma.at(items.first).Data();
        this->addSysShape(items.first, meanfilename, sigmafilename);
    }
}

RooAbsReal* Gauss::getShapeSys(std::string name)
{

    std::string outputName = category_name_;

    std::vector<SysText*>* sysvec;
    if (name=="mean") {
        sysvec = &shape_mean_sys_;
    } else if (name=="sigma") {
        sysvec = &shape_sigma_sys_;
    } else {
        cout << "shape sys: " << name << " not recognized" << endl;
        return NULL;
    }

    RooArgList bs_fiv_list;

    //loop over masses
    for (unsigned int m(0); m < masses_->size(); ++m) {
        auto fiv =  sysvec->at(m)->GetSys(Form("fiv_shape_%s_%s_%d",name.c_str(), base_name_.Data(), (int)masses_->at(m)));
        if (fiv){
            std::cout<<"adding fiv:"<<std::endl;
            fiv->Print();
            bs_fiv_list.add(*fiv);
        } else {
            log_err("no FIV!");
        }
    }

    if (bs_fiv_list.getSize() > 0){
        if (!bases_) BuildBases();
        const char * bs_fiv_name = Form("%s_shape_%s_bs_fiv",outputName.c_str(), name.c_str());
        auto bs_fiv = new RooStats::HistFactory::RooBSpline(bs_fiv_name, bs_fiv_name, bs_fiv_list, *bases_, RooArgSet());
        return bs_fiv;
    } else {
        return NULL;
    }
}
