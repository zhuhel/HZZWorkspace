// =========================================================================
//
//    Description:
//
// ==========================================================================
#include "HZZWorkspace/CBGauss.h"
#include "HZZWorkspace/RelativisticBW.h"
#include "HZZWorkspace/RelativisticBWInt.h"
#include "HZZWorkspace/RooGravitonRBWPdf.h"
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
// mH bounds currently hardcoded!!!
//-----------------------------------------------------------------------------

CBGauss::CBGauss(const char* _name,
        const char* _input,
        const char* _shape_sys,
        int _doConv, // Large width
        bool _doSys,
        bool _add_int // Add interference given large width scenario
        ) : SampleBase(_name),
    inputParameterFile(_input),
    workspace(new RooWorkspace(Form("CBGauss_%s", _name))),
    mH(new RooRealVar("mH","mH",180.,2000.)),
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
    Helper::addPoiName(mH->GetName());
}

CBGauss::~CBGauss()
{
    if(mH) delete mH;
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

bool CBGauss::setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys)
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

bool CBGauss::setParameters(const map<string, pair<double, double> >& input_parameters)
{
    if( !textInputParameterValues.empty() ){
        log_warn("CBGauss: parameters will be replaced");
        textInputParameterValues.clear();
    }

    for(const auto& paras: input_parameters){
        textInputParameterValues[paras.first] = paras.second;
    }
    return true;
}

void CBGauss::makeCBGParameterization()
{

    if (textInputParameterValues.empty()) {
        loadTextInputFile();
    }

    //The mu parameter is a polynomial in all channels
    RooAbsReal* a0cb = variable("CB_mu_p0");
    RooAbsReal* a1cb = variable("CB_mu_p1");
    RooAbsReal* a2cb = variable("CB_mu_p2");
    RooAbsReal* a3cb = variable("CB_mu_p3");
    RooAbsReal* a4cb = variable("CB_mu_p4");

    //If we are doing a convolution, then the mH is included as parameter in the truth shape, not here
    if (doConv==1){
      log_info("doing convolution of CBGA with truth.. shifting CB_mu_p1 by 1.00!");
      a1cb = new RooFormulaVar(Form("%s_CB_mu_p1_nomH",base_name_.Data()),"@0-1",RooArgList(*a1cb));
    }

    RooPolyVar meanPoly( ( base_name_ + "_mean" ).Data(), "CB #mu", *mH, RooArgList( *a0cb, *a1cb, *a2cb, *a3cb, *a4cb ) );
    workspace->import(meanPoly);

    //The alpha parameter
    RooAbsReal* a0cba = variable("CB_alpha_p0");
    RooAbsReal* a1cba = variable("CB_alpha_p1");
    RooAbsReal* a2cba = variable("CB_alpha_p2");
    RooAbsReal* a3cba = variable("CB_alpha_p3");
    RooAbsReal* a4cba = variable("CB_alpha_p4");
    RooPolyVar alphaCB( ( base_name_ + "_alphaCB" ).Data(), "CB alpha", *mH, RooArgList( *a0cba, *a1cba, *a2cba, *a3cba, *a4cba ) );
    workspace->import(alphaCB);

    //The n parameter
    RooAbsReal* a0nnCB = variable("CB_nn_p0");
    RooAbsReal* a1nnCB = variable("CB_nn_p1");
    RooAbsReal* a2nnCB = variable("CB_nn_p2");
    RooAbsReal* a3nnCB = variable("CB_nn_p3");
    RooAbsReal* a4nnCB = variable("CB_nn_p4");
    RooPolyVar nCB( ( base_name_ + "_nCB" ).Data(), "CB n", *mH, RooArgList( *a0nnCB, *a1nnCB, *a2nnCB, *a3nnCB, *a4nnCB ) );
    workspace->import(nCB);

    //The s parameter
    RooAbsReal* a0ga = variable("GA_s_p0");
    RooAbsReal* a1ga = variable("GA_s_p1");
    RooAbsReal* a2ga = variable("GA_s_p2");
    RooAbsReal* a3ga = variable("GA_s_p3");
    RooAbsReal* a4ga = variable("GA_s_p4");
    RooPolyVar sGA( ( base_name_ + "_sGA" ).Data(), "CB #mu", *mH, RooArgList( *a0ga, *a1ga, *a2ga, *a3ga, *a4ga ) );
    workspace->import(sGA);

    //The f parameter
    RooAbsReal* f0cb = variable("CB_f_p0");
    RooAbsReal* f1cb = variable("CB_f_p1");
    RooAbsReal* f2cb = variable("CB_f_p2");
    RooAbsReal* f3cb = variable("CB_f_p3");
    RooAbsReal* f4cb = variable("CB_f_p4");
    RooPolyVar fCB( ( base_name_ + "_fCB" ).Data(), "CB f", *mH, RooArgList( *f0cb, *f1cb, *f2cb, *f3cb, *f4cb ) );
    workspace->import(fCB);

    //The sCB parameter
    RooAbsReal* a0cbs = variable("CB_s_p0");
    RooAbsReal* a1cbs = variable("CB_s_p1");
    RooAbsReal* a2cbs = variable("CB_s_p2");
    RooAbsReal* a3cbs = variable("CB_s_p3");
    RooAbsReal* a4cbs = variable("CB_s_p4");
    RooPolyVar sCB( ( base_name_ + "_sCB" ).Data(), "CB #sigma", *mH, RooArgList( *a0cbs, *a1cbs, *a2cbs, *a3cbs, *a4cbs ) );
    workspace->import(sCB);

    // acceptance in truth level
    if(doConv==1){
        RooAbsReal* a0acc = variable("ACC_p0");
        RooAbsReal* a1acc = variable("ACC_p1");
        RooAbsReal* a2acc = variable("ACC_p2");
        RooAbsReal* a3acc = variable("ACC_p3");
        workspace->import(*a0acc);
        workspace->import(*a1acc);
        workspace->import(*a2acc);
        if (a3acc) workspace->import(*a3acc);

        if (addInt){
            // parameters for the amplitude of interference term
            // that's derived from "global" complex poly-nominal
            RooAbsReal* a0int = variable("INT_a0");
            RooAbsReal* a1int = variable("INT_a1");
            RooAbsReal* a2int = variable("INT_a2");
            RooAbsReal* a3int = variable("INT_a3");
            RooAbsReal* a4int = variable("INT_a4");
            RooAbsReal* b0int = variable("INT_b0");
            RooAbsReal* b1int = variable("INT_b1");
            RooAbsReal* b2int = variable("INT_b2");
            RooAbsReal* b3int = variable("INT_b3");
            RooAbsReal* b4int = variable("INT_b4");
            workspace->import(*a0int);
            workspace->import(*a1int);
            workspace->import(*a2int);
            workspace->import(*a3int);
            workspace->import(*a4int);
            workspace->import(*b0int);
            workspace->import(*b1int);
            workspace->import(*b2int);
            workspace->import(*b3int);
            workspace->import(*b4int);
        }
    }
    if(doConv==2){
        RooAbsReal* mp0 = variable("m_RBW_p0");
        RooAbsReal* mp1 = variable("m_RBW_p1");
        RooPolyVar mGB( ( base_name_ + "_m_RBW" ).Data(), "GB m", *mH, RooArgList( *mp0, *mp1 ) );
        workspace->import(mGB);

        RooAbsReal* gp0 = variable("G_RBW_p0");
        RooAbsReal* gp1 = variable("G_RBW_p1");
        RooPolyVar gGB( ( base_name_ + "_G_RBW" ).Data(), "GB g", *mH, RooArgList( *gp0, *gp1 ) );
        workspace->import(gGB);
    }
}

bool CBGauss::addShapeSys(const TString& npName)
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

RooAbsPdf* CBGauss::getPDF()
{
    makeCBGParameterization();
    RooRealVar& m=(RooRealVar&)this->obs_list_[0];

    //Get polynomial functions/variables/whatever from workspace
    RooAbsReal* mean=     (RooAbsReal*)workspace->obj((base_name_+"_mean").Data());
    RooAbsReal* sGA =     (RooAbsReal*)workspace->obj((base_name_+"_sGA").Data());
    RooAbsReal* sCB =     (RooAbsReal*)workspace->obj((base_name_+"_sCB").Data());
    RooAbsReal* alphaCB = (RooAbsReal*)workspace->obj((base_name_+"_alphaCB").Data());
    RooAbsReal* nCB =     (RooAbsReal*)workspace->obj((base_name_+"_nCB").Data());
    RooAbsReal* fCB =     (RooAbsReal*)workspace->obj((base_name_+"_fCB").Data());

    if (!mean || !sGA || !sCB || !alphaCB || !nCB || !fCB){
        std::cerr<<"ERROR! in "<<__func__<<" missing CB+GA parameters"<<std::endl;
        return NULL;
    }

    if (doSys){

        // *********************************
        //Systematics on mean
        //
        RooArgList meanArgs(*mean);

        RooAbsReal* meanSys = getShapeSys("mean");
        if (meanSys) meanArgs.add(*meanSys);

        mean = new RooProduct(Form("%s_withsys", mean->GetName()), "CB #mu w/ sys", meanArgs);

        // *********************************
        //Systematics on sigma (gaus or CB)
        //

        RooAbsReal*& sigma = (category_name_.find("4mu") != std::string::npos ? sGA : sCB); //4mu has sys on sGA, 2mu2e and 4e have sys on sCB
        RooArgList sigmaArgs(*sigma);

        RooAbsReal* sigmaSys = getShapeSys("sigma");
        if (sigmaSys) sigmaArgs.add(*sigmaSys);

        sigma = new RooProduct(Form("%s_withsys", sigma->GetName()), "", sigmaArgs);
    }


    //build pdfs
    RooCBShape* tmpcb=new RooCBShape( ( base_name_ + "_cb" ).Data(), "Crystal Ball", m, *mean, *sCB, *alphaCB, *nCB );
    RooGaussian* tmpga=new RooGaussian( ( base_name_ + "_ga" ).Data(), "Gaussian", m, *mean, *sGA );
    RooAbsPdf* finalPdf=new RooAddPdf( (base_name_  + "_cbga" ).Data(), "Crystal Ball + Gaussian", *tmpcb , *tmpga , *fCB );
   
    //LWA 
    if (doConv==1){
        // get the coefficiency for acceptance
        RooAbsReal* a0acc =   (RooAbsReal*)workspace->obj((base_name_+"_ACC_p0").Data());
        RooAbsReal* a1acc =   (RooAbsReal*)workspace->obj((base_name_+"_ACC_p1").Data());
        RooAbsReal* a2acc =   (RooAbsReal*)workspace->obj((base_name_+"_ACC_p2").Data());
        RooAbsReal* a3acc =   (RooAbsReal*)workspace->obj((base_name_+"_ACC_p3").Data());

        // width
        RooAbsReal* gamma = new RooRealVar("gamma","gamma", 0.004, 200);
        ((RooRealVar*)gamma)->setVal(75);
        RooAbsReal* gamma_rel = new RooFormulaVar("gamma_rel","gamma/mH",RooArgList(*gamma,*mH));
        Helper::addPoiName(gamma->GetName());
        RooAbsPdf* truthPdf = NULL;

        if(addInt){
            RooAbsReal* a0int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_a0").Data());
            RooAbsReal* a1int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_a1").Data());
            RooAbsReal* a2int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_a2").Data());
            RooAbsReal* a3int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_a3").Data());
            RooAbsReal* a4int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_a4").Data());
            RooAbsReal* b0int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_b0").Data());
            RooAbsReal* b1int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_b1").Data());
            RooAbsReal* b2int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_b2").Data());
            RooAbsReal* b3int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_b3").Data());
            RooAbsReal* b4int =   (RooAbsReal*)workspace->obj((base_name_+"_INT_b4").Data());

            RooAbsReal* kappa = new RooRealVar("RBW_kappa","RBW_kappa",1.,0.,100.);
            // RooAbsReal* kappa = new RooRealVar("kappa_ggF","kappa_ggF",1.,0.,100.);
            // workspace->import(*kappa);
            truthPdf = new RelativisticBWInt( (base_name_+"_rbw").Data()," Relativistic BW", m, *mH, *gamma_rel, *a0acc, *a1acc, *a2acc, *a3acc, *a0int, *a1int, *a2int, *a3int, *a4int, *b0int, *b1int, *b2int, *b3int, *b4int, *kappa);
        } else {
            truthPdf = new RelativisticBW( (base_name_+"_rbw").Data()," Relativistic BW", m, *mH, *gamma_rel, *a0acc, *a1acc, *a2acc, *a3acc);
        }

        // RooAbsPdf* truthPdf = new RelativisticBW( (base_name_+"_rbw").Data()," Relativistic BW", m, *mH, *gamma_rel);

        //Create the convolution
        // m.setRange("cache",m.getMin()-1000.,m.getMax()+1000.);
        m.setRange("cache", 180, 2000);
        m.setBins(40000,"cache");
        finalPdf = new RooFFTConvPdf( (base_name_+"_conv").Data(), "CBG x RBW", m, *truthPdf, *finalPdf, 3);
    }
    
    // Graviton
    if (doConv==2){
        // get the coefficiency for acceptance
        RooAbsReal* mGB =   (RooAbsReal*)workspace->obj((base_name_+"_m_RBW").Data());
        RooAbsReal* gGB =   (RooAbsReal*)workspace->obj((base_name_+"_G_RBW").Data());
        if (!mGB || !gGB){
            std::cerr<<"ERROR! in "<<__func__<<" missing Graviton parameters"<<std::endl;
            return NULL;
        }

        RooAbsPdf* truthPdf = NULL;
        truthPdf = new RooGravitonRBWPdf( (base_name_+"_rbwG").Data()," Graviton BW", m, *mGB, *gGB);

        //Create the convolution
        // m.setRange("cache",m.getMin()-1000.,m.getMax()+1000.);
        m.setRange("cache", 180, 2000);
        m.setBins(40000,"cache");
        finalPdf = new RooFFTConvPdf( (base_name_+"_convG").Data(), "Graviton CBG x RBW", m, *truthPdf, *finalPdf);
    }

    workspace->import(*finalPdf,RooFit::RecycleConflictNodes());
    workspace->pdf(finalPdf->GetName())->getVal();//triggers caching
    workspace->pdf(finalPdf->GetName())->Print();
    return workspace->pdf(finalPdf->GetName());

}


void CBGauss::loadTextInputFile()
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


bool CBGauss::readTextInputFile(string ParameterName, std::pair<double,double>& pars)
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

RooAbsReal* CBGauss::variable(const string& parname)
{
    pair<double,double> pars;

    bool found = readTextInputFile(parname.c_str(), pars);
    if (!found) return NULL;

    double ratio = pars.first == 0? 0:fabs(pars.second/pars.first);

    string name(Form("%s_%s", base_name_.Data(),parname.c_str()));

    // ignore the systematics that have less than per-mill effect
    if (doSys && fabs(ratio) > 1e-3)
    {
        shape_sys_names_->push_back(parname);
        RooRealVar var( Form("%s_nom", name.c_str()), Form("nominal %s", parname.c_str()), pars.first);

        vector<string> names;
        vector<double> lowValues;
        vector<double> highValues;

        names.push_back(name);
        if(ratio > 1) {
            lowValues.push_back(1E-6);
        } else {
            lowValues.push_back(1. - ratio);
        }
        if (ratio > 2){
            log_err("The error of %s is over 100%%!: %.4f", name.c_str(), ratio);
        }
        highValues.push_back(1. + ratio);

        FlexibleInterpVar* fiv = flexibleInterpVar(name, names, lowValues, highValues);

        RooProduct prod(name.c_str(), name.c_str(), RooArgList(var,*fiv));

        workspace->import(prod);
    }
    else {
        RooRealVar var( name.c_str(), name.c_str(), pars.first);
        workspace->import(var);
    }

    return (RooProduct*)workspace->obj(name.c_str());
}

FlexibleInterpVar* CBGauss::flexibleInterpVar(const string& fivName, vector<string>& names,
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


void CBGauss::BuildBases()
{
    if(masses_->size() < 1) { throw std::runtime_error("Call AddSample first please!"); }
    bases_ = new RooStats::HistFactory::RooBSplineBases(Form("bases_%s", nickname_.c_str()),
            Form("bases_%s", nickname_.c_str()), order_, *masses_, *mH);
}


void CBGauss::addSysShape(float m,
        const string& mean_input, const string& sigma_input)
{

    cout << "CBGauss: add sys mass = "<< m << " using: "<< mean_input << " " << sigma_input << endl;
    masses_->push_back(m);
    auto s2 = new SysText(mean_input.c_str());
    shape_mean_sys_.push_back(s2);
    auto s3 = new SysText(sigma_input.c_str());
    shape_sigma_sys_.push_back(s3);
}

void CBGauss::addSysShape(
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

RooAbsReal* CBGauss::getShapeSys(std::string name)
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
