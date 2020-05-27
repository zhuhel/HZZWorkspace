#include "HZZWorkspace/Coefficient.h"
#include "HZZWorkspace/IntegralHelper.h"

#include "RooStats/HistFactory/FlexibleInterpVar.h"
#include <RooFitExtensions/RooBSpline.h>
#include "RooProduct.h"
#include "RooAddition.h"
#include "RooFormulaVar.h"
#include "RooWorkspace.h"
#include "RooArgList.h"
#include "RooGlobalFunc.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#include "RooTFnBinding.h"

#include "TSystem.h"
#include "TF1.h"
#include "TF2.h"

#include <memory.h>
#include <utility>
#include <algorithm>

//-----------------------------------------------------------------------------
// Operational class to build a coefficient
// GIVE LINK TO WIKI PAGE
//-----------------------------------------------------------------------------

// =============================
// Constructor
// =============================
Coefficient::Coefficient(strmap& input) :
    m_arglist(new RooArgList()),
    m_args(input)
{
    m_global_sys.clear();
    m_fullname=m_nickname=m_channel="";
    m_builtCoef=NULL;
    m_customCutoff = std::make_pair(false,-1);
}

// ============================
// Set Name
// ============================
void Coefficient::setName(const std::string& name)
{
    m_fullname=name;

    // m_nickname=m_fullname.substr(m_fullname.find_last_of('_')+1);
    // Not safe to do so, for example: ATLAS_Background_qqZZ_all,
    // The nick name would be just *all*, which is confusing and dangerous.

    TString tempFullName(m_fullname);
    tempFullName.ReplaceAll("ATLAS_Signal_","");
    tempFullName.ReplaceAll("ATLAS_Background_","");
    tempFullName.ReplaceAll("ATLAS_","");
    m_nickname= string(tempFullName);
}

// =============================
// Set Channel
// ============================
// MG: intentionally not using with_sys?
bool Coefficient::setChannel(const char* channelName, bool with_sys)
{
    m_global_sys.clear();
    m_channel=channelName;
    m_builtCoef=NULL;
    if(m_arglist->getSize() != 0) {
        //empty out the arg list, we're starting into a new channel
        m_arglist->removeAll();
    }
    for (auto& m: m_sysHandler){
        for (auto& s : m.second)
            if(s.second) delete s.second;
        m.second.clear();
    }
    m_sysHandler.clear();

    if (with_sys && m_args.find("sys")!=m_args.end()){

        strvec sysfiles;
        Helper::tokenizeString(m_args["sys"],',',sysfiles);

        for (size_t i=0; i<sysfiles.size(); ++i){
            log_info("in %s adding sysfile %s",__func__,sysfiles[i].c_str());

            //no wildcard, just one sys handler
            if (!TString(sysfiles[i].c_str()).MaybeWildcard()){
                //only one sys handler to inform
                m_sysHandler[i][0]= new SysText();
                if (m_customCutoff.first) {
                    m_sysHandler[i][0]->SetCutoff(m_customCutoff.second);
                }
                m_sysHandler[i][0]->ReadConfig(Form("%s/%s",Helper::getInputPath().c_str(),sysfiles[i].c_str()));
                m_sysHandler[i][0]->SetChannel(channelName);
            }

            //wildcard, probably need many sys handlers
            else {

                std::cout<<"going to create wildcard-based normalization systematics spline"<<std::endl;

                TString p2 = sysfiles[i].c_str();
                if (!(p2.Contains("(") && p2.Contains(")"))){
                    log_err("Received a systematic factor with wildcards but no variable for bases? expected norm_ggF*_low.txt(mH) for example");
                    return false;
                }

                TString dep = sysfiles[i].substr(sysfiles[i].find_first_of('(')).c_str();
                p2.ReplaceAll(dep,"");
                dep.ReplaceAll("(",""); dep.ReplaceAll(")","");

                //Populate map of sys handlers and bases
                std::map<float, TString> map;
                Helper::fileList(Form("%s/%s",Helper::getInputPath().c_str(),p2.Data()),&map);

                if (map.empty()){
                    log_err("for coefficient %s, failed to produce wildcard-based normalization systematic! continuing without sys...",m_nickname.c_str());
                    return true;
                }

                std::vector<double> basevals;

                //inform all the sys handlers
                for (auto& m: map){
                    basevals.push_back(m.first);
                    m_sysHandler[i][m.first] = new SysText();
                    if (m_customCutoff.first) {
                        m_sysHandler[i][m.first]->SetCutoff(m_customCutoff.second);
                    }
                    m_sysHandler[i][m.first]->ReadConfig(m.second.Data());
                    m_sysHandler[i][m.first]->SetChannel(channelName);
                }


                //Create the dependent var
                float lo = map.begin()->first;
                float hi = map.rbegin()->first;
                m_base_var[i] = new RooRealVar(dep.Data(),dep.Data(),0.5*(hi+lo),lo,hi);
                Helper::addPoiName(dep.Data());

                //Create the bpsline bases
                const char * bsname = Form("bs_bases%lu_%s_%s",i,m_nickname.c_str(),channelName);
                m_bspline_bases[i] = new RooStats::HistFactory::RooBSplineBases(bsname, bsname, 3, basevals, *m_base_var[i]);

            }
        }
    }

    return true;
}

// =============================
// Destructor
// =============================
Coefficient::~Coefficient(){
    delete m_arglist;
    for (auto& m : m_sysHandler)
        for (auto& s : m.second)
            delete s.second;
    for (auto& v : m_base_var)
        delete v.second;
}


// =============================
// getCoefficient : the final function called!
// =============================
RooAbsReal* Coefficient::getCoefficient(std::string customname){

    if (m_builtCoef) {
        log_warn("If you code crashes soon after this message, you probably called getCoeffient, deleted the result, and called getCoefficient again! The pointer returned is unique, until setChannel is called and the args are reset!");
        return m_builtCoef;
    }

    if (!BuildCoefficient()) return NULL;

    std::cout<<"creating final coefficient for "<<m_fullname<<" using args:"<<std::endl;
    m_arglist->Print("v");

    std::string prodName(Form("nTot%s_%s", m_fullname.c_str(),m_channel.c_str()));
    if (customname!="") prodName=customname; //override the name
    m_builtCoef = new RooProduct(prodName.c_str(), prodName.c_str(), *m_arglist);
    return m_builtCoef;
}

bool Coefficient::AddSys(const TString& npName){

    //arguments to global have systematics added later
    if (m_args.find("global")!=m_args.end() &&
            TString(m_args["global"].c_str()).Contains(npName)){
        m_global_sys.push_back(npName);
        return true;
    }

    //arguments to poi have systematics added later
    if (m_args.find("poi")!=m_args.end() &&
            TString(m_args["poi"].c_str()).Contains(npName)){
        return true;
    }

    //logic returns true if any one of the m_sysHandler maps returns all true
    bool r1=false;
    for (auto&m : m_sysHandler){
        bool r2 = true;
        for (auto& s : m.second){
//             log_info("adding sys %s to coefficient reference point %.1f",npName.Data(), s.first);
            r2 = r2 && s.second->AddSys(npName);
        }
        if (r2) r1=true;
    }
    return r1;
}

// ===========================
// Build Coefficient (the main task!)
// ==========================

bool Coefficient::BuildCoefficient()
{

    // Parameter of interest (keyword 'poi')
    if (m_args.find("poi")!=m_args.end()){
        // supported syntax looks like: mu_ggF*(mu_offshell - 1 + sqrt(mu_onShell))
        //
        string& poi_syntax = m_args["poi"]; 
        if( poi_syntax.find("muZZ") != string::npos ){
          if(m_channel.find("rest") != string::npos || m_channel.find("bkg") != string::npos) poi_syntax = "muZZ_rest";
          else if(m_channel.find("ggF") != string::npos) poi_syntax = "muZZ_ggF";
          else if(m_channel.find("VBF") != string::npos) poi_syntax = "muZZ_VBF";
        }

        // first find out the POI names
        strvec products;
        strmap products_map;
        Helper::getListOfNames(poi_syntax, products, products_map);
        // sort by length, in case some parameter names share common characters, eg., mu_ggF_off, mu5, mu, ...
        // otherwise the ReplaceAll would mess it
        std::sort(products.begin(), products.end(), Helper::sort_str_len);
        // remove the sqrt, cos, tan, etc..
        auto list_end = std::remove_if(products.begin(), products.end(), Helper::isMathSyntax);

        poi_names.clear();

        unique_ptr<RooArgList> argRooFormular(new RooArgList());
        for(auto poi_name = products.begin(); poi_name != list_end; poi_name ++){
            auto poi_range = products_map[*poi_name];
            // the poi_name can be a real poi, or a nuisance
            TString p2 = poi_range.c_str();
            if(p2.Length()>0) { // nuisance parameter, fiv
              RooRealVar* factor=NULL;

              if(p2.Contains("/")) { // nuisance parameter, fiv
                strvec vals;
                Helper::tokenizeString(p2.Data(),'/',vals);
                RooArgList listA;

                RooStats::HistFactory::FlexibleInterpVar* fiv_factor=NULL;
                if (vals.size()>=1) {
                  vector<double> lowValues;
                  vector<double> highValues;
                  TString varname = "alpha_"+(*poi_name);
                  factor = new RooRealVar(varname.Data(),varname.Data(),
                    (float) atof(vals.at(2).c_str()),
                    (float) atof(vals.at(1).c_str()));
                  listA.add(*factor);
                  highValues.push_back(atof(1 + vals.at(1).c_str()));
                  lowValues. push_back(atof(1 + vals.at(2).c_str()));
                  fiv_factor = new RooStats::HistFactory::FlexibleInterpVar(
                    (*poi_name).c_str(), (*poi_name).c_str(), listA, 1., lowValues, highValues);
                  argRooFormular->add(*fiv_factor);
                  log_info("Adding fiv: %s", (*poi_name).c_str());
                }
              } else { // constant
                float val = atof(poi_range.c_str());
                factor = new RooRealVar((*poi_name).c_str(), (*poi_name).c_str(), val);
                factor->setConstant(true);
                argRooFormular->add(*factor);
                log_info("Adding const: %s", (*poi_name).c_str());
              }
            } else {
              auto poi = GetPOI(*poi_name);
              argRooFormular->add(*poi);
              poi_names.push_back( *poi_name );
              log_info("Adding poi: %s", (*poi_name).c_str());
            }
        }

        // reformulate the POI syntax, using @0, @1 instead of the name itself.
        // It will be much easier for combination!
        TString final_POI_syntax(poi_syntax);

        // firstly ignore the string between '[' and ']'
        while( final_POI_syntax.Contains("[") ) {
          auto inds = final_POI_syntax.Index("[");
          auto inde = final_POI_syntax.Index("]");
          final_POI_syntax.Remove(inds, inde+1-inds);
        }

        for(int iarg = 0; iarg < argRooFormular->getSize(); iarg++){
            final_POI_syntax.ReplaceAll( argRooFormular->at(iarg)->GetName(), Form("@%d", iarg) );
        }

        string formulaName(Form("poiFunc_%s", m_fullname.c_str()));
        if(m_fullname.find("ZZ") != string::npos){
           if(m_channel.find("rest") != string::npos || m_channel.find("bkg") != string::npos) formulaName += "_rest";
           else if(m_channel.find("ggF") != string::npos) formulaName += "_ggF";
           else if(m_channel.find("VBF") != string::npos) formulaName += "_VBF";
        }
        unique_ptr<RooFormulaVar> poiFormula(
                new RooFormulaVar(
                    formulaName.c_str(), formulaName.c_str(),
                    final_POI_syntax.Data(), *(argRooFormular.get())
                    )
                );
        poiFormula->Print();
        // when .ini file appears in 'factors', assuming POI is defined there...
        bool add_poi = true;
        if( m_args.find("factors") != m_args.end() ){
            TString factor_str(m_args["factors"]);
            if(factor_str.Contains(".ini")){
                add_poi = false;
            }
        }
        if(add_poi) AddFactor(poiFormula.get());
    }

    // Global factors (keyword 'global')
    if (m_args.find("global")!=m_args.end()){
        strvec products;
        Helper::tokenizeString(m_args["global"],'*',products);
        for (auto& p:products){
            auto f = GetGlobalFactor(p);
            AddFactor(f);
        }
    }

    // Generic factors (keyword 'factors')
    if ( m_args.find("factors")!=m_args.end() ) {
        strvec splitbycomma;
        Helper::tokenizeString(m_args["factors"],',',splitbycomma);
        if(splitbycomma.size() < 2){
            log_err("missing inputs! %s", m_args["factors"].c_str());
        }
        strvec products;
        Helper::tokenizeString(splitbycomma[0],'*',products);
        strvec terms;
        Helper::tokenizeString(splitbycomma[0],'+',terms);
        int tc=1;

        RooArgList* addlist = new RooArgList();

        for (auto& t:terms){

            strvec products;
            Helper::tokenizeString(t,'*',products);

            RooArgList* prodlist = new RooArgList();

            for (auto& p:products){

                auto f = GetGenericFactor(p, splitbycomma[1]);
                if (f) {
                    if (terms.size()>1) prodlist->add(*f);
                    AddFactor(f);
                }
            }

            if (terms.size()>1){
                auto prod = new RooProduct(Form("%s_%s_Coef_term%d",m_fullname.c_str(),m_channel.c_str(),tc++),"",*prodlist);
                delete prodlist;
                addlist->add(*prod);
            }
        }

        if (terms.size()>1){
            auto add = new RooAddition(Form("%s_%s_Coef_allTerms",m_fullname.c_str(),m_channel.c_str()),"",*addlist);
            delete addlist;
            AddFactor(add);
        }
    }

    // Systematics (keyword 'sys')
    if (m_args.find("sys")!=m_args.end()){
        auto f = GetSystematicFactor(m_args["sys"]);
        if(f) AddFactor(f);
        // gamma terms for MC stat (for SampleCount)
        auto fmc = GetMCStatFactor(m_args["sys"]);
        if(fmc) AddFactor(fmc);
    }

    return true;
}

RooAbsReal* Coefficient::GetGlobalFactor(const std::string& p){

    TString p2 = p.c_str();
    if (!(p2.Contains("(")&&p2.Contains(")"))){
        log_err("tried to add global factor for coefficient %s, but it doesn't provide a value! Received: %s",m_fullname.c_str(), p2.Data());
        log_err("Was expected a format like L(24.8/0.97/1.05)");
        return NULL;
    }
    TString varname(p2.Data(),p2.First('('));
    p2.ReplaceAll(varname,""); p2.ReplaceAll("(",""); p2.ReplaceAll(")","");

    strvec vals;
    Helper::tokenizeString(p2.Data(),'/',vals);
    RooArgList listA;

    //Add the nominal value for the global factor
    RooRealVar* factor=NULL;
    if (vals.size()>=1) {
        float val = atof(vals[0].c_str());
        factor = new RooRealVar(varname.Data(),varname.Data(),val);
        factor->setConstant(true);
    }
    //Add the systematic value to the fiv
    if (vals.size()==3) {
        for (auto& s: m_sysHandler[0]) {
            s.second->AddGlobalSys(varname.Data(),
                    (float) atof(vals.at(1).c_str()),
                    (float) atof(vals.at(2).c_str())
                    );
            if(find(m_global_sys.begin(), m_global_sys.end(), varname) != m_global_sys.end()){
                if(!s.second->AddSys(varname)){
                    cout <<"cannot add systematic for: " << varname << endl;
                }
            }
        }
    }
    if (vals.empty()||vals.size()==2||vals.size()>3){
        log_err("trying to add global var but received the wrong number of values! received %lu",vals.size());
    }

    return factor;
}

RooAbsReal* Coefficient::GetPOI(const std::string& p){

    TString p2(p);
    float lo(0),hi(1e3),nom(1.); //some intelligent guesses at poi range (always default to 1.00)
    if (p2.Contains("mu_ggF_off")){lo=-30.;hi=30.;nom=1.;}
    else if (p2.Contains("mu_VBF_off")){lo=-30.;hi=30.;nom=1.;}
    else if (p2.Contains("mu")){lo=-30.;hi=30.;nom=1.;}
    else if (p2.Contains("r_ggF")){lo=0.;hi=30.;nom=1.;}
    else if (p2.Contains("r_VBF")){lo=0.;hi=30.;nom=1.;}
    else if (p2.Contains("r")){lo=0.;hi=30.;nom=1.;}
    if (p2.Contains("br")||p2.Contains("BR")) {lo=0.;hi=1.;nom=1.;}
    if (p2.Contains("xs")||p2.Contains("XS")) {lo=0.;hi=100;nom=1.;}
    if (p2.Contains("N")) {lo=0.;hi=5000;nom=1.;}

    auto poi = new RooRealVar(p.c_str(),p.c_str(),nom,lo,hi);
    Helper::addPoiName(poi->GetName());
    poi->setConstant(true);
    return poi;
}

RooAbsReal* Coefficient::GetGenericFactor(const std::string& p, const std::string& a)
{
    if (TString(a.c_str()).Contains(".ini")) return GetGenericFactorUsingConfigDic(p,a);
    else if (TString(a.c_str()).Contains(".txt")) return GetGenericFactorUsingNormDic(p,a);
    else return NULL;
}

RooAbsReal* Coefficient::GetGenericFactorUsingConfigDic(const std::string& p, const std::string& a)
{
    log_info("AAA inside %s",__func__);
    //Get the dictionary in name [p] in the .ini file
    strdic dict;
    Helper::readConfig((Helper::getInputPath()+a).c_str(),'=',dict);
    if (dict.find(p)==dict.end()){
        log_err("Received non-matching factor and dictionary (.ini file) - must be configured under [factorname]");
        return NULL;
    }

    strmap facmap = dict[p];
    std::map<std::string, RooAbsReal*> depmap;

    //top-level args are "formula" and "data" (stores actual numbers)
    TString formula = facmap["formula"];
    std::string normDicFile = facmap["data"];

    log_info("AAA found initial formula: %s",formula.Data());
    log_info("AAA found data location: %s",normDicFile.c_str());

    //For each sub-top-level formula in the config, create the needed objects
    for (auto& arg : facmap){
        if (arg.first=="formula" || arg.first=="data") continue;
        log_info("AAA creating subformula for: %s : %s",(arg.first).c_str(),(arg.second).c_str());
        TString subname = Form("%s_%s",p.c_str(),arg.first.c_str());
        TString subformula = arg.second.c_str();

        strvec dependents;
        Helper::getListOfNames(arg.second,dependents);

        RooArgList subargs;

        //for each sub-sub-top-level item, create a POI, or re-use one that exists
        for (auto& arg : dependents){
            log_info("AAA want to add sub-sub dependent: %s",arg.c_str());
            // if (arg=="sqrt" || arg=="cos" || arg=="sin" || arg=="tan" || arg=="log" || arg=="exp") continue;
            if ( Helper::isMathSyntax(arg) ) continue;
            if (Helper::getWorkspace()->function(arg.c_str())==NULL){
                auto v = new RooRealVar(arg.c_str(),arg.c_str(),1.,-9e9,9e9);
                subargs.add(*v);
                Helper::addPoiName(arg.c_str());
                Helper::getWorkspace()->import(*v);
            }
            else
                subargs.add(*Helper::getWorkspace()->function(arg.c_str()));
        }
        RooFormulaVar* subform = new RooFormulaVar(subname.Data(),subformula.Data(),subargs);
        subform->Print();
        Helper::getWorkspace()->import(*subform,RooFit::RecycleConflictNodes());

        formula.ReplaceAll(arg.first.c_str(),subname.Data());
        depmap[arg.first] = subform;
    }

    TString augformula = formula;
    log_info("AAA augmented formula: %s",formula.Data());

    //substitute out polynomial arguments
    strmap polynomialArgSubs; // [replacementname] = [origname]
    int k(0);
    while (formula.Contains("pol(")){
        std::string formulaCopy=formula.Data();
        log_info("formula still has polynomial arguments - try to substitute them out");
        size_t l = formulaCopy.find("pol(")+3;
        size_t r = formulaCopy.find(")",l);
        auto subout = formulaCopy.substr(l,r-l+1);
        log_info("going to sub out: %s",subout.c_str());
        formula = formulaCopy.c_str();
        formula.ReplaceAll(subout.c_str(),Form("QWERTY%dQWERTY",k));
        formulaCopy=formula.Data();
        polynomialArgSubs[Form("QWERTY%dQWERTY",k)]=subout;
	k++;
        log_info("formula after this substitution: %s",formulaCopy.c_str());
    }
    log_info("AAA formula after all substitutions: %s",formula.Data());

    //split up factors in master formula
    strvec allfactors;
    Helper::getListOfNames(formula.Data(), allfactors);

    // remove sqrt, cos, tan, etc...
    auto list_end = std::remove_if(allfactors.begin(), allfactors.end(), Helper::isMathSyntax);
    allfactors.erase(list_end, allfactors.end());


    log_info("AAA split this up into:");
    for (auto& x : allfactors)
        log_info("\t%s",x.c_str());

    //substitute back in polynmoial arguments
    for (auto& x: allfactors){
        TString x2 = x.c_str();
        for (auto& r : polynomialArgSubs) x2.ReplaceAll(r.first.c_str(),r.second.c_str());
        x.assign(x2.Data());
    }

    log_info("AAA re-replaced this up into:");
    for (auto& x : allfactors)
        log_info("\t%s",x.c_str());

    //loop over factors and use GetGenericFactorsUsingNormDic to create each
    for (auto& x : allfactors){
        if (depmap.find(x)!=depmap.end()) continue;

        if (std::binary_search(poi_names.begin(), poi_names.end(), x)){
            // If the variable is in POI, just create a RooRealVar
            depmap[x] = new RooRealVar(x.c_str(), x.c_str(), 0, 100);
            continue;
        }

        log_info("AAA now I'm going over to UsingNormDic for factor: %s",x.c_str());
        depmap[x] = GetGenericFactorUsingNormDic(x,normDicFile);
        if (depmap[x]==NULL) return NULL;
    }

    //modify master formula to use new names <oldname>_<channel>
    RooArgList al;
    for (auto& x : depmap){
        if (! std::binary_search(poi_names.begin(), poi_names.end(), x.first)){
            // if not a POI, it must depend on category
            augformula.ReplaceAll(x.first.c_str(), x.second->GetName());
        }
        al.add(*x.second);
    }

    log_info("AAA final formula: %s",augformula.Data());
    log_info("AAA args:");
    al.Print("v");

    //create final formulavar
    RooFormulaVar* form = new RooFormulaVar(Form("%s_%s",p.c_str(),m_channel.c_str()),augformula.Data(),al);

    form->Print();

    return form;
}


RooAbsReal* Coefficient::GetGenericFactorUsingNormDic(const std::string& p, const std::string& a)
{
    std::map< std::string, std::map<std::string, double> > dict;
    Helper::readNormTable((Helper::getInputPath()+a).c_str(), dict);
    RooAbsReal* factor(NULL);
    TString p2 =p;
    //Polynomial
    if (p2.Contains("pol"))
    { //only special case supported

        if (!(p2.Contains("(") && p2.Contains(")"))){
            log_err("Received an coefficient factor that looks like a polynomial but with no argument? skipping..");
            return NULL;
        }
        size_t left_br = p.find_first_of('(');
        size_t right_br = p.find_first_of(')');
        TString deplist(p.substr(left_br, right_br-left_br+1).c_str());
        p2.ReplaceAll(deplist,"");
        deplist.ReplaceAll("(",""); deplist.ReplaceAll(")","");

        log_info("label: %s, parameters: %s", p2.Data(), deplist.Data());

        strvec deps;
        Helper::tokenizeString(deplist.Data(),'|',deps);

        //variables that go into the formula
        RooArgList vars;

        for (auto& d:deps){
            bool added=false;

            //look for a generic object of that name (e.g. another function) in the workspace
            if (Helper::getWorkspace()->function(d.c_str())){
                added=true;
                vars.add(*Helper::getWorkspace()->function(d.c_str()));
            }
            //look for a morphing basis variable in this coefficient
            for (auto& v: m_base_var)
                if (d == v.second->GetName()) {
                    added=true;
                    vars.add(*v.second);
                }
            //otherwise, create it as a poi
            if (!added){
                auto v = new RooRealVar(d.c_str(),d.c_str(),1.,-9e9,9e9);
                vars.add(*v);
                Helper::addPoiName(d.c_str());
            }

        }
        std::cout<<"Attempting to add polynomial coefficient! Hang on tight."<<std::endl;

        TString formula="0.0";
        //Any better way to do this (generically?)
        if (deps.size()==1){
            for (int ox(0);ox<10;++ox){
                TString pname =Form("%s_a%d",p2.Data(),ox);
                if (dict.find(pname.Data())!=dict.end())
                    formula += Form(" + (%e)*(@0)^%d",dict[pname.Data()][m_channel],ox);
            }
        }
        else if (deps.size()==2){
            for (int ox(0);ox<10;++ox) for (int oy(0);oy<10;++oy){
                TString pname =Form("%s_a%d_b%d",p2.Data(),ox,oy);
                if (dict.find(pname.Data())!=dict.end())
                    formula += Form(" + (%e)*(@0)^%d*(@1)^%d",dict[pname.Data()][m_channel],ox,oy);
            }
        }
        else if (deps.size()==3){
            for (int ox(0);ox<10;++ox) for (int oy(0);oy<10;++oy) for (int oz(0);oz<10;++oz){
                TString pname =Form("%s_a%d_b%d_c%d",p2.Data(),ox,oy,oz);
                if (dict.find(pname.Data())!=dict.end())
                    formula += Form(" + (%e)*(@0)^%d*(@1)^%d*(@2)^%d",dict[pname.Data()][m_channel],ox,oy,oz);
            }
        }
        else if (deps.size()==4){
            for (int ox(0);ox<10;++ox) for (int oy(0);oy<10;++oy) for (int oz(0);oz<10;++oz) for (int ow(0);ow<10;++ow){
                TString pname =Form("%s_a%d_b%d_c%d_d%d",p2.Data(),ox,oy,oz,ow);
                if (dict.find(pname.Data())!=dict.end())
                    formula += Form(" + (%e)*(@0)^%d*(@1)^%d*(@2)^%d*(@3)^%d",dict[pname.Data()][m_channel],ox,oy,oz,ow);
            }
        }
        else {
            log_err("TRIED TO CREATE %lu-dimensional polynomial... not supported!",deps.size());
            return NULL;
        }
        for (const int i : {0,1,2,3})
            formula.ReplaceAll(Form("*(@%d)^0",i),"");
        if (formula!="0.0")
            formula.ReplaceAll("0.0 + ","");

        log_info("Created a formula for %s_%s: %s",p2.Data(),m_channel.c_str(),formula.Data());
        // Add support of option: INT
        if( m_args.find("INT") != m_args.end()) // As long as INT is defined, don't care the value.
        {
            log_info("Found INT, consider it's with interference");
            // This part adds a factor that is the ratio of normalization of inteference over the that of signal-only

            RooFormulaVar* factor_raw = new RooFormulaVar(Form("%s_%s_raw",p2.Data(),m_channel.c_str()),formula.Data(),vars);

            float mH_down = 200.0;
            float mH_up = 2000.0;

            RooRealVar* mH = new RooRealVar("mH", "mH", mH_down, mH_up);
            mH->setVal(700.);

            // Relative width value for each mass, used for generating a TF1 map for calculating
            // the cross section! Different name w.r.t the one used in the RooRelativeBW.
            // This value is set priori to the fitting. For each benchmark width 1%, 5%, 10% and 15%,
            // before making the workspaces, one needs to change this value, recompile the code.
            RooAbsReal* gamma_frac = new RooRealVar("gamma_frac","gamma_frac",0.15);

            TF1* f_HiggsXS = new TF1("SM_xs", IntegralHelper::getHiggsXSTF1, mH_down, mH_up, 1);
            f_HiggsXS->SetNpx(2);
            RooArgList list_xs("list_xs");
            list_xs.add(*gamma_frac);
            RooAbsReal *SM_xs = RooFit::bindFunction(f_HiggsXS, *mH, list_xs);


            RooRealVar* xs = new RooRealVar("XS_ggF","XS_ggF",1.,0.00000001,100.);
            RooAbsReal* kappa_ggF = new RooFormulaVar("kappa_ggF","TMath::Sqrt(1 / ( @0*(@0>0.00001*@1)/@1+0.00001*(1-(@0>0.00001*@1))))",RooArgList(*xs,*SM_xs));

            TF1* f_SigInt = new TF1(Form("SignalIntegral_%s",m_channel.c_str()), IntegralHelper::getSignalIntegralTF1, 200., 2000., 5);
            f_SigInt->SetNpx(80);
            TF2* f_TotInt = new TF2(Form("TotalIntegral_%s",m_channel.c_str()), IntegralHelper::getTotalIntegralTF2, 200., 2000., 0.05, 500., 15);
            f_TotInt->SetNpx(80);
            f_TotInt->SetNpy(800);
            // f_TotInt->SetNpx(2);
            // f_TotInt->SetNpy(2);

            RooAbsReal* accp0 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_ACC_p0",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_ACC_p0",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* accp1 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_ACC_p1",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_ACC_p1",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* accp2 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_ACC_p2",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_ACC_p2",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* accp3 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_ACC_p3",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_ACC_p3",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* inta0 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_a0",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_a0",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* inta1 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_a1",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_a1",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* inta2 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_a2",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_a2",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* inta3 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_a3",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_a3",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* inta4 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_a4",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_a4",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* intb0 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_b0",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_b0",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* intb1 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_b1",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_b1",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* intb2 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_b2",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_b2",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* intb3 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_b3",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_b3",m_channel.c_str()),1.,-10000.,10000.);
            RooAbsReal* intb4 = new RooRealVar(Form("ATLAS_Signal_ggF_%s_INT_b4",m_channel.c_str()),Form("ATLAS_Signal_ggF_%s_INT_b4",m_channel.c_str()),1.,-10000.,10000.);
            RooArgList list_sig("list_sig");
            list_sig.add(*gamma_frac);
            list_sig.add(*accp0);
            list_sig.add(*accp1);
            list_sig.add(*accp2);
            list_sig.add(*accp3);
            RooArgList list_int("list_int");
            list_int.add(*gamma_frac);
            //	list_int.add(*kappa_ggF);
            list_int.add(*accp0);
            list_int.add(*accp1);
            list_int.add(*accp2);
            list_int.add(*accp3);
            list_int.add(*inta0);
            list_int.add(*inta1);
            list_int.add(*inta2);
            list_int.add(*inta3);
            list_int.add(*inta4);
            list_int.add(*intb0);
            list_int.add(*intb1);
            list_int.add(*intb2);
            list_int.add(*intb3);
            list_int.add(*intb4);
            RooAbsReal *signal_integral = RooFit::bindFunction(f_SigInt, *mH, list_sig);
            RooAbsReal *total_integral  = RooFit::bindFunction(f_TotInt, *mH, *kappa_ggF, list_int);

            //        Helper::getWorkspace()->import(*xs);
            Helper::getWorkspace()->import(*SM_xs, RooFit::RecycleConflictNodes());
            Helper::getWorkspace()->import(*mH, RooFit::RecycleConflictNodes());
            Helper::getWorkspace()->import(*gamma_frac, RooFit::RecycleConflictNodes());
            Helper::getWorkspace()->import(*kappa_ggF, RooFit::RecycleConflictNodes());

            //Helper::getWorkspace()->var("gamma_frac")->setConstant(true);

            /*        Helper::getWorkspace()->import(*accp0);
                      Helper::getWorkspace()->import(*accp1);
                      Helper::getWorkspace()->import(*accp2);
                      Helper::getWorkspace()->import(*inta0);
                      Helper::getWorkspace()->import(*inta1);
                      Helper::getWorkspace()->import(*inta2);
                      Helper::getWorkspace()->import(*inta3);
                      Helper::getWorkspace()->import(*inta4);
                      Helper::getWorkspace()->import(*intb0);
                      Helper::getWorkspace()->import(*intb1);
                      Helper::getWorkspace()->import(*intb2);
                      Helper::getWorkspace()->import(*intb3);
                      Helper::getWorkspace()->import(*intb4);
                      */
            factor = new RooFormulaVar(Form("%s_%s",p2.Data(),m_channel.c_str()),"(@0)*(@1/@2)",RooArgList(*factor_raw,*total_integral,*signal_integral));
            factor->Print();
        } else {
            factor = new RooFormulaVar(Form("%s_%s",p2.Data(),m_channel.c_str()),formula.Data(),vars);
        }
    } else {
        //Not a polynomial, just a generic factor
        if (dict.find(p)!=dict.end() && dict[p].find(m_channel)!=dict[p].end())
        {
            float val = dict[p][m_channel];
            TString name = (p + "_" + m_channel).c_str();
            factor = new RooRealVar(name.Data(),name.Data(),val);
            ((RooRealVar*)factor)->setConstant(true);
        } else {
            log_err("trying to add a factor %s but couldn't find it in my dictionary. It will be skipped!",p.c_str());
            log_err("the available dictionary was:");
            for (auto& p: dict) std::cout<<p.first<<std::endl;
            if (dict[p].find(m_channel)!=dict[p].end()) for (auto& q : dict[p]) std::cout<< q.first<<" = "<<q.second<<std::endl;
        }
    }
    return factor;
}

RooAbsReal* Coefficient::GetMCStatFactor(const std::string& p)
{
    TString sysarg(p.c_str());

    RooArgList prodset;

    strvec sysfiles;
    Helper::tokenizeString(p,',',sysfiles);

    for (size_t i=0; i<sysfiles.size(); ++i){
        log_info("in %s adding sysfile (%lu/%lu) %s ",__func__,i,sysfiles.size(),sysfiles[i].c_str());

        //Form a single set of systematics
        if (!TString(sysfiles[i].c_str()).MaybeWildcard()){
            // MC stat
            const char* mcstat_name = Form("gamma_stat%lu_%s_%s",i,m_nickname.c_str(),m_channel.c_str());
            auto mcstat = m_sysHandler[i][0]->GetSys(mcstat_name);
            if (!mcstat){
                log_err("unable to get mc stat for this file: %s", sysfiles[i].c_str());
                continue;
            }
            prodset.add(*mcstat);
            gamma_.add(m_sysHandler[i][0]->GetGammas());
        }
     }
    if (prodset.getSize()==0) return NULL;

    const char* fivname = Form("gamma_stat_%s_%s",m_nickname.c_str(),m_channel.c_str());
    return new RooProduct(fivname,fivname,prodset);
}

RooAbsReal* Coefficient::GetSystematicFactor(const std::string& p)
{
    TString sysarg(p.c_str());

    RooArgList prodset;

    strvec sysfiles;
    Helper::tokenizeString(p,',',sysfiles);

    for (size_t i=0; i<sysfiles.size(); ++i){
        log_info("in %s adding sysfile (%lu/%lu) %s ",__func__,i,sysfiles.size(),sysfiles[i].c_str());

        //Form a single set of systematics
        if (!TString(sysfiles[i].c_str()).MaybeWildcard()){
            const char* fivname = Form("fiv%lu_%s_%s",i,m_nickname.c_str(),m_channel.c_str());
            m_sysHandler[i][0]->Print();
            auto fiv = m_sysHandler[i][0]->GetSys(fivname);
            if (!fiv){
                log_err("unable to create FIV for this file: %s", sysfiles[i].c_str());
                continue;
            }
            prodset.add(*fiv);
            np_.add(m_sysHandler[i][0]->GetNPs());
        }

        else{ //For a wildcard (multiple) set of systematics.. make a bspline out of them

            if (!m_bspline_bases[i]){ log_err("coefficient %s has no bspline_bases! no systematic will be added", m_nickname.c_str()); return NULL; }

            RooArgList bs_fiv_list;
            for (auto & s: m_sysHandler[i]) {
                // Loop over all the files that fit the WildCard...

                TString numToString = Form("%.2f",s.first);
                numToString.ReplaceAll(".","p");
                const char* fivname = Form("fiv%lu_%s_%s_%s",i,numToString.Data(),m_nickname.c_str(),m_channel.c_str());

                auto fiv = s.second->GetSys(fivname);
                if (!fiv){
                    log_err("unable to create FIV for %.1f point of spline, all normalization systematics will be skipped for this sample!",s.first);
                    continue;
                }
                bs_fiv_list.add(*fiv);
                np_.add(s.second->GetNPs());
            }
            const char* bsname = Form("bs_fiv%lu_%s_%s",i,m_nickname.c_str(),m_channel.c_str());
            if (bs_fiv_list.getSize() > 0){ // In case of no systematics added
                auto bs_fiv = new RooStats::HistFactory::RooBSpline(bsname, bsname, bs_fiv_list, *m_bspline_bases[i], RooArgSet());
                prodset.add(*bs_fiv);
            } else {
                log_warn("Likely no systematics added in %s channel", m_channel.c_str());
            }
        }
     }
    if (prodset.getSize()==0) return NULL;

    const char* fivname = Form("fiv_%s_%s",m_nickname.c_str(),m_channel.c_str());
    return new RooProduct(fivname,fivname,prodset);
}

bool Coefficient::AddFactor(RooAbsReal* f){
  if (!f) return false;
  log_info("adding %s into factors", f->GetName());
  // use addClone, take the ownership of the objects
  // and release it when starting a new channel,
  // this can avoid potential memory leakage.
  return (m_arglist->addClone(*f) != NULL);
}


void Coefficient::SetCutoff (float in) {
    m_customCutoff.first=true; m_customCutoff.second=in;
}
