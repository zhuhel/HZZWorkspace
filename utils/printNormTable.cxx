#include "TString.h"
#include "TFile.h"
#include "RooAbsReal.h"
#include "RooAbsArg.h"
#include "HZZWorkspace/Helper.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooAddition.h"
#include "RooRealVar.h"
#include <stdio.h>
#include <iostream>
#include <map>
#include <algorithm>

void extractMatch(TString in, std::vector<TString>& options, TString& out){
    for (auto& o : options){
        if (in.Contains(o)) out=o;
    }
}

void calculateNormErrors(TString varname, RooWorkspace* wsp, float& nom, float& err){

    //std::cout<<"working on "<<varname<<std::endl;

    RooAbsReal* norm = wsp->function(varname.Data());
    nom = norm->getVal();

    std::vector<float> effects;
    RooArgSet* allNuis = (RooArgSet*)wsp->allVars().selectByName("alpha_*");
    TIterator* iter(allNuis->createIterator());
    for (RooAbsArg* a = (RooAbsArg*)iter->Next(); a!=0; a=(RooAbsArg*)iter->Next()){
        RooRealVar* alpha = dynamic_cast<RooRealVar*>(a);

        alpha->setVal(1.);
        float up = norm->getVal();
        alpha->setVal(-1.);
        float down = norm->getVal();
        alpha->setVal(0.);

        up = fabs(up-nom)/nom;
        down = fabs(down-nom)/nom;
        effects.push_back(0.5*(up+down));
    }

    err=0.;
    for (auto& e: effects) err+= (e*e);
    err = nom*sqrt(err);


}

int main(int argc, char** argv){

    std::vector<TString> sources;
    std::vector<TString> categories;
    std::vector<TString> skip;
    TString fin = "combined.root";

    for (int a(1);a<argc;++a){
        TString s = argv[a];
        if (s=="--help") {
            std::cout<<"optional arguments (only the first is defaulted):"<<std::endl;
            std::cout<<"--workspace=combined.root (specify path to your workspace)"<<std::endl;
            std::cout<<"--sources=qqZZ,ggZZ,ggH (specify sources and order to be included in table)"<<std::endl;
            std::cout<<"--categories=4e,4mu,2mu2e (specify categories to be included in table)"<<std::endl;
            std::cout<<"--skip=VVV (specify sources to skip in table)"<<std::endl;
            std::cout<<"\nNote that names to not need to be explicit - they will be wildcard matched. If the categories/sources have overlapping (ambiguous) names, your table will be a disaster. If no names are provided, they will be determined automatically."<<std::endl;
            return -1;
        }
        if (s.Contains("--workspace=")){
            s.ReplaceAll("--workspace=","");
            fin = s;
        }
        if (s.Contains("--sources=")){
            s.ReplaceAll("--sources=","");
            Helper::tokenizeString(s.Data(),',',sources);
        }
        if (s.Contains("--categories=")){
            s.ReplaceAll("--categories=","");
            Helper::tokenizeString(s.Data(),',',categories);
        }
        if (s.Contains("--skip=")){
            s.ReplaceAll("--skip=","");
            Helper::tokenizeString(s.Data(),',',skip);
        }
    }

    //some customizable instructions
    //std::vector<TString> sources = {"qqZZ","ggZZ","VVV","ZJets","all","total expected"};
    //std::vector<TString> categories = {"4mu","4e","2mu2e","2e2mu","comb"};
    //std::vector<TString> skip = {"all"};
    //std::vector<TString> skip = {};

    std::cout<<"trying to print norm table from "<<fin<<std::endl;

    TFile* file = new TFile(fin.Data(),"READ");
    RooWorkspace* wsp = (RooWorkspace*)file->Get("combined");

    if (categories.empty()){//determine categories automatically
        std::cout<<"no categories provided - determining them on-the-fly.."<<std::endl;
        RooCategory* cat = wsp->cat("channelCat");
        if (!cat) {
            std::cout<<"trying to determine categories on the fly, but didn't find channelCat in your workspace.."<<std::endl;
            return -1;
        }
        for (int i(0);i<999;++i){
            if (!cat->isValidIndex(i)) continue;
            cat->setIndex(i,false);
            TString candidateName=cat->getLabel();
            candidateName.ReplaceAll("Cat","");
            categories.push_back(candidateName);
        }
    }

    if (sources.empty()){
        std::cout<<"no sources provided - determining them on-the-fly.."<<std::endl;

        std::set<TString> uniquesources;
        
        RooArgSet* allNorms = (RooArgSet*)wsp->allFunctions().selectByName("nTot*");
        TIterator* iter(allNorms->createIterator());
        for (RooAbsArg* a = (RooAbsArg*)iter->Next(); a!=0; a=(RooAbsArg*)iter->Next()){
            TString candidateName = a->GetName();
            candidateName.ReplaceAll("nTotATLAS_Bkg_","");
            candidateName.ReplaceAll("nTotATLAS_Signal_","");
            candidateName.ReplaceAll("nTotATLAS_","");
            for (auto& c : categories) candidateName.ReplaceAll(Form("_%s",c.Data()),"");
            uniquesources.insert(candidateName);
        }
        for (auto & s : uniquesources) sources.push_back(s);
    }

    std::cout<<"Going to make table out of sources:"<<std::endl;
    for (auto& s : sources) std::cout<<s<<std::endl;
    std::cout<<"And categories:"<<std::endl;
    for (auto& c : categories) std::cout<<c<<std::endl;
    std::cout<<"But skipping:"<<std::endl;
    for (auto& s : skip) std::cout<<s<<std::endl;
    std::cout<<std::endl;
    sources.push_back("total expected");
    categories.push_back("combined");

    std::map<TString, std::map<TString, std::pair<float,float> > > normResults; //eg. normResults["qqZZ"]["4mu"]= <10., 1.>

    std::map< TString, RooArgList > listCombProd;
    std::map<TString, RooArgList > listCombChan;

    RooArgSet* allNorms = (RooArgSet*)wsp->allFunctions().selectByName("nTot*");
    TIterator* iter(allNorms->createIterator());
    for (RooAbsArg* a = (RooAbsArg*)iter->Next(); a!=0; a=(RooAbsArg*)iter->Next()){
        RooAbsReal* n = dynamic_cast<RooAbsReal*>(a);

        TString thissource,thiscat;
        extractMatch(n->GetName(),sources,thissource);
        extractMatch(n->GetName(),categories,thiscat);

        //skip if requested
        if (std::find(skip.begin(),skip.end(),thissource)!=skip.end()) continue;
        //or if not part of sources
        if (thissource=="") continue;
        //or categories
        if (thiscat=="") continue;


        //calculate values from workspace
        float nom,err;
        calculateNormErrors(n->GetName(),wsp,nom,err);

        //store single channel results
        normResults[thissource][thiscat] = std::make_pair(nom,err);

        listCombProd[thissource].add(*wsp->function(n->GetName()));
        listCombChan[thiscat].add(*wsp->function(n->GetName()));
    }

    for (auto& s : listCombProd){
        RooAddition* add = new RooAddition(Form("combined_yield_%s",s.first.Data()),"",s.second);
        listCombChan["combined"].add(*add);
        wsp->import(*add,RooFit::RecycleConflictNodes(),RooFit::Silence());
        float nom,err;
        calculateNormErrors(Form("combined_yield_%s",s.first.Data()),wsp,nom,err);
        normResults[s.first]["combined"] = std::make_pair(nom,err);
    }

    for (auto& s : listCombChan){
        RooAddition* add = new RooAddition(Form("combined_yield_%s",s.first.Data()),"",s.second);
        wsp->import(*add,RooFit::RecycleConflictNodes(),RooFit::Silence());
        float nom,err;
        calculateNormErrors(Form("combined_yield_%s",s.first.Data()),wsp,nom,err);
        normResults["total expected"][s.first] = std::make_pair(nom,err);
    }

    //////// SCARY CODE THAT WORKS :o

    //Sort the results as requested
    std::vector<std::pair< TString, std::vector< std::pair< TString , std::pair<float,float > > > > > normResultsSorted;
    for (auto& namekey : sources){
        if (normResults.count(namekey)){
            std::vector< std::pair< TString , std::pair<float,float > > > resChan;
            //put the sorted ones in first
            for (auto& catkey: categories){
                if (normResults[namekey].count(catkey)){
                    resChan.push_back(std::make_pair( catkey, normResults[namekey][catkey] ) );
                }
            }
            //then anyone who isn't sorted
            for (auto& remainderB : normResults[namekey]){
                if (std::find(categories.begin(),categories.end(),remainderB.first)==categories.end())
                    resChan.push_back( std::make_pair(remainderB.first, remainderB.second) );
            }
            //then insert it into the top vector
            normResultsSorted.push_back(std::make_pair(namekey, resChan));
        }
    }
    for (auto& remainderA : normResults){
        if (std::find(sources.begin(),sources.end(),remainderA.first)==sources.end()) {
            std::vector< std::pair< TString , std::pair<float,float > > > resChan;
            //put the sorted ones in first
            for (auto& catkey: categories){
                if (remainderA.second.count(catkey)){
                    resChan.push_back(std::make_pair( remainderA.first , remainderA.second[catkey] ) );
                }
            }
            //then anyone who isn't sorted
            for (auto& remainderB : remainderA.second){
                if (std::find(categories.begin(),categories.end(),remainderB.first)==categories.end())
                    resChan.push_back( std::make_pair(remainderB.first, remainderB.second) );
            }
            //then insert it into the top vector
            normResultsSorted.push_back(std::make_pair(remainderA.first, resChan));
        }
    }
            
    ////// END SCARY CODE


    for (auto& a : normResultsSorted){
        for (auto& b : a.second){
            std::cout<<a.first<<" "<<b.first<<" "<<Form("%.2f +- %.2f",b.second.first,b.second.second)<<std::endl;
        }
    }

    std::cout<<"\n\n\n"<<std::endl;

    //Printing table
    std::cout<<"Category ";
    for (unsigned int s(0);s < normResultsSorted.size();++s) std::cout<<" & "<<normResultsSorted[s].first<<" ";
    std::cout<<" \\\\ "<<std::endl;
    for (unsigned int c(0);c < normResultsSorted[0].second.size();++c){
        std::cout<<normResultsSorted[0].second[c].first<<" ";
        for (unsigned int s(0);s < normResultsSorted.size();++s) {
            std::cout<<" &  $"<<Form("%.2f\\pm%.2f",normResultsSorted[s].second[c].second.first, normResultsSorted[s].second[c].second.second)<<"$ ";
        }
        std::cout<<" \\\\ "<<std::endl;
    }

}

