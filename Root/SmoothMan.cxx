// 
//    Description:  
// 
#include "Hzzws/SmoothMan.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <exception>

#include "RooRealVar.h"
#include "TVectorD.h"

#include "Hzzws/Helper.h"

using namespace std;
SmoothMan::SmoothMan(const char *configFile) {
    readConfig(configFile);
}

SmoothMan::~SmoothMan() {

}

void SmoothMan::readConfig(const char *configFile) {
    Helper::readConfig(configFile, '=', m_dic);
    Helper::printDic(m_dic);
}

void SmoothMan::process() {
  string filedir = m_dic["main"]["filedir"];
  string outputname = m_dic["main"]["outputname"];
  string outdir = "./";
  try{
    outdir = m_dic["main"]["outdir"];
  }catch (out_of_range& oor){

  }

  bool makeValidPlots=false;
  try{
    std::string s = m_dic["main"]["validPlots"];
    if (s=="1" || s=="true" || s=="yes") makeValidPlots=true;
  }catch (out_of_range& oor){}

  Helper::getInputPath(filedir);

  std::cout<<"got filedir as: "<<filedir<<std::endl;
  std::cout<<"got outdir as:"<<outdir<<std::endl;
  std::cout<<"got output name as:"<<outputname<<std::endl;

  std::vector<string> samples;
  Helper::tokenizeString(m_dic["main"]["samples"], ',', samples);

  cout << "Smoothing samples." << endl;
  for (auto &s: samples) {
    if (!(m_dic["main"].count(s))) {
      cout << "No sample " << s << " defined, skipping." << endl;
      continue;
    }
    cout << "Sample: " << s << endl;

    string outname = Form("%s/%s_%s.root",outdir.c_str(),outputname.c_str(),s.c_str());
    std::cout<<"outfile will be called: "<<outname<<std::endl;

    vector<string> sampleArgs;
    Helper::tokenizeString(m_dic["main"][s], ',', sampleArgs);
    if (sampleArgs.size() < 2) {
      std::cout<<"sample "<<s<<" didn't provide args! must provide [rho], optionally [smoothing]"<<std::endl;
      continue;
    }

    std::vector<float> rho;
    std::string mirror="";

    //interpret args
    for (int i(1);i<sampleArgs.size();++i){
      TString s2 = sampleArgs[i].c_str();
      if (s2.IsFloat()) rho.push_back(s2.Atof());
      else if (s2.Contains("Mirror")) mirror = s2.Data();
      //add another arg interpretation here if needed
    }



    string infile = filedir + "/" + sampleArgs[0];

    Smoother *sm = new Smoother(outname.c_str(), rho, mirror);
    sm->setMakeValidPlots(makeValidPlots);
    sm->setNickname(s);
    processSmoother(sm, infile);

    delete sm;
  }
}

void SmoothMan::processSmoother(Smoother *sm, const string& infile_name) 
{
  std::cout<<__func__<<" on "<<infile_name<<std::endl;
    string oname, treename;
    RooArgSet treeobs;
    if (m_dic["main"].count("treename") && m_dic["main"].count("observables")) 
    {
        getObs("main", oname, treename, treeobs);
    }

    vector<string> categories;
    Helper::tokenizeString(m_dic["main"]["categories"], ',', categories);
    for (auto &c : categories) {
        if (!(m_dic.count(c))) {
            cout << "No category " << c << " defined, skipping." << endl;
            continue;
        }
        cout << "Category: " << c << endl;

        string cut = m_dic[c]["cut"];
        
        if ((m_dic[c].count("treename") && m_dic[c].count("observables")) || (m_dic["main"].count("treename") && m_dic[c].count("observables"))) {
            string onametemp, treenametemp;
            RooArgSet treeobstemp;
            getObs(c, onametemp, treenametemp, treeobstemp);
            sm->smooth(infile_name, onametemp + c, treenametemp, treeobstemp, cut);
        }
        else {
            sm->smooth(infile_name, oname + c, treename, treeobs, cut);
        }
    }
}

void SmoothMan::readObservable(const string& str, vector<string>& obs_str, string& branch_name)
{
    string tmp_str(str);
    Helper::tokenizeString(tmp_str, ',', obs_str);

    vector<string> nick_vec;
    Helper::tokenizeString(obs_str.at(0), ':', nick_vec);
    if (nick_vec.size()==1) branch_name = obs_str.at(0) = nick_vec.at(0);
    else if (nick_vec.size()==2) { branch_name = nick_vec.at(0); obs_str.at(0) = nick_vec.at(1); }
    else {
      std::cout<<"ERROR: invalid observable argument!! must be \"branchname,nBins,low,high\" or \"branchname:nickname,nBins,low,high\""<<std::endl;
      return;
    }
} 

void SmoothMan::getObs(string cat, string &oname, string &treename, RooArgSet &treeobs) 
{
    if(m_dic[cat].count("treename"))  treename = m_dic[cat]["treename"];
    else  treename = m_dic["main"]["treename"];

    vector<string> branch;  
    Helper::tokenizeString(m_dic[cat]["observables"], ';', branch);

    for (int i=0; i<(int)branch.size(); ++i){
        vector<string> tmp_obs;
        string tmp_branch;
        readObservable(branch[i], tmp_obs, tmp_branch);

        if (tmp_obs.size() != 4) {
            cout << "Wrong number of parameters for observable, skipping." << endl;
            return;
        }

        RooRealVar *v = new RooRealVar(tmp_branch.c_str(), tmp_branch.c_str(), atof(tmp_obs.at(2).c_str()), atof(tmp_obs.at(3).c_str()));
        v->setBins( atoi(tmp_obs.at(1).c_str()) ); 
        treeobs.add(*v);
        oname += tmp_obs.at(0) + "_";
    }
}

