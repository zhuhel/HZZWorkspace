#include "Hzzws/Helper.h"
#include "Hzzws/Constants.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>

#include "TSystem.h"
#include <TString.h>
#include <TRegexp.h>
#include <TMath.h>
#include "RooStats/HistFactory/RooBSplineBases.h"
#include "RooStats/HistFactory/RooBSpline.h"

#include <algorithm>

namespace Helper{


bool readConfig(const char* filename, char delim,
        map<string, map<string, string> >& all_dic)
{
    log_info("reading: %s",filename);
    ifstream file(filename, ifstream::in);
    if(file.fail() || file.bad()){
        log_err("file: %s is bad.",filename);
        return false;
    }
    string line;
    int lineCount = 0;
    map<string, string> section_dic;
    string section_name;
    string tagName;
    string token;
    while ( getline(file, line) )
    { 
        boost::algorithm::trim(line);
        if (line.empty()) continue;
        if( line[0] == '[' ){
            if( lineCount < 1 ){
                section_name = string(line.begin()+1, line.end()-1);
            }else{
                all_dic[section_name] = section_dic;
                section_dic.clear();
                section_name = string(line.begin()+1, line.end()-1);
            }
        }else if( line[0] == '#' ){
            continue;
        }else{
            size_t delim_pos = line.find(delim);
            if (delim_pos != string::npos) {
                tagName = line.substr(0, delim_pos-1);
                token = line.substr(delim_pos+1, line.size());
                boost::algorithm::trim(tagName);
                boost::algorithm::trim(token);
                section_dic[tagName] = token;
            } else { // allow multiple lines for one tagName
                token = line;
                boost::algorithm::trim(token);
                section_dic[tagName] += " "+token;
                //cerr << line << " does not have delimeter '" << delim << "', ignored" << endl;
            }
        }
        lineCount ++ ;
    }
    all_dic[section_name] = section_dic;  //pick up the last section
    file.close();
    return true;
}

void tokenizeString(const string& str, char delim, vector<string>& tokens)
{
    tokens.clear();
    istringstream iss(str);
    string token;
    while ( getline(iss, token, delim) ){
        boost::algorithm::trim(token);
        tokens.push_back(token);
    }
}

void tokenizeString(const char* str, char delim, vector<string>& tokens)
{
    string tmp_str(str);
    tokenizeString(tmp_str, delim, tokens);
}

void tokenizeString(const char* str, char delim, vector<TString>& tokens){
  string tmp_str(str);
  strvec tokens2;
  tokenizeString(tmp_str, delim, tokens2);
  for (auto & s: tokens2) tokens.push_back(s.c_str());
}


//this function produces a vector list of files, and (optionally) maps them by wildcard
//
tstrvec fileList(const char* pattern, std::map<float, TString>* map)
{
  tstrvec files;
  TString pipe = gSystem->GetFromPipe(Form("ls -1 %s",pattern));
  pipe.ReplaceAll(" ",""); //remove empty space
  tokenizeString(pipe.Data(),'\n',files);
  for (auto& s : files) s = gSystem->GetFromPipe(Form("readlink -f %s",s.Data()));

  //attempt to map by floats in the files
  if (map){
    //check if the string contains more than one wildcard (if so, bail out!)

    for (auto& s : files) {
        Ssiz_t pos_slash = s.Last('/');
        TString file_name(s(pos_slash+1, s.Length()));
        TRegexp match_number("[0-9]+.[0-9]?[0-9]?");
        TString val(file_name(match_number));
        if (!val.IsFloat()) log_err("file name: %s and %s not a float.", file_name.Data(), val.Data());

        (*map)[val.Atof()] = s;
    }

    cout << "produced map from file list:" << endl;
    for (auto& m:*map) cout << m.first << " : " << m.second << endl;
  } else {
      std::cout<<"produced vector from file list:"<<std::endl;
      for (auto& f: files) cout << f << endl;
  }

  return files;
}


void readAcceptancePoly(std::vector<double>& params, const char* prod, const char* chan, const char* sys) { 

    if (!params.empty()){
    cerr <<"ERROR: "<< __func__ <<" given params vector must be empty!"<<endl;
    return;
  }

  const char* inPolyFile = Form("%s/polyNorm.txt", getInputPath().c_str());
  string line;
  string paramString;
  ifstream myfile(inPolyFile);
  if (myfile.is_open()){
    //loop until we find prod and chan
    while (getline(myfile, line) ){
      if (line.find(prod)!=string::npos && line.find(chan)!=string::npos){
        //loop until we find right sys
        while (getline(myfile, line) ){
          if (line.find(sys)!=string::npos){
            paramString = line; 
            break;
          }
        }
        break;
      }
    }
    myfile.close();
  }
  else {
    cerr<<"ERROR: file "<<inPolyFile<<" is not open!"<<endl;
    return;
  }

  if (paramString.empty()){
    cerr<<"ERROR: could not find in "<<inPolyFile<<": ["<<prod<<" "<<chan<<"] "<<sys<<endl;
    return;
  }

  vector<string> tokens;
  tokenizeString(paramString.c_str(),' ',tokens);
  for (int i(1);i<(int)tokens.size();++i)
    params.push_back(stod(tokens[i]));
}


RooRealVar* createNuisanceVar(const char* npName)
{
    string npVarName(Form("alpha_%s",npName));
    RooRealVar* npVar = new RooRealVar(npVarName.c_str(), npVarName.c_str(), 0.0, -5., 5.);
    return npVar;
}

RooRealVar* createGlobalVar(const char* npName)
{
    string npVarName(Form("nom_%s",npName));
    RooRealVar* npVar = new RooRealVar(npVarName.c_str(), npVarName.c_str(), 0.0, -5., 5.);
    npVar->setConstant();
    return npVar;
}

RooAbsPdf* createConstraint(const char* npName)
{
    // TODO: may add other constraint functions here
    auto* var = createNuisanceVar(npName);
    auto* mean = createGlobalVar(npName);
    auto* sigma = new RooRealVar("sigma", "sigma", 1.0,1.0,1.0);
    sigma->setConstant(true);
    string _pdfname(Form("alpha_%sConstraint", npName ));
    RooGaussian* gauss = new RooGaussian(_pdfname.c_str(), _pdfname.c_str(), *var, *mean, *sigma);
    return gauss;
}

void readNormTable(const char* file_name, 
        dbldic& all_norm_dic,
        double multiplier )
{
    cout << "reading normalization table" << endl;
    cout << "input: " << file_name << endl;
    ifstream file(file_name, ifstream::in);
    if (!file.is_open()) {
        cerr << "ERROR: cannot open " << file_name << endl;
        exit(1);
    }
    string line;
    int lineCount = 0;
    map<string, double> sample_dic;
    vector<string> category_names;
    string section_name;
    while ( !file.eof() && file.good() )
    {
        getline(file, line);
        line = TString(line).ReplaceAll("\\","").Data();
        if (lineCount++ == 0){
           tokenizeString(line, '&', category_names);
        } else if( line[0] == '#' ){
            cout << line <<" ignored" << endl;
            continue;
        }else{
            istringstream iss(line);
            string cat_name;
            iss >> cat_name;
            if (cat_name == "") continue;
            int total = (int) category_names.size();
            double yield;
            int index = 0;
            while ( index < total ){
                char ch;
                iss >> ch >> yield;
                sample_dic[category_names.at(index)] = yield * multiplier;
                index ++;
            }
            all_norm_dic[cat_name] = sample_dic;
            sample_dic.clear();
        }
    } 
    file.close();
    cout << "end of reading normalization table: " << file_name << endl;
}

void readScaleFile(const char* file_name, map<string, double>& all_dic)
{
   ifstream file(file_name, ifstream::in);
   string sample_name;
   double scale_value;
   while (file >> sample_name >> scale_value){
        all_dic[sample_name] = scale_value;
   }
  file.close(); 
}

TChain* loader(const string& inFile_name, const string& chain_name)
{
    TString fullname = (gSystem->IsAbsoluteFileName(inFile_name.c_str())) ?  inFile_name.c_str() :  Form("%s/%s",getInputPath().c_str(),inFile_name.c_str());

    log_info("loader wants to find input called %s",inFile_name.c_str());
    log_info("loader going to load : %s", fullname.Data());

    TChain* chain = new TChain(chain_name.c_str());
    if(fullname.Contains(".root")) {
        int res = chain->Add(fullname.Data());
        Long64_t nentries = chain->GetEntries();
        if (res == 0 || nentries == 0){
            log_err("%s is Missing or contains no event", fullname.Data());
            delete chain;
            return NULL;
        } else {
            cout << "total events: " << nentries << " in " << fullname << endl;
            return chain;
        }
    }

    fstream input(fullname.Data(), fstream::in);
    string file_name;
    int ncounter = 0;
    while (input >> file_name){
      if (file_name[0]=='#'){
        std::cout<<"skipping: "<<file_name<<std::endl;
        continue;
      }

      cout << "adding: " << file_name << endl;
      chain->Add(file_name.c_str());
      ncounter ++;
    }
    cout << "total events: " << chain->GetEntries() << " in " << ncounter << " files." << endl;
    input.close();
    return chain;
}

TChain* loader(const vector<string>& inFile_name, const string& chain_name)
{
    TChain* chain = new TChain(chain_name.c_str());

    for(size_t i = 0; i < inFile_name.size(); i++)
    {
        TString fullname(inFile_name[i]);
        if(fullname.Contains(".root")) {
            chain->Add(fullname);
        }
        else
        {
            cout<<"ERROR: Trying to add file "<<fullname<<endl;
            cout<<"Not supported"<<endl;
            exit(1);
        }
    }
    cout << "total events: " << chain->GetEntries() << endl;
    return chain;
}

bool IsGoodTH1(TH1* h1){ return IsSafeTH1(h1) && (h1->Integral() != 0); }
bool IsSafeTH1(TH1* h1){ return (h1 != NULL) && (!TMath::IsNaN(h1->Integral())); }
bool TH1FHasEmptyBin(TH1F* h) {
  for (int b(1);b<=h->GetNbinsX();++b){
    if (h->GetBinContent(b)==0) return true;
  }
  return false;
}


void printStopwatch(TStopwatch& timer)
{
  double kestRealTime = timer.RealTime();
  double kestCpuTime  = timer.CpuTime();
  int real_h = 0, real_m = 0, real_s =0;
  int cpu_h = 0, cpu_m = 0, cpu_s =0;
  real_h = (int) floor(kestRealTime/3600.) ;
  real_m = (int) floor((kestRealTime - real_h * 3600)/60.);
  real_s = kestRealTime - real_h*3600 - real_m*60;

  cpu_h = (int) floor(kestCpuTime/3600.) ;
  cpu_m = (int) floor((kestCpuTime - real_h * 3600)/60.);
  cpu_s = kestCpuTime - real_h*3600 - real_m*60;

  printf("RealTime=%dH%dM%ds, CPU=%dH%dM%ds\n",real_h,real_m,real_s,cpu_h,cpu_m,cpu_s);
}

const std::string& getInputPath(std::string i){
  static const std::string path(i);
  return path;
}
const std::string& addPoiName(std::string i){
  static std::string pois;
  if (i!=""){
    if (pois=="") pois=i;
    else pois=pois+","+i;
  }
  return pois;
}
RooWorkspace* getWorkspace( RooWorkspace* i){
    static  RooWorkspace* w=i;
    return w;
}

float getSysCutoff(std::string type, float setValue){
    static std::map<std::string, float> cutoff;
    if (cutoff.find(type)==cutoff.end()) cutoff[type]=setValue;
    return cutoff[type];
}


RooArgSet& getDisconnectedArgs(){
    static RooArgSet keepargs;
    return keepargs;
}

void getListOfNames(const string& cut, strvec& name_list, strmap& name_map) {
    TString cutStr_org(cut); // backup the original one
    TString cutStr(cut);
    // firstly ignore the string between '[' and ']'
    while( cutStr.Contains("[") ) {
      auto inds = cutStr.Index("[");
      auto inde = cutStr.Index("]");
      cutStr.Remove(inds, inde+1-inds);
    }
    // strvec operation_list {"=", ">", "<", "&", "|", "!", "(", ")"," ", "*","/"," ","+","-"};
    for(auto& op : constant::OPERATION_LIST) {
        cutStr.ReplaceAll(op, ",");
    }
    tokenizeString(cutStr.Data(), ',', name_list);
    name_list.erase(std::remove(name_list.begin(),name_list.end(),""), name_list.end());

    // remove numbers from the list
    for(unsigned i = 0; i < name_list.size(); i++)
    {
        TString currStr(name_list[i]);
        if(currStr.IsFloat()) name_list.erase(name_list.begin() + i);
    }

    // remove duplicated names. 
    // unique--> only remove duplication if two objects are nearby
    // therefore need to sort vector first, then call unique.. 
    // But I don't want to do that. 
    // Orders matter!
    strvec new_list;
    new_list.clear();
    for(auto& x: name_list){
        //if( std::binary_search(new_list.begin(), new_list.end(), x) ) continue;
        // remove all duplication
        if( std::find(new_list.begin(), new_list.end(), x) != new_list.end() ) continue;
        new_list.push_back(x);
    }
    name_list.assign(new_list.begin(), new_list.end());

    // order strings by length, in case some parameter names share common characters, eg., mu_ggF_off, mu5, mu, ...
    std::sort(new_list.begin(), new_list.end(), sort_str_len);
    // find non-poi names, i.e., the ones like kgg[1.0/0.9/1.1], or mu10[10]
    for(auto& x: new_list){
      if(cutStr_org.Contains(x)) {
        auto indx = cutStr_org.Index(x) + x.length();
        // check if a '[' follows
        bool found_range=false;
        while(indx!=kNPOS) {
          if(cutStr_org(indx)=='[') {
            auto ind_end = cutStr_org.Index("]", indx);
            if(ind_end!=kNPOS) {
              name_map[x]=cutStr_org(indx+1, ind_end-(indx+1));
              found_range=true;
              break; // found it
            }
            else {
              exit(1);
            }
          } else indx = cutStr_org.Index(x, indx+1); // find the next one
        } 
        if(!found_range) name_map[x]="";
      } 
    }

}

void readObservable(const string& str, vector<string>& obs_str, string& branch_name)
{
    string tmp_str(str);
    Helper::tokenizeString(tmp_str, ',', obs_str);
    if(obs_str.at(0).find(":") != string::npos){
      vector<string> tmp_vec;
      Helper::tokenizeString(obs_str.at(0), ':', tmp_vec);
      branch_name = tmp_vec.at(0); 
      obs_str.at(0) = tmp_vec.at(1);
      }
    else{ // uses leftmost part of branch name: varName_aaaa_bbbb
      branch_name = obs_str.at(0);
      vector<string> tmp;
      Helper::tokenizeString(obs_str.at(0), '_', tmp);
      obs_str.at(0) = tmp.at(0);
      } 
} 

std::string removeSpaces(const std::string& s){TString s2=s.c_str();s2.ReplaceAll(" ","");return s2.Data();}

std::string extractDecay(const std::string& s){
    const char* n[4]={"4mu","4e","2e2mu","2mu2e"}; 
    for (int i(0);i<4;++i) if (s.find(n[i])!=std::string::npos) return n[i];
    return "";
}


bool isMathSyntax(const string& str){
    return std::find(constant::MATH_SYNTAX_LIST.begin(), constant::MATH_SYNTAX_LIST.end(), str) != constant::MATH_SYNTAX_LIST.end();
}

int isBoolean(const string& str){
    if (!(str=="true" || str=="false" || 
                str=="1" || str=="0" || 
                str=="t" || str=="f"))
    {
        return -1;
    }
    return (int) (str=="true" || str=="1" || str=="t");
}

    bool sort_str_len(const std::string A, const std::string B){
       return A.size() > B.size();
    }

// END OF NAMESPACE
}
