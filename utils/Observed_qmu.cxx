/*
 * =====================================================================================
 *
 *       Filename:  Observed_qmu.cxx
 *
 *    Description:  Get NLL value for observed data 
 *
 *        Version:  1.0
 *
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <string>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include "Hzzws/RooStatsHelper.h"
#include "Hzzws/Helper.h"
#include "TMath.h"


using namespace std;
int main(int argc, char** argv)
{
    if (argc < 6){
        cout <<argv[0]<< " input_ws.root out_nll.root muName poi_input mass combined obsData otherPOIZeroAndConstant" << endl;
        exit(1);
    }
  string input_name(argv[1]);
  string out_name(argv[2]);
  string muName(argv[3]);
  string poi_input(argv[4]); //example: "6.5,7.5,8.5";
  float mass = atof(argv[5]); // 1200.0

  cout << "Input: " << input_name << endl;
  cout << "Onput: " << out_name << endl;
  const char* mcName = "ModelConfig";
  string  wsName("combined");
  string  dataName("obsData");
  if (argc > 6){
      wsName = string(argv[6]);
  }
  if (argc > 7){
    dataName = string(argv[7]);
  }
  bool is_other_poi_const = false;
  if (argc > 8){
      is_other_poi_const = (bool) atoi(argv[8]);
  }
  cout << "Input: " << input_name << endl;
  cout << "Onput: " << out_name << endl;
  cout << "POI Name:" << muName << endl;
  cout << "ws Name:" << wsName << endl;
  cout << "data Name:" << dataName << endl;
  cout << "Other POI zero and constant:" << is_other_poi_const << endl;

  auto* file_in = TFile::Open(input_name.c_str(), "read");
  auto* workspace = (RooWorkspace*) file_in->Get(wsName.c_str());
  auto mc = (RooStats::ModelConfig*) workspace ->obj("ModelConfig");

  // set other POIs to zero
  auto pois = const_cast<RooArgSet*>(mc->GetParametersOfInterest());
  if(is_other_poi_const) {
      float other_poi_value = 0;
      RooStatsHelper::setOtherPOIs(pois, muName, other_poi_value, is_other_poi_const);
  }else{
      RooStatsHelper::setOtherPOIs(pois, muName, 1., is_other_poi_const);
  }
  pois->Print("v");

/**
  TIter itr_nuis(mc->GetNuisanceParameters()->createIterator());
  RooRealVar* tmp_obj;
  while ( (tmp_obj = (RooRealVar*) itr_nuis()) ){
      tmp_obj->setConstant(true);
  }
  ***/

  auto* fout = TFile::Open(out_name.c_str(),"recreate");
  TTree* physics = new TTree("physics","physics");
  map<string, double> res;
  res["nll_SB_free"] = -1;
  res["poi_SB_free"] = -1;
  res["status_SB_free"] = -1;
  res["nll_SB_fixed"] = -1;
  res["poi_SB_fixed"] = -1;
  res["status_SB_fixed"] = -1;
  res["nll_B_fixed"] = -1;
  res["poi_B_fixed"] = -1;
  res["status_B_fixed"] = -1;
  res["mH"] = 0;
  res["mu"] = 0;
  for( auto& dic : res){
    physics->Branch(dic.first.c_str(), 
        &(dic.second), Form("%s/D",dic.first.c_str()));
  }

  RooStatsHelper::setDefaultMinimize();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

  vector<double> poi_list;
  vector<string> items;
  Helper::tokenizeString(poi_input.c_str(), ',', items);
  cout <<"poi input: " << poi_input << endl;
  for(auto& str_in : items){
      cout <<"adding: " << str_in << endl;
      poi_list.push_back(atof(str_in.c_str()));
  }

  if(workspace->var("mH")){
      workspace->var("mH")->setVal(mass);
      workspace->var("mH")->setConstant();
  }

  for (int i = 0; i < (int)poi_list.size(); i++)
  {
      double poival = poi_list.at(i);
      res["mH"]= mass;
      res["mu"] = poival;
      std::cout<<"\n\n mu="<<poival<<std::endl;
      RooStatsHelper::fitData(workspace, mcName, dataName.c_str(), muName.c_str(), poival, res);
      physics->Fill();
  }

  fout->cd();
  physics->Write();
  fout->Close();
  file_in->Close();
}
