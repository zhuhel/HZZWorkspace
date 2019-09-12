/*
 * =====================================================================================
 *
 *       Filename:  Expected_qmu.cxx
 *
 *    Description:  Get NLL values for background-only toys
 *
 *        Version:  1.0
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *
 * =====================================================================================
   Requires 9 arguments
   - Input workspace ROOT file
   - Output ROOT file name for nll
   - muName: POI name
   - poi_input: value of POI
   - Mass: mass point
   - Combined: assumed name of workspace
   - obsData: assumed data name (could swap for eg. asimovData)
   - ntoys: number of toys to use in this test
   - seed: initial value for randomizer

Assumes all POI are associated with signal.
Tests background-only hypothesis, signal + background hypothesis. Requires generating events.
------------------------------------------------------------------------------
 */

#include <stdlib.h>
#include <string>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include "HZZWorkspace/RooStatsHelper.h"
#include "HZZWorkspace/Helper.h"
#include "TMath.h"

using namespace std;

int main(int argc, char** argv)
{
    if (argc < 6){
        cout <<argv[0]<< " input_ws.root out_nll.root muName poi_input mass combined obsData ntoys seed" << endl;
        exit(1);
    }
  string input_name(argv[1]);
  string out_name(argv[2]);
  string muName(argv[3]);
  string poi_input(argv[4]); //example: "6.5,7.5,8.5";
  float mass = atof(argv[5]); // 1200.0

  const char* mcName = "ModelConfig";
  string  wsName("combined");
  string  dataName("obsData");
  if (argc > 6){
      wsName = string(argv[6]);
  }
  if (argc > 7){
    dataName = string(argv[7]);
  }
  int ntoys = 100;
  if (argc > 8){
      ntoys = atoi(argv[8]);
  }
  int seed = 1;
  if (argc > 9){
      seed = atoi(argv[9]);
  }
  cout << "Input: " << input_name << endl;
  cout << "Onput: " << out_name << endl;
  cout << "POI Name:" << muName << endl;
  cout << "ws Name:" << wsName << endl;
  cout << "data Name:" << dataName << endl;
  cout << "seed:" << seed << endl;

  auto* file_in = TFile::Open(input_name.c_str(), "read");
  auto* workspace = (RooWorkspace*) file_in->Get(wsName.c_str());
  auto mc = (RooStats::ModelConfig*) workspace ->obj("ModelConfig");
  auto obs_data = dynamic_cast<RooDataSet*>(workspace->obj(dataName.c_str()));
  auto pois = const_cast<RooArgSet*>(mc->GetParametersOfInterest());
  RooRealVar* poi_val;
  TIter poi_itr(pois->createIterator());

  // set POIs to zero
  if(obs_data){
      // profile NPs to b-only fit
      poi_itr.Reset();
      while( (poi_val = (RooRealVar*)poi_itr())) {
          poi_val->setVal(0);
          poi_val->setConstant(true);
      }
      unique_ptr<RooNLLVar> nll(RooStatsHelper::createNLL(obs_data, mc));
      RooStatsHelper::minimize(nll.get());
  } else {
      std::cout<<"could not find dataset:"<<dataName<<". Skipping conditional step!"<<std::endl;
  }
  workspace->saveSnapshot("condGO",*mc->GetGlobalObservables());
  workspace->saveSnapshot("condNP",*mc->GetNuisanceParameters());
  workspace->saveSnapshot("condPOI",*mc->GetParametersOfInterest());
  pois->Print("v");

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

  for(int itoy = 0; itoy < ntoys; itoy++){
      int curr_seed = itoy + seed;
      poi_itr.Reset();
      while( (poi_val = (RooRealVar*)poi_itr())) {
          poi_val->setVal(0);
          poi_val->setConstant(false);
      }
      RooDataSet* toyData = dynamic_cast<RooDataSet*>(RooStatsHelper::generatePseudoData(workspace, muName.c_str(), curr_seed));
      for (int i = 0; i < (int)poi_list.size(); i++)
      {
          double poival = poi_list.at(i);
          res["mH"]= mass;
          res["mu"] = poival;
          std::cout<<"\n\n mu="<<poival<<std::endl;
          RooStatsHelper::fitData(workspace, mcName, toyData, muName.c_str(), poival, res);
          physics->Fill();
      }
      delete toyData;
  }

  fout->cd();
  physics->Write();
  fout->Close();
  file_in->Close();
}
