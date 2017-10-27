#include "Hzzws/AnalyticHMBkg.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "TFile.h"
#include "RooDataSet.h"

using namespace RooFit;

int main(){

  strmap args;
  args["global"]="N(10)";

  auto c = new Coefficient(args);

  Helper::getInputPath("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/Workspaces/Analytical/Prod_v04/20160512/input_v0/");
  auto p = new AnalyticHMBkg("ATLAS_Signal_ggF", "analyticHM_params_qqZZ.txt", "");
  p->addCoefficient(c);


  RooRealVar x("m4l_constrained","m4l_constrained",140,1200);
  RooRealVar event_type("event_type","event_type",0,3);
  RooRealVar weight("weight","weight",0,999);
  RooArgSet obs(x);

  TFile* file = new TFile("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc/Nominal/combined/mc15_qq2ZZ.root","READ");
  TTree* tree = (TTree*)file->Get("tree_incl_all");
  RooDataSet * ds4mu = new RooDataSet("data4mu","data4mu",tree,RooArgSet(event_type,x,weight),"m4l_constrained>140 && m4l_constrained<1200 && event_type==0","weight");
  RooDataSet * ds4e = new RooDataSet("data4e","data4e",tree,RooArgSet(event_type,x,weight),"m4l_constrained>140 && m4l_constrained<1200 && event_type==1","weight");
  RooDataSet * ds2mu2e = new RooDataSet("data2mu2e","data2mu2e",tree,RooArgSet(event_type,x,weight),"m4l_constrained>140 && m4l_constrained<1200 && (event_type==2 || event_type==3)","weight");

  ds4mu->Print();
  ds4e->Print();
  ds2mu2e->Print();
  

  RooPlot* frame4mu = x.frame(Title("4mu"));
  ds4mu->plotOn(frame4mu);
  p->setChannel(obs, "ggF_4mu_13TeV", true);
  p->getPDF()->plotOn(frame4mu,LineColor(kRed));

  RooPlot* frame4e = x.frame(Title("4e"));
  p->setChannel(obs, "ggF_4e_13TeV", true);
  ds4e->plotOn(frame4e);
  p->getPDF()->plotOn(frame4e,LineColor(kCyan));

  RooPlot* frame2e2mu = x.frame(Title("2e2mu"));
  ds2mu2e->plotOn(frame2e2mu);
  p->setChannel(obs, "ggF_2mu2e_13TeV", true);
  p->getPDF()->plotOn(frame2e2mu,LineColor(kGreen+1));

  TCanvas* can = new TCanvas("c","c",1200,800);
  can->Divide(3);
  //can->cd(1)->SetLogy(); frame4mu->Draw();
  //can->cd(2)->SetLogy(); frame4e->Draw();
  //can->cd(3)->SetLogy(); frame2e2mu->Draw();
  can->cd(1); frame4mu->Draw();
  can->cd(2); frame4e->Draw();
  can->cd(3); frame2e2mu->Draw();
  can->Print("plots/test_AnalyticHMBkg.png");
}
