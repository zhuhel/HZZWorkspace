#include "Hzzws/CBGaussSum.h"
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

  Helper::getInputPath("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/2016_06_17/input00/");
  auto p = new CBGaussSum("ATLAS_Signal_VBF", "ggH_para_VBFCat.txt");
  p->addCoefficient(c);

  RooRealVar x("m4l_constrained_HM","m4l_constrained_HM",140,1300);
  RooRealVar dijet_invmass("dijet_invmass","dijet_invmass",0,5000);
  RooRealVar dijet_deltaeta("dijet_deltaeta","dijet_deltaeta",0,10);
  RooRealVar weight("weight","weight",0,999);
  RooArgSet obs(x);

  TFile* file = new TFile("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v02/mc/Nominal/mc15_13TeV.341278.PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_ggH600NW_ZZ4lep.root","READ");
  TTree* tree = (TTree*)file->Get("tree_incl_all");
  RooDataSet * ds = new RooDataSet("data","data",tree,RooArgSet(x,dijet_invmass,dijet_deltaeta,weight),"m4l_constrained_HM>140 && m4l_constrained_HM<1300 && dijet_invmass>400 && dijet_deltaeta>3.3","weight");

  ds->Print();
  

  RooPlot* frame = x.frame(Title("CBGausSum validation"),Range(500,700),Bins(40));
  ds->plotOn(frame);
  p->setChannel(obs, "dummyCat", true);
  p->getPDF()->plotOn(frame,LineColor(kRed));


  TCanvas* can = new TCanvas("c","c",1200,800);
  frame->Draw();
  can->Print("plots/test_CBGaussSum.png");
}
