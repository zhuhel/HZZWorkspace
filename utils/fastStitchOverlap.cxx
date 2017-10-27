
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGaxis.h"

#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooNDKeysPdf.h"

using namespace std;
using namespace RooFit;

int main() {

  bool linearPlots=true;
  
  string outputNames[3]={"4mu","2mu2e","4e"};
  float rho_lo[3]={0.2, 0.2, 0.2}; // was 0.25 for all for Higgs approval
  float rho_hi[3]={0.6, 0.6, 0.6}; // was 0.35 for all for Higgs approval
  float rho_highest[3]={0.25, 0.25, 0.25}; 
  
  int m4lLow=140.;
  int m4lHigh=1200.;
  
  float binWidth=5.;
  int m4lBins = (int)((m4lHigh-m4lLow)/binWidth);

  float splitVal=275.; // was 225 for Higgs approval
  float overlap=75.; // was 30 for first Higgs approval
  float splitValHi=500.;
  float overlapHi=100.;
  int m4lBinsLow=(int)((splitVal+overlap-m4lLow)/binWidth);
  int m4lBinsHigh=(int)((splitValHi+overlapHi-splitVal+overlap)/binWidth);
  int m4lBinsHighest=(int)((m4lHigh-splitValHi+overlapHi)/binWidth);
  
  RooRealVar m4l("m4l","m4l",m4lBins,m4lLow,m4lHigh);
  RooRealVar m4lLowV("m4lLow","m4l",m4lBinsLow,m4lLow,splitVal+overlap);
  RooRealVar m4lHighV("m4lHigh","m4l",m4lBinsHigh,splitVal-overlap,splitValHi+overlapHi);
  RooRealVar m4lHighestV("m4lHighest","m4l",m4lBinsHighest,splitValHi-overlapHi,m4lHigh);

  // input tree
  TChain* tree=new TChain("tree_incl_all");
  //tree->Add("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v02/mc/Nominal/combined/mc15_qq2ZZ.root");
  tree->Add("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc_15b/Nominal/mc15_13TeV.343232.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4_m4l_500_13000.root");
  tree->Add("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v03/mc_15b/Nominal/mc15_13TeV.361603.PowhegPy8EG_CT10nloME_AZNLOCTEQ6L1_ZZllll_mll4.root");
  float tree_m4l, tree_weight;
  int tree_event_type;
  tree->SetBranchAddress("m4l_constrained", &tree_m4l);
  tree->SetBranchAddress("weight", &tree_weight);
  tree->SetBranchAddress("event_type", &tree_event_type); //  _4mu, _4e, _2mu2e, _2e2mu
  
  TH1F** hists=new TH1F*[3];
  TH1F** histsLow=new TH1F*[3];
  TH1F** histsHigh=new TH1F*[3];
  TH1F** histsHighest=new TH1F*[3];
  
  for (int ich=0; ich<3; ich++) {
    string histname  ="m4l_ggF_"+outputNames[ich]+"_raw";
    string histnameLo="m4l_ggF_"+outputNames[ich]+"_raw_lo";
    string histnameHi="m4l_ggF_"+outputNames[ich]+"_raw_hi";
    string histnameHighest="m4l_ggF_"+outputNames[ich]+"_raw_highest";
    hists    [ich]=new TH1F(histname  .c_str(),histname  .c_str(),m4lBins    , m4lLow          , m4lHigh);
    histsLow [ich]=new TH1F(histnameLo.c_str(),histnameLo.c_str(),m4lBinsLow , m4lLow          , splitVal+overlap);
    
    histsHigh[ich]=new TH1F(histnameHi.c_str(),histnameHi.c_str(),m4lBinsHigh, splitVal-overlap, splitValHi+overlapHi);
    
    histsHighest[ich]=new TH1F(histnameHighest.c_str(),histnameHighest.c_str(),m4lBinsHighest, splitValHi-overlapHi, m4lHigh);
  }
  
  for (int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (tree_m4l > m4lHigh || tree_m4l < m4lLow) continue;
    

    float nnlo_corr = 1.0;
    float weight_nnlo = tree_weight * nnlo_corr;

    int ich=-1;
    if (tree_event_type==0) ich=0;
    else if (tree_event_type==1) ich=2;
    else ich=1;

    hists[ich]->Fill(tree_m4l, weight_nnlo);
    if (tree_m4l<(splitVal+overlap)) 
      histsLow[ich]->Fill(tree_m4l, weight_nnlo);
    if (tree_m4l>(splitVal-overlap) && tree_m4l<(splitValHi+overlapHi))
      histsHigh[ich]->Fill(tree_m4l, weight_nnlo);
    if (tree_m4l>(splitValHi-overlapHi))
      histsHighest[ich]->Fill(tree_m4l, weight_nnlo);      
  }

  TH1F** rats=new TH1F*[3];

  TCanvas canv;
  canv.Divide(1,3);

  TFile* outFile=new TFile("qqZZbackground.root","RECREATE");
  RooWorkspace* myws=new RooWorkspace("qqZZbackground");
  
  for (int ich=0; ich<3; ich++) {

    string outputName = outputNames[ich];
    
    // Smoothing for high mass
    TH1F* hist=hists[ich];
    TH1F* histLo=histsLow[ich];
    TH1F* histHi=histsHigh[ich];
    TH1F* histHighest=histsHighest[ich];

    // low mass region
    m4lLowV.setRange(m4lLow, splitVal+overlap);
    m4lLowV.setBins(m4lBinsLow);

    RooArgSet obstempLo(m4lLowV);
    RooNDKeysPdf *keyspdflo = new RooNDKeysPdf((outputName + "_keyspdf_lo").c_str(), "keyspdf", obstempLo, *histLo, "m", rho_lo[ich]);

    TH1F* smoothedhistLo = (TH1F*) keyspdflo->createHistogram((outputName + "_TH1_smoothedLo").c_str(), m4lLowV, Binning(m4lBinsLow));
    smoothedhistLo->Scale(histLo->Integral() / smoothedhistLo->Integral());

    // high mass region
    m4lHighV.setRange(splitVal-overlap, splitValHi+overlapHi);
    m4lHighV.setBins(m4lBinsHigh);

    RooArgSet obstempHi(m4lHighV);
    RooNDKeysPdf *keyspdfhi = new RooNDKeysPdf((outputName + "_keyspdf_hi").c_str(), "keyspdf", obstempHi, *histHi, "ma", rho_hi[ich]);

    TH1F *smoothedhistHi = (TH1F*) keyspdfhi->createHistogram((outputName + "_TH1_smoothedHi").c_str(), m4lHighV, Binning(m4lBinsHigh));
    smoothedhistHi->Scale(histHi->Integral() / smoothedhistHi->Integral());

    // highest mass region
    m4lHighestV.setRange(splitValHi-overlapHi, m4lHigh);
    m4lHighestV.setBins(m4lBinsHighest);

    RooArgSet obstempHighest(m4lHighestV);
    RooNDKeysPdf *keyspdfhighest = new RooNDKeysPdf((outputName + "_keyspdf_highest").c_str(), "keyspdf", obstempHighest, *histHighest, "m", rho_highest[ich]);

    TH1F *smoothedhistHighest = (TH1F*) keyspdfhighest->createHistogram((outputName + "_TH1_smoothedHighest").c_str(), m4lHighestV, Binning(m4lBinsHighest));
    smoothedhistHighest->Scale(histHighest->Integral() / smoothedhistHighest->Integral());
    
    string smhistname="m4l_ggF_"+outputName+"_13TeV";
    TH1F *smoothedhist = new TH1F(smhistname.c_str(),smhistname.c_str(), m4lBins, m4lLow, m4lHigh);
    int overlapbins=(int)(overlap/binWidth+.001);
    int overlapbinsHi=(int)(overlapHi/binWidth+.001);
    
    // low
    for (int ibin=0; ibin<m4lBinsLow-overlapbins; ibin++) {
      smoothedhist->SetBinContent(ibin+1,smoothedhistLo->GetBinContent(ibin+1));
    }

    // high
    for (int ibin=0; ibin<m4lBinsHigh-overlapbins-overlapbinsHi; ibin++) {
      smoothedhist->SetBinContent(m4lBinsLow-overlapbins+ibin+1,smoothedhistHi->GetBinContent(overlapbins+ibin+1));
    }

    // highest
    for (int ibin=0; ibin<=m4lBinsHighest-overlapbinsHi; ibin++) {
      smoothedhist->SetBinContent(m4lBinsLow-2*overlapbins-overlapbinsHi+m4lBinsHigh+ibin+1,smoothedhistHighest->GetBinContent(overlapbinsHi+ibin+1));
    }

    smoothedhist->Scale(hist->Integral() / smoothedhist->Integral());

    smoothedhist->Write();
    hist->Write();
    
    m4l.setRange(m4lLow, m4lHigh);
    m4l.setBins(m4lBins);

    rats[ich]=(TH1F*)smoothedhist->Clone("rat");
    
    canv.cd(ich+1);

    if (!linearPlots) {
      TPad* pad1=new TPad("pad1","pad1",0,0.3,1,1.0);
      pad1->SetBottomMargin(0);
      pad1->SetGridx();
      pad1->Draw();
      pad1->cd()->SetLogy(1);
    }
    
    hist->SetStats(0);
    hist->SetLineColor(kBlack);
    hist->Draw();
    smoothedhist->SetLineColor(kRed);
    smoothedhist->Draw("same");

    hist->GetYaxis()->SetLabelSize(0.);
    TGaxis* axis = new TGaxis(-5,20,-5,220,20,220,510,"");
    axis->SetLabelFont(43);
    axis->SetLabelSize(15);
    axis->Draw();

    if (!linearPlots) {
      canv.cd(ich+1);
      TPad* pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
      pad2->SetTopMargin(0);
      pad2->SetBottomMargin(0.2);
      pad2->SetGridx();
      pad2->Draw();
      pad2->cd();
      
      TH1F* rat=rats[ich];
      rat->SetMinimum(0.);
      rat->SetMaximum(2.);
      rat->Sumw2();
      rat->SetStats(0);
      
      rat->Divide(hist);
      rat->SetMarkerStyle(7);
      rat->Draw("ep");
      
      hist->GetYaxis()->SetTitleSize(8);
      hist->GetYaxis()->SetTitleFont(43);
      hist->GetYaxis()->SetTitleOffset(1.55);
      
      rat->SetTitle("");
      
      rat->GetYaxis()->SetTitle("ratio smoothed hist/raw hist");
      rat->GetYaxis()->SetNdivisions(5);
      rat->GetYaxis()->SetTitleSize(8);
      rat->GetYaxis()->SetTitleFont(43);
      rat->GetYaxis()->SetTitleOffset(1.55);
      rat->GetYaxis()->SetLabelFont(43);
      rat->GetYaxis()->SetLabelSize(8);
      
      rat->GetXaxis()->SetTitleSize(20);
      rat->GetXaxis()->SetTitleFont(43); 
      rat->GetXaxis()->SetTitleOffset(4.);
      rat->GetXaxis()->SetLabelFont(43);
      rat->GetXaxis()->SetLabelSize(8);
    }
  }
  if (linearPlots)
    canv.Print("hists_linear.eps");
  else
    canv.Print("hists.eps");
    
  myws->Write();
  outFile->Close();
  
}
