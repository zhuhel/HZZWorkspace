#include "TCanvas.h"
#include "THStack.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "HZZWorkspace/PlotHelp.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"

#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>

//Important flags to change
float lumi=36.1;
bool doNWA=true;
std::string workspace="/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/Prod_v12/20170505_NWA/combined_scaled_corrected_buggy_sherpa_qqZZ.root";
std::string localdata="/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v12/data/Nominal/data_13TeV.root";
// ------------------------------


// a bit of plotting options too:


const char* rangenick="";
int hbins(55);
int hbinsVBF(7); //must be a divisor of hbins!
float hmin(130);
float hmax(1230);


/*
const char* rangenick="_peak";
int hbins(14);
int hbinsVBF(7); //must be a divisor of hbins!
float hmin(160);
float hmax(300);
*/


TGraphErrors* HistToGraph(TH1F*);
using namespace std;

int main(){

  using namespace RooFit;

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TH1::SetDefaultSumw2(true);

  std::vector<std::string> types;
  types.push_back("data");
  types.push_back("ZJets");
  types.push_back("VVV");
  types.push_back("ZZ");

  std::map<std::string, int> colmap;
  colmap["data"]=kBlack;
  colmap["ZZ"]=TColor::GetColor("#ee0000");
  colmap["VVV"]=TColor::GetColor("#EDE810");
  colmap["ZJets"]=kViolet+1;

  std::map<std::string, int> linecolmap;
  linecolmap["data"]=kBlack;
  linecolmap["ZZ"]=TColor::GetColor("#A10000");
  linecolmap["VVV"]=TColor::GetColor("#CCC70E");
  linecolmap["ZJets"]=kViolet+5;

  std::vector<std::string> cats;
  cats.push_back("ggF_2mu2e_13TeV");
  cats.push_back("ggF_4e_13TeV");
  cats.push_back("ggF_4mu_13TeV");
  if (doNWA) cats.push_back("VBF_incl_13TeV");

  std::map<std::string, std::string> catlabels;
  catlabels["ggF_2mu2e_13TeV"]="ggF-enriched 2#mu2e";
  catlabels["ggF_4e_13TeV"]="ggF-enriched 4e";
  catlabels["ggF_4mu_13TeV"]="ggF-enriched 4#mu";
  catlabels["VBF_incl_13TeV"]="VBF-enriched";

  std::map<std::string, std::string> cutmap;
  cutmap["ggF_2mu2e_13TeV"]=Form("(pass_vtx4lCut==1&&%f<m4l_constrained_HM&&m4l_constrained_HM<%f&&(event_type==3||event_type==2) && !(dijet_invmass>400 && dijet_deltaeta>3.3))",hmin,hmax);
  cutmap["ggF_4mu_13TeV"]=Form("(pass_vtx4lCut==1&&%f<m4l_constrained_HM&&m4l_constrained_HM<%f&&(event_type==0) && !(dijet_invmass>400 && dijet_deltaeta>3.3))",hmin,hmax);
  cutmap["ggF_4e_13TeV"]=Form("(pass_vtx4lCut==1&&%f<m4l_constrained_HM&&m4l_constrained_HM<%f&&(event_type==1) && !(dijet_invmass>400 && dijet_deltaeta>3.3))",hmin,hmax);
  if (doNWA){
    cutmap["VBF_incl_13TeV"]=Form("(pass_vtx4lCut==1&&%f<m4l_constrained_HM&&m4l_constrained_HM<%f && (dijet_invmass>400 && dijet_deltaeta>3.3))",hmin,hmax);
  }
  else {
    cutmap["ggF_2mu2e_13TeV"]="(pass_vtx4lCut==1&&140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==3||event_type==2))";
    cutmap["ggF_4mu_13TeV"]="(pass_vtx4lCut==1&&140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==0))";
    cutmap["ggF_4e_13TeV"]="(pass_vtx4lCut==1&&140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==1))";
  }

  std::vector<float> uncertainty;
  uncertainty.push_back(0.08057); //2mu2e
  uncertainty.push_back(0.0787); //4mu
  uncertainty.push_back(0.0814); //4e
  if (doNWA) uncertainty.push_back(0.5490); //VBF
  uncertainty.push_back(0.08031); //inclusive
  uncertainty.push_back(0.08031); //inclusive-ggF

  std::vector<TString> shapeNPs;
  shapeNPs.push_back("alpha_ATLAS_Bkg_qqZZ_fit_2mu2e_b4");
  shapeNPs.push_back("alpha_ATLAS_Bkg_qqZZ_fit_4mu_b4");
  shapeNPs.push_back("alpha_ATLAS_Bkg_qqZZ_fit_4e_b4");
  shapeNPs.push_back("alpha_ATLAS_Bkg_qqZZ_fit_VBF_b4");

  std::map<std::string, std::map<std::string, TH1F*> > hists;

  //Loop to populate all histograms
  for (unsigned int i(0);i<types.size();++i){
    std::string sample = types[i];

    for (unsigned int c(0);c<cats.size();++c){

      const char* hname = Form("hist_%s_%s",sample.c_str(),cats[c].c_str());
      TH1F* hist = new TH1F(hname,"",hbins,hmin,hmax);

      if (sample=="data"){
        //read from minitree
        TChain* chain = new TChain("tree_incl_all","tree_incl_all");
        chain->Add(localdata.c_str());

        chain->Draw(Form("m4l_constrained_HM>>%s",hname),cutmap[cats[c]].c_str());
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.25);
      }
      else {
        //read from workspace
        TFile* file = new TFile(workspace.c_str(),"READ");
        RooWorkspace* w = (RooWorkspace*)file->Get("combined");

        std::vector<std::string> pdfname, normname; 
        if (sample=="VVV"||sample=="ZJets") {
          pdfname.push_back(Form("ATLAS_Bkg_%s_%s_m4l",sample.c_str(),cats[c].c_str()));
          normname.push_back(Form("nTotATLAS_Bkg_%s_%s",sample.c_str(),cats[c].c_str()));
        }
        if (sample=="ZZ"){
          pdfname.push_back(Form("ATLAS_Bkg_ggZZ_%s_ana_shape",cats[c].c_str()));
          pdfname.push_back(Form("ATLAS_Bkg_qqZZ_%s_ana_shape",cats[c].c_str()));
          normname.push_back(Form("nTotATLAS_Bkg_ggZZ_%s",cats[c].c_str()));
          normname.push_back(Form("nTotATLAS_Bkg_qqZZ_%s",cats[c].c_str()));
        }

        for (unsigned int i(0);i<pdfname.size();++i){
          if (!w->pdf(pdfname[i].c_str())){
              cout <<pdfname[i] << " does not exits" << endl;
              continue;
          }
          TH1* hi = w->pdf(pdfname[i].c_str())->createHistogram(Form("%s%d",hname,i),*w->var("m4l"),Binning(hbins,hmin,hmax));
          TH1* hifull = w->pdf(pdfname[i].c_str())->createHistogram("temp",*w->var("m4l"),Binning(1000,w->var("m4l")->getMin(),w->var("m4l")->getMax()));
          hi->Scale(w->function(normname[i].c_str())->getVal() / hi->Integral());
          hi->Scale(hifull->Integral(hifull->FindBin(hmin),hifull->FindBin(hmax))/hifull->Integral());
          //shape uncertainties
          if (1){
            for (auto& np: shapeNPs){
              if (!w->var(np.Data())) continue;

              w->var(np.Data())->setVal(1);
              TH1* up = w->pdf(pdfname[i].c_str())->createHistogram("tempup",*w->var("m4l"),Binning(1000,w->var("m4l")->getMin(),w->var("m4l")->getMax()));
              w->var(np.Data())->setVal(-1);
              TH1* down = w->pdf(pdfname[i].c_str())->createHistogram("tempdown",*w->var("m4l"),Binning(1000,w->var("m4l")->getMin(),w->var("m4l")->getMax()));
              w->var(np.Data())->setVal(0.);
              TH1* nom = w->pdf(pdfname[i].c_str())->createHistogram("tempnom",*w->var("m4l"),Binning(1000,w->var("m4l")->getMin(),w->var("m4l")->getMax()));

              for (int b(1);b<=hi->GetNbinsX();++b){
                int b2 = nom->FindBin(hi->GetBinCenter(b));
                float old = hi->GetBinError(b);
                float added = hi->GetBinContent(b) *0.5*(up->GetBinContent(b2)-down->GetBinContent(b2))/nom->GetBinContent(b2);
                added*=2.0; //multiply effect to be conservative
                if (added!=0 && added==added) hi->SetBinError(b , sqrt(old*old + added*added));
              }
              delete up;
              delete down;
              delete nom;
            }
          }
          delete hifull;
          hist->Add(hi);
        }
      }

      hist->SetFillColor(colmap[sample]);
      hist->SetLineColor(linecolmap[sample]);
      hists[sample][cats[c]]=hist;

    }
  }


  //Combine things and plot
  TCanvas* can = new TCanvas("c","c",800,800);
  can->SetMargin(0.09,0.04,0.09,0.03);
  for (unsigned int c(0);c<=cats.size()+1;++c){
    can->Clear();

    int rebin = (c<cats.size()&&cats[c]=="VBF_incl_13TeV" ? hbins/hbinsVBF : 1);

    THStack ths("ths",Form(";m_{4l} [GeV];Events/%d GeV",(int)(rebin*(hmax-hmin)/(hbins))));
    std::map<std::string, TH1F*> inclusivehists;
    TH1F summed("summed","summed",hbins/rebin,hmin,hmax);
    if (c>=cats.size()){
      for (unsigned int p2(0);p2<types.size();++p2) { 
        inclusivehists[types[p2]] = new TH1F(Form("incl_%s",types[p2].c_str()), "",hbins,hmin,hmax); 
        for (unsigned int c2(0);c2<cats.size();++c2) { if (c==cats.size()+1 && cats[c2]=="VBF_incl_13TeV") continue; inclusivehists[types[p2]]->Add(hists[types[p2]][cats[c2]]); }
        inclusivehists[types[p2]]->SetLineColor(kBlack);
        if (types[p2]=="data") {
          inclusivehists[types[p2]]->SetMarkerStyle(20); 
          inclusivehists[types[p2]]->SetMarkerSize(1.25);
          continue;
        }
        inclusivehists[types[p2]]->SetFillColor(colmap[types[p2]]);
        inclusivehists[types[p2]]->SetLineColor(linecolmap[types[p2]]);
        summed.Add(inclusivehists[types[p2]]);
        ths.Add(inclusivehists[types[p2]]);
      }
    }
    else {
      for (unsigned int p2(0);p2<types.size();++p2) { 
        if (types[p2]=="data") continue;
        TH1F* hnew = (TH1F*)hists[types[p2]][cats[c]]->Rebin(rebin,Form("%s_rebinned",hists[types[p2]][cats[c]]->GetName()));
        ths.Add(hnew); summed.Add(hnew);
      }
    }
    ths.Draw("hist");

    //add errors to histogram
    for (int b(1);b<=summed.GetNbinsX();++b) {
      float old = rebin*summed.GetBinError(b);//old errors might have gotten shrunk by rebinning - scale them back!
      float added = summed.GetBinContent(b)*uncertainty[c];
      summed.SetBinError(b, sqrt(old*old+added*added));
    }
    auto * graph = HistToGraph(&summed);
    graph->SetFillStyle(3254);
    graph->SetFillColor(kBlack);
    graph->Draw("2 same");

    TH1F* hdata;
    if (c>=cats.size()) hdata = inclusivehists["data"];
    else hdata = hists["data"][cats[c]];
    if (rebin!=1) hdata = (TH1F*)hdata->Rebin(rebin,Form("data_%s_rebinned",cats[c].c_str()));
    hdata->Draw("PEX0 same");

    float pmax=1600*lumi/(float)hbins/(float)rebin;
    float pmin=0.2/(float)rebin;
    ths.SetMinimum(pmin);
    ths.SetMaximum(pmax);
    ths.GetYaxis()->SetTitleOffset(1.2);
    can->Update();

    //Decorate plot
    ATLASLabel(0.12,0.90,"Internal",kBlack);
    if (c==cats.size()) myText(0.12,0.85,kBlack,"H #rightarrow ZZ* #rightarrow 4l, inclusive",0.035);
    else if (c==cats.size()+1) myText(0.12,0.85,kBlack,"H #rightarrow ZZ* #rightarrow 4l, ggF-enriched",0.035); //inclusive ggF-enriched
    else myText(0.12,0.85,kBlack,Form("H #rightarrow ZZ* #rightarrow 4l, %s",catlabels[cats[c]].c_str()),0.035);
    myText(0.12,0.80,kBlack,Form("13TeV, %.1f fb^{-1}",lumi),0.035);

    float legx(0.70),legy(0.92);
    myMarkerText(legx,legy,kBlack,20,1.25,"Data",0.035,true);
    myBoxText(legx,legy-0.05,0.05, colmap["ZZ"], 1001, linecolmap["ZZ"], "ZZ*", 0.035);
    myBoxText(legx,legy-0.10,0.05, colmap["VVV"], 1001, linecolmap["VVV"], "t#bar{t}+V, VVV", 0.035);
    myBoxText(legx,legy-0.15,0.05, colmap["ZJets"], 1001, linecolmap["ZJets"], "Z+Jets, t#bar{t}", 0.035);
    myBoxText(legx,legy-0.20,0.05, kBlack, 3254, kWhite, "Uncertainty", 0.035);

    std::cout<<"printint for cat="<<c<<std::endl;
    const char* channick;
    if (c==cats.size()) channick="inclusive";
    else if (c==cats.size()+1) channick="inclusive_ggFonly";
    else channick = cats[c].c_str();
    can->SetLogy();
    can->Print(Form("plots/results%s_%s%s.eps",rangenick,channick,(doNWA?"":"_LWA")));
    can->SaveAs(Form("plots/results%s_%s%s.root",rangenick,channick,(doNWA?"":"_LWA")));
    can->SetLogy(false);
    ths.SetMinimum(0);
    ths.SetMaximum(1.45*std::max(hdata->GetBinContent(hdata->GetMaximumBin()), summed.GetBinContent(summed.GetMaximumBin())));
    can->Print(Form("plots/results%s_%s%s_lin.eps",rangenick,channick,(doNWA?"":"_LWA")));
    can->SaveAs(Form("plots/results%s_%s%s_lin.root",rangenick,channick,(doNWA?"":"_LWA")));
  }
}


TGraphErrors* HistToGraph(TH1F* hist){
  std::vector<float> x,y,xe,ye;
  for (int i(1);i<=hist->GetNbinsX();++i){
    x.push_back(hist->GetBinCenter(i));
    y.push_back(hist->GetBinContent(i));
    xe.push_back(0.5*hist->GetBinWidth(i));
    ye.push_back(hist->GetBinError(i));
  }
  std::cout<<"making a tgraph out of "<<std::endl;
  for (unsigned int i(0);i<x.size();++i){
    std::cout<<x[i]<<"\t"<<y[i]<<"\t"<<ye[i]<<std::endl;
  }
  return new TGraphErrors(x.size(), &(x[0]), &(y[0]), &(xe[0]), &(ye[0]));
}
