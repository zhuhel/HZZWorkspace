#include "TSystem.h"
#include "TFile.h"
#include "RooProduct.h"
#include "HZZWorkspace/PlotHelp.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Rtypes.h"
#include <map>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TLatex.h"

#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

void myText(Double_t x,Double_t y,Color_t color,const char *text, float tsize, int font=42);
void ATLASLabel(Double_t x,Double_t y,const char* text,Color_t color);

int main(){

  using namespace RooFit;

  TFile file("combined.root","READ");
  RooWorkspace* wsp = (RooWorkspace*)file.Get("combined");

  std::vector<std::string> prod;
  //prod.push_back("ggF");
  //prod.push_back("VBF");
  prod.push_back("Signal");

  std::vector<std::string> categories;
  std::vector<std::string> catlabels;
  std::vector<int> color;
  /*
  categories.push_back("ggF_4mu_13TeV"); color.push_back(kViolet); catlabels.push_back("ggF-enriched 4#mu");
  categories.push_back("ggF_4e_13TeV"); color.push_back(kOrange); catlabels.push_back("ggF-enriched 4e");
  categories.push_back("ggF_2mu2e_13TeV"); color.push_back(kGreen+2); catlabels.push_back("ggF-enriched 2#mu 2e");
  categories.push_back("VBF_incl_13TeV"); color.push_back(kRed); catlabels.push_back("VBF-enriched");
  */
  categories.push_back("ggF_4mu_13TeV"); color.push_back(kViolet); catlabels.push_back("4#mu");
  categories.push_back("ggF_4e_13TeV"); color.push_back(kOrange); catlabels.push_back("4e");
  categories.push_back("ggF_2mu2e_13TeV"); color.push_back(kGreen+2); catlabels.push_back("2#mu 2e");
  categories.push_back("ggF_2e2mu_13TeV"); color.push_back(kCyan); catlabels.push_back("2e 2#mu");

  std::vector<std::string> sys;
  std::vector<std::string> sysnames;
  std::vector<int> linestyles;
  sys.push_back("alpha_ATLAS_EL_EFF_ISO_TOTAL_1NPCOR_PLUS_UNCOR"); linestyles.push_back(2); sysnames.push_back("Electron ID Efficiency");
  sys.push_back("alpha_ATLAS_EL_EFF_RECO_TOTAL_1NPCOR_PLUS_UNCOR"); linestyles.push_back(3); sysnames.push_back("Electron Reco Efficiency");
  sys.push_back("alpha_ATLAS_MUON_TTVA_SYS"); linestyles.push_back(8); sysnames.push_back("Muon Track-to-Vertex Association");
  sys.push_back("alpha_ATLAS_MUON_ISO_SYS"); linestyles.push_back(4); sysnames.push_back("Muon Isolation");


  TCanvas can("can","",600,600);
  can.SetLeftMargin(0.20);

  //
  // LOOP OVER PRODUCTION MODES
  //

  for (unsigned int p(0);p<prod.size();++p){
    std::cout<<prod[p]<<std::endl;

    //wsp->var(Form("XS_%s",prod[p].c_str()))->setVal(1./wsp->var("ATLAS_LUMI")->getVal());

    //
    // LOOP OVER CATEGORIES
    //
    TCanvas can("can","",600,600);

    for (unsigned int c(0);c<categories.size();++c){
        std::cout<<categories[c]<<std::endl;

        wsp->var(Form("BRRatio_%s",categories[c].c_str()))->setVal(1.);

        RooPlot* frame = wsp->var("mH")->frame(Range(120,130));
        
        //std::string funcname = Form("nTotATLAS_Signal_%s_%s", prod[p].c_str(), categories[c].c_str());
        std::string funcname = Form("ACpol_%s", categories[c].c_str());
        if (!wsp->function(funcname.c_str())) {
            std::cout<<"ERROR! No "<<funcname<<" found! :("<<std::endl;
            continue;
        }

        RooAbsReal* acpol = wsp->function(Form("ACpol_%s", categories[c].c_str()));
        RooAbsReal* sysfactor = wsp->function(Form("fiv_all_%s", categories[c].c_str()));

        RooProduct* acwsys = new RooProduct(Form("prod_%s",categories[c].c_str()), "", RooArgList(*acpol,*sysfactor));

        acpol->Print();
        sysfactor->Print();
        acwsys->Print();

        acwsys->plotOn(frame, LineColor(color[c]),Precision(0.0001));


        for (unsigned int s(0);s<sys.size();++s){ //begin loop over systematics
            if (!wsp->var(sys[s].c_str())) {
                std::cout<<"failed to find systematic: "<<sys[s]<<std::endl;
                continue;
            }

            if (categories[c]=="ggF_4mu_13TeV" && sys[s].find("_EL_")!=std::string::npos) continue;
            if (categories[c]=="ggF_4e_TeV" && sys[s].find("_MU_")!=std::string::npos) continue;

            wsp->var(sys[s].c_str())->setVal(1.0);
            acwsys->plotOn(frame, LineStyle(linestyles[s]), LineColor(kBlack), LineWidth(1), Precision(0.0001));
            wsp->var(sys[s].c_str())->setVal(-1.0);
            acwsys->plotOn(frame, LineStyle(linestyles[s]), LineColor(kBlack), LineWidth(1), Precision(0.0001));
            wsp->var(sys[s].c_str())->setVal(0.0);

        } //end loop over systematics

        frame->GetXaxis()->SetTitle("m_{H}");
        frame->GetXaxis()->SetLabelSize(0.035);
        frame->GetYaxis()->SetTitle(Form("Acceptance"));
        frame->GetYaxis()->SetTitleOffset(1.4);
        frame->GetYaxis()->SetLabelSize(0.035);
        //frame->SetMaximum(0.45);
        wsp->var("mH")->setVal(120.);
        frame->SetMinimum(0.8*acwsys->getVal());
        frame->SetTitle("");
        frame->Draw();

        //float x=0.18;
        //float y=0.80;
        //if (categories[c]=="ggF_2mu2e"){
        //  x=0.55;
        //  y=0.40;
        //}
        float x = 0.20;
        float y = 0.85;

        //Decorate the plot
        for (unsigned int s(0),o(0);s<sys.size();++s){ //begin loop over systematics
            if (!wsp->var(sys[s].c_str())) {
                std::cout<<"failed to find systematic: "<<sys[s]<<std::endl;
                continue;
            }
            if (categories[c]=="ggF_4mu_13TeV" && sys[s].find("_EL_")!=std::string::npos) continue;
            if (categories[c]=="ggF_4e_13TeV" && sys[s].find("_MUON_")!=std::string::npos) continue;
            myLineText(x,y-0.05*(o+1),kBlack,linestyles[s],1,Form("#pm1#sigma %s",sysnames[s].c_str()),0.03);
            ++o;
        }
        myLineText(x,y,color[c],0,1,Form("%s #rightarrow %s",prod[p].c_str(),catlabels[c].c_str()),0.03);

        can.Print(Form("plots/normsys_%s_%s.eps",prod[p].c_str(),categories[c].c_str()));
    }//end loop over categories

  } //end loop over prods

}

