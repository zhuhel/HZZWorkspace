#include <vector>
#include <string>

#include <stdio.h>
#include <iostream>

#include "Hzzws/PlotHelp.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "TFile.h"
#include "TString.h"
#include "RooRealVar.h"

int main(){

  TFile* file = new TFile(Form("combined.root"),"READ");
  RooWorkspace* ws = (RooWorkspace*)file->Get("combined");
  //ws->Print();


  using namespace RooFit;


  RooArgSet nset(*ws->var("m4l"));

  std::vector<const char*> chan { "ggF_2mu2e_13TeV", "ggF_4mu_13TeV", "ggF_4e_13TeV", "VBF_incl_13TeV" };

  TCanvas* c = new TCanvas("c","c",700,700);

  for (int i(0);i<4;++i){
    c->Clear();
    RooPlot* frame = ws->var("m4l")->frame();

    ws->var("XS_ggF")->setVal(0.);
    ws->var("XS_VBF")->setVal(0.);

    ws->pdf(Form("modelunc_%s",chan[i]))->plotOn(frame, Normalization(ws->pdf(Form("modelunc_%s",chan[i]))->expectedEvents(nset)),LineColor(kRed));

    std::vector<const char*> prod{"ggF","VBF"};
    std::vector<EColor> cols{kMagenta,kCyan};

    for (int j(0);j<2;++j){
      ws->var(Form("XS_%s",prod[j]))->setVal(1.);

      const char* nick = (i==3 ? "cbga_sum" : "cbga");

      std::vector<float> mass{200., 300., 400., 500., 600., 700., 800., 900., 1000.};

      for (int m(0);m<9;++m){

        ws->var("mH")->setVal(mass[m]);
        ws->pdf(Form("ATLAS_Signal_%s_%s_%s",prod[j],chan[i],nick))->plotOn(frame, Normalization(ws->function(Form("nTotATLAS_Signal_%s_%s",prod[j],chan[i]))->getVal()), LineColor(cols[j]));
      }

    }

    frame->SetTitle(Form(";m_{4l};f(m_{4l}) [GeV^{-1}]"));
    frame->Draw();

    myText(0.4,0.8,kBlack,Form("cat: %s",chan[i]),0.04);
    myLineText(0.5,0.74,kRed,kSolid,3,"f_{B}",0.05);
    myLineText(0.5,0.68,kMagenta,kSolid,3,"f_{ggF}",0.05);
    myLineText(0.5,0.62,kCyan,kSolid,3,"f_{VBF}",0.05);
    c->Print(Form("plots/pdfvalid_%s.eps",chan[i]));

  }
}

