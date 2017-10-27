#include "TString.h"
#include "TColor.h"

void plot_mH_morphing(){

    TString channels[4] = {"ggF_4mu_13TeV","ggF_4e_13TeV","ggF_2mu2e_13TeV","ggF_2e2mu_13TeV"};
    int cols[] = {kSpring, kSpring-9, kYellow-3, kOrange, kOrange+2, kRed, kRed-9, kPink+6, kMagenta+1, kViolet, kViolet+8, kBlue-1, kBlue-3, kBlue+1, kAzure, kAzure+10,kCyan+3,kCyan-7,kTeal,kTeal-6,kTeal+3,kTeal+9,kGreen+1,kGreen-6};

    TFile* file = new TFile("combined.root");
    RooWorkspace* w = (RooWorkspace*) file->Get("combined");

    w->var("mH")->setRange(0,200);

    TCanvas* can = new TCanvas("c","c",600,800);
    can->SetLeftMargin(0.13);

    for (int c(0);c<4;++c){

        int k=0;

        RooPlot* frame = w->var("m4l")->frame(RooFit::Title(Form("")));

        for (float m(120.);m<=130.;m+=0.5){
            w->var("mH")->setVal(m);
            w->pdf(Form("ATLAS_Signal_all_%s_Morph",channels[c].Data()))->plotOn(frame,RooFit::LineColor(cols[k++]));
        }
        frame->SetTitle(";m_{4l};P(m_{4l}|m_{H})");
        TLatex* txt1 = new TLatex(0.15,0.86,Form("%s",channels[c].Data()));
        txt1->SetNDC();
        txt1->SetTextFont(42);
        txt1->SetTextSize(0.035);
        frame->addObject(txt1);
        frame->GetYaxis()->SetTitleOffset(1.4);
        frame->GetYaxis()->SetTitleSize(0.04);
        frame->GetXaxis()->SetTitleSize(0.04);
        frame->Draw();
        can->Print(Form("plots/test_mH_morphing_%s.eps",channels[c].Data()));
    }
}
