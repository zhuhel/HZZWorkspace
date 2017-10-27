#include "Hzzws/RooHistParamPdf.h"

#include "TCanvas.h"
#include "RooPlot.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "TString.h"
#include "RooDataHist.h"
#include "RooRealVar.h"

#include <RooArgSet.h>

#include <iostream>
#include <stdlib.h>

using namespace RooFit;
using namespace std;

int main(int argc, char** argv){
    if ((argc > 1 && string(argv[1]) == "help") || argc < 2){
        cout << argv[0] << " kappa" << endl;
        return 0;
    }
    Double_t kappa_val = atof(argv[1]);
    TFile* f1 = TFile::Open("/afs/cern.ch/user/b/bistapf/public/Hllvv_Interference/ggH1200_wH15.root");
    TH1F* h_sig = (TH1F*) f1->Get("mT_mm_signal");
    TH1F* h_hb = (TH1F*) f1->Get("mT_mm_HB");
    TH1F* h_hH = (TH1F*) f1->Get("mT_mm_hH");

    // Add background
    TFile* f_bkg = TFile::Open("/afs/cern.ch/user/b/bistapf/public/Hllvv_Interference/ggZZ.root");
    TH1F* h_bkg = (TH1F*) f_bkg->Get("mT_mm");
    h_bkg->Scale(10.6/h_bkg->Integral());
    cout << "ggZZ bkg yileds in mm: " << h_bkg->Integral() << endl;


    RooRealVar kappa("kappa", "kappa", 0.10, 0., 1.);
    kappa.setVal(kappa_val);

    // Use the same range and binning as the histogram..
    RooRealVar mT("mT", "mT", 300., 0., 2000.);
    mT.setBins(40);

    RooHistParamPdf* pdf = new RooHistParamPdf("pdf", "pdf", mT, kappa, *h_sig, *h_hb, *h_hH, *h_bkg);

    pdf->Print();

    RooPlot* frame = mT.frame(Title("mumu"));

    TH1F* h_all = (TH1F*) h_sig->Clone("h_all");
    Double_t kv = kappa.getVal();
    h_all->Scale(kv);
    h_hb->Scale(TMath::Sqrt(kv));
    h_hH->Scale(TMath::Sqrt(kv));

    h_all->Add(h_hb);
    h_all->Add(h_hH);
    h_all->Add(h_bkg);

    h_all->SetLineColor(4);
    h_sig->SetLineColor(8);
    h_sig->Scale(kv);

    Double_t yield = h_all->Integral();
    // h_all->Scale(yield);
    cout << "Total: " << yield << endl;

    pdf->plotOn(frame, LineColor(2), Normalization(yield, RooAbsReal::NumEvent));
    frame->addTH1(h_all, "EP");
    frame->addTH1(h_sig, "EP");
    
    // datahist->plotOn(frame);

    TCanvas* can = new TCanvas("c","c",600, 600);
    can->SetLogy();
    frame->Draw();
    can->SaveAs(Form("test_kappa_%.2f_Log.eps", kappa_val));

    delete frame;
    delete pdf;
    f1->Close();

    return 0;
}
