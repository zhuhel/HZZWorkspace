/*
 */
#include <stdlib.h>
#include <vector>

#include <RooRealVar.h>
#include <RooRealSumPdf.h>
#include <RooParamKeysPdf.h>
#include <RooProduct.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <TCanvas.h>
#include <TChain.h>

#include "Hzzws/ParametrizedSample.h"
#include "Hzzws/SampleHist.h"
#include <RooStats/HistFactory/RooBSpline.h>
#include <RooStats/HistFactory/RooBSplineBases.h>

void generate_para_ws();
void plot_para_norm();

const char* ws_file_name = "combined_para_test.root";
const char* ch_name = "all";

using namespace std;
int main(int argc, char** argv){
    generate_para_ws();
    // plot_para_norm();
    return 1;
}

void generate_para_ws()
{
    auto* workspace = new RooWorkspace("combined");

    const char* path = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v01/mc/Nominal/";
    const char* tree_name = "tree_incl_all";
    const char* f1 = Form("%s/mc15_13TeV.341504.PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_ggH124_ZZ4lep_noTau.root", path);
    const char* f2 = Form("%s/mc15_13TeV.341505.PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_ggH125_ZZ4lep_noTau.root", path);
    const char* f3 = Form("%s/mc15_13TeV.341506.PowhegPythia8EvtGen_CT10_AZNLOCTEQ6L1_ggH126_ZZ4lep_noTau.root", path);
    RooRealVar* m4l = new RooRealVar("m4l_constrained", "m4l", 110., 140.); 
    RooRealVar* mH = new RooRealVar("mH", "mH", 125.0, 110., 140.); 
    RooRealVar* one = new RooRealVar("one", "one", 1.0);
    double rho = 1.0;
    cout << f1 << endl;
    cout << f2 << endl;
    cout << f3 << endl;
    TChain* ch1 = new TChain(tree_name, tree_name);
    ch1->Add(f1);
    cout << "total entries: " << ch1->GetEntries() << endl;
    RooDataSet* ds1 = new RooDataSet("ds1", "ds1", RooArgSet(*m4l), RooFit::Import(*ch1));
    ds1->Print();
    cout << "h1.0" << endl;
    auto* key1 = new RooParamKeysPdf("key1", "key1", *m4l, *mH, 124.0, *one, 
            *ds1, RooParamKeysPdf::NoMirror,  rho);
    cout << "h1" << endl;
    key1->Print();
    workspace->import(*key1);
    
    TChain* ch2 = new TChain(tree_name, tree_name);
    ch2->Add(f2);
    RooDataSet ds2("ds2", "ds1", RooArgSet(*m4l), RooFit::Import(*ch2));
    auto* key2 = new RooParamKeysPdf("key2", "key1", *m4l, *mH, 125.0, *one, 
            ds2, RooParamKeysPdf::NoMirror,  rho);
    key2->Print();
    workspace->import(*key2);

    cout << "h2" << endl;

    TChain* ch3 = new TChain(tree_name, tree_name);
    ch3->Add(f3);
    RooDataSet ds3("ds3", "ds1", RooArgSet(*m4l), RooFit::Import(*ch3));
    auto* key3 = new RooParamKeysPdf("key3", "key1", *m4l, *mH, 126.0, *one, 
            ds3, RooParamKeysPdf::NoMirror,  rho);
    key3->Print();
    workspace->import(*key3);

    cout << "h3" << endl;
    vector<double>* masses = new vector<double>();
    masses->push_back(124.0);
    masses->push_back(125.0);
    masses->push_back(126.0);
    RooArgList* cps = new RooArgList();
    cps->add(*key1);
    cps->add(*key2);
    cps->add(*key3);
    auto* bases = new RooStats::HistFactory::RooBSplineBases("bases", "bases", 3, *masses, *mH);
    auto* bs_pdf = new RooStats::HistFactory::RooBSpline("bspline", "bspline", *cps, *bases, RooArgSet());

    workspace->import(*ds1);
    workspace->import(ds2);
    workspace->import(ds3);
    workspace->import(*bs_pdf, RooFit::RecycleConflictNodes());
    workspace->writeToFile(ws_file_name);

    delete workspace;
    delete bs_pdf;
    delete bases;
    delete cps;
    delete ch1;
    delete ch2;
    delete ch3;
}

void plot_para_norm()
{
    auto* file_in = TFile::Open(ws_file_name, "read");
    auto* workspace = (RooWorkspace*) file_in->Get("combined");
    RooRealSumPdf* pdf = (RooRealSumPdf*) workspace->obj("bspline");
    auto* bases = (RooStats::HistFactory::RooBSplineBases*) workspace->obj("bases");
    const vector<double>& tvalues = bases->getTValues();
    for(auto& value : tvalues){
        cout <<"value: " << value << endl;
    }
    
    RooRealVar* mH = (RooRealVar*) workspace->var("mH");
    RooRealVar* m4l = (RooRealVar*) workspace->var("m4l_constrained");
    mH->Print();
    mH->setVal(125);

    auto* canvas = new TCanvas("c1", "c1", 600, 600);
    // m4l->setRange(140., 700.);
    auto* m4l_frame = m4l->frame();
    int color = 2;
    for(double ini_mass=124.0; ini_mass <= 126.0; ini_mass += 0.5)
    {
        mH->setVal(ini_mass);
        pdf->plotOn(m4l_frame, RooFit::LineStyle(7), 
                RooFit::LineColor(color++),
                RooFit::LineWidth(2));
    }
    m4l_frame->Draw();
    canvas ->SaveAs("test_pdf.eps");
    delete canvas;

    delete m4l_frame;

    file_in->Close();
}
