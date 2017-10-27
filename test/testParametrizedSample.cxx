/*
 */
#include <stdlib.h>
#include <RooRealVar.h>
#include <RooRealSumPdf.h>
#include <RooProduct.h>
#include <RooPlot.h>
#include <RooWorkspace.h>
#include <TCanvas.h>

#include "Hzzws/ParametrizedSample.h"
#include "Hzzws/SampleHist.h"
#include "Hzzws/SampleKeys.h"
#include <RooStats/HistFactory/RooBSpline.h>

void generate_para_ws();
void generate_keys_ws();
void plot_para_norm();

const char* ws_file_name = "combined_para_test.root";
const char* ch_name = "ggF_4e_13TeV";
const char* obs_name = "m4l";
int main(int argc, char** argv){
    // generate_keys_ws();
    // generate_para_ws();
    plot_para_norm();
    return 1;
}

void generate_para_ws()
{
    /***
    const char* path = "/afs/cern.ch/user/x/xju/work/h4l/workspace/mc15_13TeV_v1/low_mass/para/";
    const char* norm_sys = "norm_ggF_125_Low.txt";
    const char* shape_sys = "ggF_125_Low_Shape.root";
    ***/
    const char* path = "/afs/cern.ch/user/x/xju/work/h4l/workspace/mc15_13TeV_v2/himass/";
    const char* norm_sys = "norm_ggF_125_Low.txt";
    const char* shape_sys = "ggF_125_Low_Shape.root";
    auto* sample_h124 = new SampleHist("ATLAS_Signal_ggH124",  "test_ggH_200.root", shape_sys);
    auto* sample_h125 = new SampleHist("ATLAS_Signal_ggH125",  "test_ggH_300.root", shape_sys);
    auto* sample_h126 = new SampleHist("ATLAS_Signal_ggH126",  "test_ggH_400.root", shape_sys);
    auto* sample_h127 = new SampleHist("ATLAS_Signal_ggH127",  "test_ggH_500.root", shape_sys);
    auto* sample_h128 = new SampleHist("ATLAS_Signal_ggH128",  "test_ggH_600.root", shape_sys);
    sample_h124->set_mass(200.0);
    sample_h125->set_mass(300.0);
    sample_h126->set_mass(400.0);
    sample_h127->set_mass(500.0);
    sample_h128->set_mass(600.0);
    auto* sample_ggF = new ParametrizedSample("ATLAS_Signal_ggH", "ggF");
    // sample_ggF->SetParaRange("nominal", 200, 600);
    sample_ggF->AddSample(sample_h124);
    sample_ggF->AddSample(sample_h125);
    sample_ggF->AddSample(sample_h126);
    sample_ggF->AddSample(sample_h127);
    sample_ggF->AddSample(sample_h128);

    auto* obs = new RooRealVar("m4l", "m4l", 140, 3000);
    bool with_sys = false;
    sample_ggF->setChannel(RooArgSet(*obs), ch_name, with_sys);
    //sample_ggF->addNormSys("ATLAS_MU_ESCALE");
    //sample_ggF->addShapeSys("ATLAS_MU_ESCALE");

    RooRealSumPdf* pdf = (RooRealSumPdf*) sample_ggF->getPDF();
    RooProduct* coeff = (RooProduct*) sample_ggF->getCoefficient();
    pdf->Print();
    coeff->Print();
    
    auto* workspace = new RooWorkspace("combined");
    workspace->import(*pdf, RooFit::RecycleConflictNodes());
    workspace->import(*coeff, RooFit::RecycleConflictNodes());
    workspace->writeToFile(ws_file_name);

    delete workspace;
    delete pdf;
    delete coeff;
    delete obs;
    delete sample_ggF;
}

void generate_keys_ws()
{
    auto* workspace = new RooWorkspace("combined");
    const char* path = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2015/MiniTrees/Prod_v01/mc/Nominal/";
    // const char* tree_name = "tree_incl_all";
    const char* f1 = Form("%s/mc15_13TeV.341946.Pythia8EvtGen_A14NNPDF23LO_ZH124_ZZ4l.root", path);
    const char* f2 = Form("%s/mc15_13TeV.341947.Pythia8EvtGen_A14NNPDF23LO_ZH125_ZZ4l.root", path);
    const char* f3 = Form("%s/mc15_13TeV.341948.Pythia8EvtGen_A14NNPDF23LO_ZH126_ZZ4l.root", path);

    const char* norm_sys = "norm_ggF_125_Low.txt";

    auto* sample_h124 = new SampleKeys("ATLAS_Signal_ggH124", 124, 110, 140,
            f1,"");
    auto* sample_h125 = new SampleKeys("ATLAS_Signal_ggH125", 125, 110, 140, 
            f2,"");
    auto* sample_h126 = new SampleKeys("ATLAS_Signal_ggH126", 126, 110, 140,
            f3,"");

    auto* sample_ggF = new ParametrizedSample("ATLAS_Signal_ggH");
    // sample_ggF->SetParaRange("nominal", 200, 600);
    sample_ggF->AddSample(sample_h124);
    sample_ggF->AddSample(sample_h125);
    sample_ggF->AddSample(sample_h126);

    auto* obs = new RooRealVar(obs_name, "m4l", 110, 140);
    bool with_sys = false;
    sample_ggF->setChannel(RooArgSet(*obs), ch_name, with_sys);

    RooRealSumPdf* pdf = (RooRealSumPdf*) sample_ggF->getPDF();
    RooProduct* coeff = (RooProduct*) sample_ggF->getCoefficient();
    pdf->Print();
    coeff->Print();

    // sample_h124->SaveDataSet(workspace);
    // sample_h125->SaveDataSet(workspace);
    // sample_h126->SaveDataSet(workspace);
    workspace->import(*pdf, RooFit::RecycleConflictNodes());
    // workspace->import(*coeff, RooFit::RecycleConflictNodes());
    workspace->writeToFile(ws_file_name);

    delete workspace;
    delete pdf;
    delete coeff;
    delete obs;
    delete sample_ggF;
}

void plot_para_norm()
{
    auto* file_in = TFile::Open(ws_file_name, "read");
    auto* workspace = (RooWorkspace*) file_in->Get("combined");
    const char* prod = "ZH";
    // RooRealSumPdf* pdf = (RooRealSumPdf*) workspace->obj(Form("ATLAS_Signal_ggH_%s_Para", ch_name));
    // auto* bases = (RooStats::HistFactory::RooBSplineBases*) workspace->obj("bases_ggF");
    RooRealSumPdf* pdf = (RooRealSumPdf*) workspace->obj(Form("ATLAS_Signal_%s_%s_Para", prod,ch_name));
    auto* bases = (RooStats::HistFactory::RooBSplineBases*) workspace->obj(Form("bases_%s", prod));
    const vector<double>& tvalues = bases->getTValues();
    for(auto& value : tvalues){
        cout <<"value: " << value << endl;
    }
    
    RooRealVar* mH = (RooRealVar*) workspace->var("mH");
    RooRealVar* m4l = (RooRealVar*) workspace->var(obs_name);
    mH->Print();
    mH->setRange(110,140.);
    mH->setVal(125);

    auto* canvas = new TCanvas("c1", "c1", 600, 600);
    m4l->setRange(110., 140.);
    auto* m4l_frame = m4l->frame();
    int color = 2;
    for(double ini_mass=124.0; ini_mass <= 126.0; ini_mass += 0.5)
    {
        // cout << "mass: " << ini_mass << endl;
        mH->setVal(ini_mass);
        // cout << "yield: " << coeff->getVal() << endl;
        pdf->plotOn(m4l_frame, RooFit::LineStyle(7), 
                RooFit::LineColor(color++),
                RooFit::LineWidth(2));
    }
    m4l_frame->Draw();
    canvas ->SaveAs("test_pdf_vary.eps");
    delete canvas;

/*
    auto* c2 = new TCanvas("c2", "c2", 600, 600);
    auto* mh_frame = mH->frame();
    coeff->plotOn(mh_frame, RooFit::LineStyle(7), RooFit::LineWidth(2));
    mh_frame->Draw();
    c2->SaveAs("test_base.eps");
    delete c2;
    delete mh_frame;

    auto* c3 = new TCanvas("c3", "c2", 600, 600);
    mH->setVal(125.5);
    color = 2;
    RooRealVar* muon_escale = (RooRealVar*) workspace->var("alpha_ATLAS_MU_ESCALE");
    cout << "initial np: " << muon_escale->getVal() << endl;
    for(double ini_np = -3; ini_np < 3; ini_np += 0.2){
        muon_escale->setVal(ini_np);
        pdf->plotOn(m4l_frame, RooFit::LineStyle(7), RooFit::LineWidth(2),
                RooFit::LineColor(color++));
    }
    m4l_frame->Draw();
    c3->SaveAs("test_np.eps");
    delete c3;
    */
    delete m4l_frame;

    file_in->Close();
}
