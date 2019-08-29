#pragma once 

/// This is an odd assortment of various methods, most 
/// of which don't actually load anything.
/// Ask the original authors why the have to go in here. 
/// I take NO responsibility for this... 

#include <TString.h>
#include <exception>
#include <TH1F.h>
#include <TVector2.h>
#include <TFile.h>
#include <TLegend.h>
#include <THStack.h>
#include <TF1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TChain.h>
#include <TCut.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TColor.h>

#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooTFnPdfBinding.h>
#include <RooBinning.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

using namespace std;
using namespace RooFit;

#include "HZZWorkspace/AtlasStyle.h"
#include "HZZWorkspace/AtlasUtils.h"


/// Put these methods, which used to be in a .c macro,
/// into a class as static methods, to make them accessible 
/// from within pyROOT. 

class  Loader{
public:
    typedef struct BranchInfo {
        string name_;
        int n_;
        float low_;
        float high_;
        BranchInfo(){
            name_ = "";
            n_ = 100;
            low_ = 0;
            high_ = 1000;
        }
    } BranchInfo;


    static TChain* loader(const char* inFile_name, const char* chain_name = "physics");
    static TH1F* draw_hist_from_chain(TChain* chain, const char* branch_name, 
            const TCut& cut, const char* hist_name, 
            int n_bins, float low_value, float high_value);

    static TH1F* draw_hist_from_chain(TChain* chain, const TCut& cut, 
            const char* hist_name, const BranchInfo& br);

    static TH1F* draw_hist_from_file(const char* file_name, const char* chain_name, 
            const TCut& cut, const char* hist_name, const BranchInfo& br);

    static TH1F* create_hist(const char* file_name, const char* chain_name,
            const char* branch_name, const TCut& cut, const char* hist_name,
            int n_bins, float low_value, float high_value, int color);

    static void save_hist(TH1F* h1, const char* out_file_name, 
            const char* hist_name = "met_all");

    static TPad* add_ratio_pad(TH1* h_signal, const TList& h_bkgs);

    static TPad* add_ratio_pad(TH1* h_signal, TH1* h_bkg);

    static void norm_hist(TH1* h1);

    static void CheckNull(TObject* obj);

    static void SetAtlasStyleCanvas(TCanvas* canvas, bool for_2d = false);

    static void SetAtlasStyleHist(TH1* h1);

    static void SetAtlasOpt();

    static void add_hist(TList* objarray, TH1F* h1);

    static TH1* SumHistsWithSysUncertainties(const TList& hists, 
            const vector<double>& systematics, bool only_sys) ;

    static void add_hist(TList* objarray, const string& hist_name);

    static double get_significance(double s, double b);

    static double get_significance_with_sysB(double s, double b, double sigmaB);

    static int get_roc(TH1F* sig, TH1F* bkg, bool reverse = false);

    static void print_after_cut(const string& name, TH1F* h1, int cutbin);

    static void get_ratio_and_error(float a, float b, float& f, float& error);

    static TH1F* generate_th(TH1F* h_s1, const char* hist_name, float low_, float hi_ = 0);

    static void print_correlation(TH2* h2);

    static double get_mT(double v1_x, double v1_y,  double v2_x, double v2_y);

    static void compare_two_hists(TH1* h1, TH1* h2,
            const char* x_title, const char* h1_name, const char* h2_name, bool is_log);

    static double GetMinOfHist(TH1* h1);

    static void GetMaxMin(const TList& hists, double& max, double& min);

    static void compare_hists(const TList& histograms, const vector<string>& tags, 
            const char* x_title, bool shape_only, bool is_log, bool color_me) ;

    static TH1* merge_hist_list(const TList& hists);

    static TH1* merge_stack(const THStack& stack);

    static TH1F * DrawOverflow(TH1F *h);

    static void clean_TList(TList& list);

    static void compare_hists_from_files(const vector<string>& file_names, 
            const vector<string>& tag_name);

    static void print_graph(TGraph* gr);

    static double sum_graph_entries(RooCurve* gr, double low, double hi, int nbins);

};

