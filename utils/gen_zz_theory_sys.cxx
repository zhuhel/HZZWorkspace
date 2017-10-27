// produce shape histograms for QCD/PDF uncertainties of ggZZ/qqZZ
//
#include <stdlib.h>
#include <utility>
#include <string>
#include <vector>
#include <iostream>
#include <map>


#include "TFile.h"
#include "TH1F.h"
#include "TString.h"

#include "Hzzws/Helper.h"
using namespace std;

map<string, double> pdf_qqZZ(double m4l){
    double var = 0.0035 * sqrt(m4l - 30);
    double up =  var;
    double down = - var;
    map<string, double> out;
    out["up"] = up;
    out["down"] = down;
    return out;
}
map<string, double> pdf_ggZZ(double m4l){
    double var = 0.0066 * sqrt(m4l - 10);
    double up =  var;
    double down = - var;
    map<string, double> out;
    out["up"] = up;
    out["down"] = down;
    return out;
}
map<string, double> qcd_qqZZ(double m4l){
    double var = 0.01 * sqrt((m4l - 20)/13.);
    double up =   var;
    double down =  - var;
    map<string, double> out;
    out["up"] = up;
    out["down"] = down;
    return out;
}
map<string, double> qcd_ggZZ(double m4l){
    double var = 0.10 * sqrt((m4l + 40)/ 40.);
    double up = .04 + var;
    double down = .04 - var;
    map<string, double> out;
    out["up"] = up;
    out["down"] = down;
    return out;
}
// QCDscale_ggVV, pdf_VV_highMass
// QCDscale_VV,   pdf_VV_highMass

int main(int argc, char* argv[]){
    map<string, map<string, string> > all_dic;
    if( argc < 2){
        cout << "please provide a config file" << endl;
        return 1;
    }
    Helper::readConfig(argv[1], '=', all_dic);
    Helper::printDic<string>(all_dic);
    vector<string> m4l_range;
    Helper::tokenizeString(all_dic["main"]["m4l"], ',', m4l_range);
    vector<string> categories;
    Helper::tokenizeString(all_dic["main"]["categories"], ',', categories);
    int n_bins = atoi(m4l_range.at(0).c_str());
    double min_x = atof(m4l_range.at(1).c_str());
    double max_x = atof(m4l_range.at(2).c_str());
    cout << "m4l range: [" << min_x << "," << max_x << "] with " 
        << n_bins << " bins " << endl;

    double step = (max_x - min_x) / n_bins;
    for (const auto& key : all_dic){
        string sampleName = key.first;
        if (sampleName == "main") continue;
        auto& sample_dic = key.second;
        TFile* file = TFile::Open(Form("%s_Shape_th.root",sampleName.c_str()), "recreate");
        TH1F* h_pdf_up = new TH1F("pdf_up", "pdf", n_bins, min_x, max_x);
        TH1F* h_pdf_down = new TH1F("pdf_down", "pdf", n_bins, min_x, max_x);
        TH1F* h_qcd_up = new TH1F("qcd_up", "pdf", n_bins, min_x, max_x);
        TH1F* h_qcd_down = new TH1F("qcd_down", "pdf", n_bins, min_x, max_x);
        for(int i = 0; i < n_bins; i++){
            double m4l = min_x + i * step;
            map<string, double> pdf_var;
            map<string, double> qcd_var;
            if (sampleName == "ggZZ"){
                pdf_var = pdf_ggZZ(m4l);
                qcd_var = qcd_ggZZ(m4l);
            } else if (sampleName == "qqZZ") {
                pdf_var = pdf_qqZZ(m4l);
                qcd_var = qcd_qqZZ(m4l);
            }
            h_pdf_down->Fill(m4l, 1+pdf_var["down"]);
            h_pdf_up->Fill(m4l, 1+pdf_var["up"]);
            h_qcd_down->Fill(m4l, 1+qcd_var["down"]);
            h_qcd_up->Fill(m4l, 1+qcd_var["up"]);
        }
        string pdf_name = sample_dic.at("pdfName");
        string qcd_name = sample_dic.at("qcdName");
        for (auto& category : categories){
            TH1F* h_pdf_up_copy = (TH1F*) h_pdf_up->Clone(Form("m4l-%s-%s-up", pdf_name.c_str(), category.c_str()));
            TH1F* h_pdf_down_copy = (TH1F*) h_pdf_down->Clone(Form("m4l-%s-%s-down", pdf_name.c_str(), category.c_str()));
            TH1F* h_qcd_up_copy = (TH1F*) h_qcd_up->Clone(Form("m4l-%s-%s-up", qcd_name.c_str(), category.c_str()));
            TH1F* h_qcd_down_copy = (TH1F*) h_qcd_down->Clone(Form("m4l-%s-%s-down", qcd_name.c_str(), category.c_str()));
            h_pdf_up_copy->Write();
            h_pdf_down_copy->Write();
            h_qcd_up_copy->Write();
            h_qcd_down_copy->Write();
        }
        file->Close();
    }
}
