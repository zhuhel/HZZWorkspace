#include <fstream>
#include <iostream>
#include "TRegexp.h"
#include "HZZWorkspace/Helper.h"
#include <sstream>
#include <algorithm>
#include <exception>
#include <sys/types.h>
#include "TStyle.h"
#include <boost/algorithm/string.hpp>
#include <dirent.h>
#include <errno.h>
#include "TSystem.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TText.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TVectorD.h"


int main(int argc, char* argv[]){

    if (argc==1){
        std::cout<<"provide workspace input path! \"checkNP /path/to/files/\""<<std::endl;
        return -1;
    }

    std::map<TString, float> wildMatch;
    const float normFactor=125000; //for MeV (scale, res sys)
    //const float normFactor=100; //for % (norm sys)
    wildMatch[""]=0.;
    wildMatch["MUON_SCALE"]=0.;
    wildMatch["EG_SCALE"]=0.;
    wildMatch["EG_RESOLUTION"]=0.;
    wildMatch["EL_EFF"]=0.;
    wildMatch["PRW_DATASF"]=0.;
    wildMatch["MUON_MS"]=0.;
    wildMatch["MUON_ID"]=0.;

    TString path = argv[1];
    float thresh=-1;
    bool doMini=false;
    if (argc>2) thresh = TString(argv[2]).Atof();

    std::cout<<"input directory: "<<path<<std::endl;
    std::cout<<"threshold set at "<<thresh<<std::endl;

    void* dirp = gSystem->OpenDirectory(path.Data());

    std::string str;
    const char* entry;
    while ((entry=(char*)gSystem->GetDirEntry(dirp))) {
        str = entry;
        size_t pos = str.find(".txt");

        if (pos==std::string::npos) continue;

        TString file = path+"/"+str.c_str();

        std::cout<<"Going to try to make NP plot out of "<<str<<std::endl;

        strdic p_dic;
        Helper::readConfig(file.Data(),'=',p_dic);


        if (p_dic.empty()) continue;

        for (auto& cat : p_dic){

            std::cout<<"Analyszing file "<<file<<" in category "<<cat.first<<std::endl;
            for (auto& w : wildMatch) w.second=0.;

            float max=0.001;

            std::map<float, std::pair<std::string, std::pair<float,float> > > sys_map;

            for (auto& sys : cat.second){
                std::string npname = sys.first;
                strvec vals;
                Helper::tokenizeString(sys.second,' ',vals);
                if (vals.size()!=2) continue;
                float down = TString(vals[0].c_str()).Atof();
                float up = TString(vals[1].c_str()).Atof();
                float favg = 0.5*fabs(up-1.0) + 0.5*fabs(down-1.0);
                for (auto & w : wildMatch) if (TString(npname.c_str()).Contains(w.first.Data())) w.second+= favg*favg;

                if ( fabs(up-1.0)>thresh ||
                        fabs(down-1.0)>thresh ||
                        fabs(up-down)>thresh) {
                    //std::cout<<"adding sys above threshold: "<<npname<<", down="<<down<<", up="<<up<<std::endl;
                    sys_map[fabs(up-1.0)+fabs(down-1.0)] = std::make_pair( npname, std::make_pair(down-1.0,up-1.0) );
                }
                else {
                    //std::cout<<"skipping sys below threshold: "<<npname<<std::endl;
                }

                while (100*fabs(up-1.0)>max || 100*fabs(down-1.0)>max) max*=2;
            }

            gStyle->SetHistMinimumZero();

            TH1F* hist_down = new TH1F(Form("%s_down",cat.first.c_str()),"", sys_map.size(),0,sys_map.size());
            TH1F* hist_up = new TH1F(Form("%s_up",cat.first.c_str()),"", sys_map.size(),0,sys_map.size());

            hist_down->SetFillColor(kAzure+1);
            hist_down->SetFillStyle(3145);
            hist_down->SetBarWidth(0.75);
            hist_down->SetStats(0);
            hist_down->SetMinimum(-1.*max);
            hist_down->SetMaximum(max);
            hist_down->SetBarOffset(0.125);

            hist_up->SetFillStyle(3154);
            hist_up->SetFillColor(kBlue+2);
            hist_up->SetBarWidth(0.75);
            hist_up->SetBarOffset(0.125);

            int i(0);
            for (auto& sys: sys_map){
                ++i;
                TString lab = sys.second.first.c_str();
                lab.ReplaceAll("ATLAS_","");
                hist_down->GetXaxis()->SetBinLabel(i,lab.Data());

                hist_down->SetBinContent(i  , 100*sys.second.second.first);
                hist_up->SetBinContent(i    , 100*sys.second.second.second);
            }

            TCanvas* cv = new TCanvas("can","can",600,doMini?400:800);
            gPad->SetRightMargin(0.03);
            gPad->SetLeftMargin(0.5);
            gPad->SetBottomMargin(0.10);
            hist_down->GetYaxis()->SetNdivisions(6,10,0);
            hist_down->Draw("hbar");
            hist_up->Draw("hbar same");
            hist_down->GetXaxis()->SetLabelSize(doMini?0.06:0.03);
            hist_down->GetYaxis()->SetTitle("Down/up [%]");
            hist_down->GetYaxis()->SetLabelSize((doMini?0.04:0.03));
            hist_down->GetYaxis()->SetTitleSize((doMini?0.045:0.03));
            hist_down->GetYaxis()->SetTitleOffset(1);
            hist_down->GetYaxis()->CenterTitle();
            TLegend *leg=new TLegend(.79,.905,.99,doMini?0.98:.95);
            leg->SetFillColor(0);
            leg->SetBorderSize(0);
            leg->SetTextFont(42);
            leg->SetTextSize(doMini?0.05:0.035);
            leg->AddEntry(hist_down, "Down", "F");
            leg->AddEntry(hist_up,   "Up", "F");
            leg->Draw(); 
            TLatex* lm = new TLatex();
            lm->SetTextSize(doMini?0.06:0.035);
            lm->SetTextFont(42);
            lm->DrawLatexNDC(.45,0.92, cat.first.c_str());
            lm->Draw();

            TString savename = str.c_str();
            savename.ReplaceAll(".txt","");
            cv->Print(Form("plots/NPcheck_%s_%s.eps",savename.Data(),cat.first.c_str()));
            std::cout<<"\n===================="<<std::endl;
            std::cout<<savename.Data()<<" in category "<< cat.first <<std::endl;
            for (auto& w : wildMatch)
                std::cout<<(w.first==""?"total":w.first)<<" & "<<Form("%.1f",sqrt(w.second)*normFactor)<<std::endl;
            std::cout<<"====================\n"<<std::endl;
            delete leg;
            delete hist_down;
            delete hist_up;
            delete cv;
        }
    }
}

