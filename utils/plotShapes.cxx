/*
 * Plot the shape
 */

#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooSimultaneous.h>
#include <RooCategory.h>
#include <RooStats/ModelConfig.h>
#include <RooCurve.h>
#include "RooMinimizer.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TIterator.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TColor.h>


#include <stdlib.h>
#include <string>

#include "Hzzws/Helper.h"
#include "Hzzws/RooStatsHelper.h"

using namespace std;
int main(int argc, char** argv)
{
    if ((argc > 1 && string(argv[1]) == "help") ||
            argc < 6)
    {
        cout << argv[0] << " combined.root ws_name mu_name data_name mc_name fitmode with_data do_visual_error min,max strategy color var:value,var:value np1,np2" << endl;
        return 0;
    }
    string input_name(argv[1]);
    string wsName(argv[2]);
    string muName(argv[3]);
    string dataName(argv[4]);
    string mcName(argv[5]);

    int opt_id = 6;

    int fit_mode = -1;
    if (argc > opt_id) {
        fit_mode = atoi(argv[opt_id]);
    }
    opt_id ++;

    bool with_data = false;
    if (argc > opt_id) {
        with_data = (bool) atoi(argv[opt_id]);
    }
    opt_id ++;

    bool do_visual_error = true;
    if (argc > opt_id) {
        do_visual_error = (bool) atoi(argv[opt_id]);
    }
    opt_id ++;

    double max_obs = 2000, min_obs = 200;
    if (argc > opt_id) {
        string options(argv[opt_id]);
        vector<string> tokens;
        Helper::tokenizeString(options, ',', tokens);
        if (tokens.size() < 2) {
            cout << "option: " << options << " invalid"<<endl;
            cout << "provide two numbers, e.g. 110,140" << endl;
            exit(1);
        }
        min_obs = (double) atof(tokens.at(0).c_str());
        max_obs = (double) atof(tokens.at(1).c_str());
    }
    opt_id ++;

    int strategy = 1;
    if (argc > opt_id) {
        strategy = atoi(argv[opt_id]);
    }
    opt_id ++;
    
    int fill_color = kGreen; // (416)
    if (argc > opt_id){
        fill_color = atoi(argv[opt_id]);
    }
    opt_id ++;
    
    string fix_variables = "";
    if (argc > opt_id)
    {
        fix_variables = argv[opt_id];
    }
    opt_id ++;

    vector<string> np_names;
    if (argc > opt_id) {
        Helper::tokenizeString(string(argv[opt_id]), ',', np_names);
    }
    opt_id ++;
    if (np_names.size() > 0){
        cout << "total NPs: " << np_names.size() << endl;
    }

    gSystem->Load("/afs/cern.ch/user/x/xju/public/src/HggTwoSidedCBPdf_cc.so");
    gSystem->Load("/afs/cern.ch/user/x/xju/public/src/HggScalarLineShapePdf_cc.so");
    gSystem->Load("/afs/cern.ch/user/x/xju/public/src/HggGravitonLineShapePdf_cc.so");
    gSystem->Load("/afs/cern.ch/user/x/xju/public/src/FlexibleInterpVarMkII_cc.so");
    
   
    auto* file_in = TFile::Open(input_name.c_str(), "read");
    auto workspace = (RooWorkspace*) file_in->Get(wsName.c_str());
    auto* mc = (RooStats::ModelConfig*) workspace->obj(mcName.c_str());
    if(!mc) {
        cout << "ERROR: no ModelConfig" << endl;
        exit(1);
    }
    RooSimultaneous* simPdf = NULL;
    try {
        simPdf = (RooSimultaneous*) mc->GetPdf();
    } catch (...){
        cout << "ERROR: no combined pdf in ModelConfig" << endl;
        exit(1);
    }

    RooRealVar* obs = (RooRealVar*) mc->GetObservables()->first();
    auto* mu = (RooRealVar*) workspace->var(muName.c_str()); 
    auto* data =(RooDataSet*) workspace->data(dataName.c_str());
    if(!data) {
        log_err("data(%s) does not exist", dataName.c_str());
        if (dataName.find("asimov") != string::npos) {
            bool do_profile = false;
            data = RooStatsHelper::makeAsimovData(workspace, 0.0, 0.0, mu->GetName(), mcName.c_str(), "combData", do_profile);
        }
        exit(2);
    }
    int nbins = (int) (max_obs-min_obs)/3.;
    if(!obs || !mu) {
        file_in->Close();
        return 0;
    } else {
        mu->Print();
        obs->Print();
        // obs->setRange(min_obs, max_obs);
        /*need to pre-define binning!*/
        // obs->setBins(nbins);
    }

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(strategy);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

    // summary of options
    cout<<" Input: " << input_name << endl;
    cout<<" wsName: " << wsName << endl;
    cout<<" muName: " << muName << endl;
    cout<<" dataName: " << dataName << endl;
    cout<<" fit mode: " << fit_mode << endl;
    cout<<" withData: " << with_data << endl;
    cout<<" Visual Error: " << do_visual_error << endl;
    cout<<" NP size: " << np_names.size() << endl;
    cout<<" Range of obs: [" << min_obs << "-" << max_obs << "] " << endl;
    cout<<" strategy: " << strategy << endl;
    cout<<" Fix variables: " << fix_variables << endl;


    /* unconditional fit*/
    if(fit_mode == 0) 
    { // background only
        mu->setVal(0);
        mu->setConstant();
        RooStatsHelper::fixTermsWithPattern(mc, "ATLAS_");
        // but free the following two, to include sprious signal uncertainties
        auto hgg_bias = workspace->var("ATLAS_Hgg_BIAS");
        auto mRes = workspace->var("ATLAS_mRes");
        if(hgg_bias) hgg_bias->setConstant(0);
        if (mRes) mRes ->setConstant(0);
    } else if(fit_mode == 1) 
    { // signal + background fit
        mu->setConstant(0);
        auto mG = workspace->var("mG");
        auto kappa = workspace->var("GkM");
        auto mX = workspace->var("mX");
        auto wX = workspace->var("wX");
        if (mG) mG->setConstant(0);
        if (kappa) kappa->setConstant(0);
        if (mX) mX->setConstant(0);
        if (wX) wX->setConstant(0);
    } else {}

    /* fix parameters given in option*/
    if(fix_variables.find("gamma") != string::npos) {
        auto nll = RooStatsHelper::createNLL(data, mc);
        RooStatsHelper::minimize(nll, workspace, true);
        RooStatsHelper::fixTermsWithPattern(mc, "gamma_stat");
        delete nll;
    }
    RooStatsHelper::fixVariables(workspace, fix_variables, NULL);

    // print out the message
    mc->GetParametersOfInterest()->Print("v");
    mc->GetNuisanceParameters()->Print("v");


    auto nll = RooStatsHelper::createNLL(data, mc);
    RooFitResult* fit_results = NULL;
    if (do_visual_error && fit_mode >= 0) {
        fit_results = RooStatsHelper::minimize(nll, workspace, true);
    }
    // fit_results = NULL;
    auto* canvas = new TCanvas("c1", "c1", 600, 600);
    gStyle->SetMarkerSize(0.5);
    bool is_log = false;
    if(is_log) canvas->SetLogy();

    const RooCategory& category = *dynamic_cast<const RooCategory*>(&simPdf->indexCat());
   
    RooCatType* obj;
    TList* data_lists = NULL;
    if(data) {
        data_lists = data->split(category, true);
    }
    obs->setRange(min_obs, max_obs);
    obs->setBins(nbins);
    auto* out_file = TFile::Open("out_hist.root", "recreate");

    vector <bool> comp; 
    comp.push_back(true); comp.push_back(false); 
    for( int cp=0; cp<(int)comp.size(); cp++){
      float oldval=mu->getVal();
      if(!comp.at(cp)) {
        std::cout<<"Switching to background only" <<std::endl;
        mu->setVal(0);
      }
      TIter cat_iter(category.typeIterator());
      obj=NULL;
      while( (obj= (RooCatType*)cat_iter()) )
        {
          const char* label_name = obj->GetName();
          RooAbsPdf* pdf = simPdf->getPdf(label_name);
          pdf->Print();
          auto* obs_frame = obs->frame(min_obs, max_obs, nbins);
          obs_frame->SetMarkerSize(0.015);
          int color = 2;
          
          double bkg_evts = pdf->expectedEvents(RooArgSet(*obs));
          string component="bck";
          if(!comp.at(cp)) component="sb";
          auto* hist_bkg = (TH1F*) pdf->createHistogram(Form("hist_%s_%s",component.c_str(),label_name), *obs, RooFit::Binning(nbins));
          std::cout<<"Saving to "<<Form("hist_%s_%s",component.c_str(),label_name)<<std::endl;
          hist_bkg->Scale(bkg_evts/hist_bkg->Integral());
          log_info("number of events: %.2f", bkg_evts);
          // RooCmdArg add_arg = (fit_results==NULL)?RooCmdArg::none():RooFit::VisualizeError(*fit_results, 1.0, kFALSE);
          //RooCmdArg add_arg = (fit_results==NULL)?RooCmdArg::none():RooFit::VisualizeError(*fit_results);
          RooCmdArg add_arg = RooFit::VisualizeError(*fit_results);
          add_arg.Print();
          // possibly with band
          pdf->plotOn(obs_frame, RooFit::LineStyle(1), 
                      RooFit::LineColor(1),
                      RooFit::LineWidth(2),
                      RooFit::FillColor(fill_color),
                      RooFit::Normalization(bkg_evts, RooAbsReal::NumEvent),
                      add_arg
                      );
          // no band
          pdf->plotOn(obs_frame, RooFit::LineStyle(1), 
                      RooFit::LineColor(1),
                      RooFit::LineStyle(2),
                      RooFit::LineWidth(1),
                      RooFit::Normalization(bkg_evts, RooAbsReal::NumEvent)
                      );
          
          /* deal with nuisance parameters */
        if(np_names.size() > 0)
          {
            double sigma_level = 2.0;
            for(auto itr = np_names.begin(); itr != np_names.end(); itr++)
              {
                string np_name(*itr);
                auto* np_var = (RooRealVar*) workspace->var(np_name.c_str());
                if(!np_var) continue;
                np_var->setVal(sigma_level);

                double splusb_evts = pdf->expectedEvents(RooArgSet(*obs));
                log_info("number of events: %.2f with %s up %.f sigma", splusb_evts, np_name.c_str(), sigma_level);
                auto hist_sb_2s_up = (TH1F*) pdf->createHistogram(Form("hist_%s_%s_2up_%s", component.c_str(),np_name.c_str(), label_name), *obs, RooFit::Binning(nbins));
                if(hist_sb_2s_up){
                    hist_sb_2s_up->Scale(splusb_evts/hist_sb_2s_up->Integral());
                    out_file->cd();
                    hist_sb_2s_up->Write();
                    delete hist_sb_2s_up;
                }
                pdf->plotOn(obs_frame, RooFit::LineStyle(7), 
                        RooFit::LineColor(color++),
                        RooFit::LineWidth(2),
                        RooFit::Normalization(splusb_evts, RooAbsReal::NumEvent)
                        );

                np_var->setVal(-1*sigma_level);
                splusb_evts = pdf->expectedEvents(RooArgSet(*obs));
                auto* hist_np_2sdown = (TH1F*) pdf->createHistogram(Form("hist_%s_%s_2down_%s", component.c_str(),np_name.c_str(), label_name), *obs, RooFit::Binning(nbins));
                if(hist_np_2sdown){
                    hist_np_2sdown->Scale(splusb_evts/hist_np_2sdown->Integral());
                    out_file->cd();
                    hist_np_2sdown->Write();
                    delete hist_np_2sdown;
                }
                pdf->plotOn(obs_frame, RooFit::LineStyle(7), 
                        RooFit::LineColor(color++),
                        RooFit::LineWidth(2),
                        RooFit::Normalization(splusb_evts, RooAbsReal::NumEvent)
                        );
                np_var->setVal(0);
            }
        }
        out_file->cd();
        hist_bkg->SetName(Form("hist_%s_%s", component.c_str(),label_name));
        hist_bkg->Write();

        if(data && with_data)
        {
            auto* data_ch = (RooDataSet*) data_lists->At(obj->getVal());
            // double num_data = data_ch->sumEntries();
            data_ch->plotOn(obs_frame, 
                    RooFit::LineStyle(1), 
                    RooFit::LineColor(1),
                    RooFit::LineWidth(2),
                    RooFit::DrawOption("ep")
                    // RooFit::Normalization(num_data, RooAbsReal::NumEvent)
                    );
            cout <<"Data: " << data_ch->sumEntries() << endl;
            auto* hist_data = data_ch->createHistogram(Form("hist_data_%s", label_name), *obs, RooFit::Binning(nbins));
            if(hist_data){
              out_file->cd();
              hist_data->Write();
              delete hist_data;
            }
        }
        if (is_log) obs_frame->GetYaxis()->SetRangeUser(1E-2, 1E2);
        else obs_frame->GetYaxis()->SetRangeUser(0, 10);
        obs_frame->Print();
        obs_frame->Print("v");
        obs_frame->Draw();
        TString out_pdf_name(input_name);
        out_pdf_name.ReplaceAll("root", "pdf");
        canvas->SaveAs(Form("%s_%s", label_name, out_pdf_name.Data()));
        out_file->cd();

        // get error band
        RooCurve* error_band = (RooCurve*) obs_frame->getObject(0);
        RooCurve* nominal = (RooCurve*) obs_frame->getObject(1);
        error_band->SetName(Form("error_band_%s_%s",component.c_str(),label_name));
        nominal ->SetName(Form("nominal_%s_%s",component.c_str(),label_name));
        error_band->Write();
        nominal->Write();
        
        obs_frame->SetName(Form("frame_obs_%s_%s",component.c_str(),label_name));
        obs_frame->Write();
        delete obs_frame;
        }
      mu->setVal(oldval);
    }
    out_file->Close();
    delete canvas;
    
    delete nll;
    file_in->Close();
    return 1;
}
