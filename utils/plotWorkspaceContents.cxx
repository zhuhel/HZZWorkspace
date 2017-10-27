#include "RooWorkspace.h"
#include "TFile.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <stdio.h>
#include "TObject.h"

int main(int argc, const char** argv){

  using namespace RooFit;

  std::string workspaceFile = "combined.root";
  std::string workspaceName = "combined";
  std::string dsName = "obsData";
  std::string outFileName = "workspaceContents.root";

  bool doPdfs(false);
  bool doFunctions(false);
  bool doVars(false);

  for (int i(1);i<argc; ++i){
    TString arg(argv[i]);

    if (arg=="pdfs" || arg=="p") doPdfs=true;
    if (arg=="functions" || arg=="f") doFunctions=true;
    if (arg=="vars" || arg=="v") doVars=true;
  }

  if (!(doPdfs||doFunctions||doVars)){
    std::cout<<"give valid arguments! options are \"pdfs\", \"functions\", \"vars\" (or \"p\",\"f\",\"v\")"<<std::endl;
    return -1;
  }
  

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  auto file = new TFile(workspaceFile.c_str());
  if (!file) { std::cout<<"no file name "<<workspaceFile<<std::endl; return -1; }

  auto w = (RooWorkspace*)file->Get(workspaceName.c_str());
  if (!w) { std::cout<<"no workspace named "<<workspaceName<<std::endl; return -1;}

  auto data = w->data(dsName.c_str());
  if (!data) { std::cout<<"no dataset named "<<dsName<<std::endl; return -1;}

  auto fileout = new TFile(Form("plots/%s",outFileName.c_str()),"RECREATE");


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (doPdfs){

    fileout->mkdir("pdfs")->cd();

    auto allPdfs = w->allPdfs();

    auto it1 = allPdfs.createIterator();
    while (auto pdf = (RooAbsPdf*)it1->Next()){

      if (TString(pdf->GetName())=="simPdf") continue;

      auto allObs = pdf->getObservables(data);

      if (allObs->getSize()!=0)
        fileout->mkdir(Form("pdfs/%s",pdf->GetName()))->cd();

      auto it2 = allObs->createIterator();
      while (auto obs = it2->Next()){

        if (!dynamic_cast<RooRealVar*>(obs)) continue;

        std::cout<<"going to plot pdf:"<<pdf->GetName()<<" on observable "<<obs->GetName()<<std::endl;
        auto frame = ((RooRealVar*)obs)->frame();
        pdf->plotOn(frame);
        fileout->cd(Form("pdfs/%s",pdf->GetName()));
        frame->Write(Form("pdf__%s__%s",pdf->GetName(),obs->GetName()));
        delete frame;

        auto allPars = pdf->getParameters(data)->selectByName("alpha*");
        if (allPars->getSize()!=0)
          fileout->mkdir(Form("pdfs/%s/sysDependence",pdf->GetName()))->cd();
        
        auto it3 = allPars->createIterator();
        while (auto par = it3->Next()){

          if (!dynamic_cast<RooRealVar*>(par)) continue;

          std::cout<<"going to plot pdf:"<<pdf->GetName()<<" dependence on parameter "<<par->GetName()<<" on observable "<<obs->GetName()<<std::endl;
          auto frame = ((RooRealVar*)obs)->frame();
          pdf->plotOn(frame, Normalization(pdf->expectedEvents(RooArgSet(*((RooRealVar*)obs)))));
          float orig = ((RooRealVar*)par)->getVal();
          ((RooRealVar*)par)->setVal(2.0);
          pdf->plotOn(frame,Normalization(pdf->expectedEvents(RooArgSet(*((RooRealVar*)obs)))),LineColor(kRed));
          ((RooRealVar*)par)->setVal(-2.0);
          pdf->plotOn(frame,Normalization(pdf->expectedEvents(RooArgSet(*((RooRealVar*)obs)))),LineColor(kCyan));
          ((RooRealVar*)par)->setVal(orig);
          fileout->cd(Form("pdfs/%s/sysDependence",pdf->GetName()));
          frame->Write(Form("pdf__%s__%s__%s",pdf->GetName(),obs->GetName(),par->GetName()));
          delete frame;
        }

      }
    }
  }
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (doFunctions){

    fileout->mkdir("functions")->cd();

    auto allFuncs = w->allFunctions();

    std::vector<float> vals;
    std::vector<const char*> names;

    auto it1 = allFuncs.createIterator();
    while (auto func = (RooAbsReal*)it1->Next()){

      auto allParams = func->getParameters(data);

      if (allParams->getSize()!=0)
        fileout->mkdir(Form("functions/%s",func->GetName()))->cd();

      //store values of only non-unit functions
      if (func->getVal()!=1.0){
        vals.push_back(func->getVal());
        names.push_back(func->GetName());
      }

      auto it2 = allParams->createIterator();
      while (auto par = it2->Next()){

  
        auto rrv = dynamic_cast<RooRealVar*>(par);
        if (!rrv) continue;
        
        double def = rrv->getVal();
        double low = rrv->getBinning().lowBound();
        double high = rrv->getBinning().highBound();

        if (low<-1e25 || std::isnan(low) || high>1e25 || std::isnan(high)) continue;

        rrv->setVal(low);
        float a = func->getVal();
        rrv->setVal(high);
        float b = func->getVal();
        rrv->setVal(def);
        float c = func->getVal();

        float min = std::min(std::min(a,b),c);
        float max = std::max(std::max(a,b),c);

        if (a==b) continue;

        std::cout<<"going to plot function:"<<func->GetName()<<" on parameter "<<par->GetName()<<std::endl;
        auto frame = rrv->frame();
        func->plotOn(frame);
        frame->SetMinimum(min-0.25*fabs(max-min));
        frame->SetMaximum(max+0.25*fabs(max-min));
        fileout->cd(Form("functions/%s",func->GetName()));
        frame->Write(Form("function__%s__%s",func->GetName(),par->GetName()));
        delete frame;
      }
    }
    //create summary histogram for non-1 functions

    TH1F* histfunc = new TH1F("_functions_summary","",vals.size(),0,1);
    for (int i(0);i<vals.size();++i){
      histfunc->SetBinContent(i+1,vals[i]);
      histfunc->GetXaxis()->SetBinLabel(i+1,names[i]);
    }
    histfunc->Write(0,TObject::kWriteDelete);
  }
  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (doVars){

    fileout->mkdir("vars")->cd();
    
    auto allVars = w->allVars();
    allVars.remove(*data->get());

    std::vector<float> vals;
    std::vector<const char*> names;
    auto it1 = allVars.createIterator();
    while (auto var = (RooRealVar*)it1->Next()) {

      if (var->getVal()!=0.){
        std::cout<<"storing variable "<<var->GetName()<<" with val "<<var->getVal()<<std::endl;
        vals.push_back(var->getVal());
        names.push_back(var->GetName());
      }
    }

    //create summary histogram for non-0 vars

    TH1F* histvars = new TH1F("_vars_summary","",vals.size(),0,1);
    std::cout<<"size = "<<vals.size()<<std::endl;
    for (int i(0);i<vals.size();++i){
      std::cout<<i<<" "<<names[i]<<" "<<vals[i]<<std::endl;
      histvars->SetBinContent(i+1,vals[i]);
      histvars->GetXaxis()->SetBinLabel(i+1,names[i]);
    }
    histvars->Write(0,TObject::kWriteDelete);
  }
      


  fileout->Write();
  fileout->Close();
  return 1;
}


