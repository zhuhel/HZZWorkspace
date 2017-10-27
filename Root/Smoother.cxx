#include "Hzzws/Smoother.h"

#include "TCanvas.h"
#include <RooDataSet.h>
#include <RooNDKeysPdf.h>
#include <RooKeysPdf.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <RooCmdArg.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include "TMath.h"
#include "RooPlot.h"
#include "TLatex.h"

#include <algorithm>

#include "Hzzws/Helper.h"


Smoother::Smoother(const string& outname, std::vector<float> rho, std::string mirror){
    outfile = TFile::Open(outname.c_str(), "RECREATE");
    m_rho = rho;
    m_makeValidPlots=false;
    m_mirror = mirror;
}

Smoother::~Smoother() {
    if(outfile) outfile->Close();
}


void Smoother::smooth(const string& input_name, 
        const string& oname, const string& treename, 
        const RooArgSet &treeobs, const string& cut) const 
{

  TChain* tcut = Helper::loader(input_name, treename);
  std::cout<<"got tree with "<<tcut->GetEntries()<<" total entries, "<<tcut->GetEntries(cut.c_str())<<" passing cut"<<std::endl;
  if (tcut->GetEntries()==0) {
    std::cout<<"no entries! returning null for sample "<<oname<<std::endl;
    return;
  }

  //Add variables to cut on
  RooArgSet obsAndCut(treeobs);
  strvec name_list;
  Helper::getListOfNames(cut, name_list);
  name_list.push_back("weight");
  for(auto name : name_list){
    if (!obsAndCut.find(name.c_str())){
      obsAndCut.add(*(new RooRealVar(name.c_str(), name.c_str(), -1E6, 1E6)));
    }
  }

  std::cout<<"variables to be read:"<<std::endl;
  obsAndCut.Print("v");

  RooDataSet *ds = new RooDataSet(Form("%s_RooDataSet", oname.c_str()), "dataset", 
      obsAndCut, RooFit::Import(*tcut),
      RooFit::Cut(cut.c_str()),
      RooFit::WeightVar("weight")
      );

  ds->Print();

  RooArgList obsList(treeobs);
  RooCmdArg arg2 = RooCmdArg::none();
  RooRealVar *x = (RooRealVar*) obsList.at(0); 

  // 1D pdf
  TH1F *h1 = NULL;
  if (obsList.getSize() == 1){
    RooKeysPdf::Mirror mirror = RooKeysPdf::NoMirror;  // NoMirror, MirrorBoth
    ReadMirror(mirror);

    RooKeysPdf *keyspdf = new RooKeysPdf(Form("%s_RooKeysPdf", oname.c_str()), "keyspdf", *x, *ds, mirror, m_rho[0]);

    h1 =(TH1F*) keyspdf->createHistogram(Form("%s", oname.c_str()), *x, RooFit::Binning(x->getBinning()));
    h1->SetMinimum(0);
    h1->Scale(ds->sumEntries()/h1->Integral());
    h1->SetName(oname.c_str());

    RooPlot* frame=NULL;
    if (m_makeValidPlots){ 
      int bins = ds->numEntries()/60;
      if (bins>100) bins=100;
      if (bins<10) bins=10;
      while (bins>1){
        TH1F* htemp = (TH1F*)ds->createHistogram("htemp",*x,RooFit::Binning(bins));
        bool bad=Helper::TH1FHasEmptyBin(htemp);
        delete htemp;
        if (bad) --bins;
        else break;
      }
      frame= x->frame(RooFit::Bins(bins));
      ds->plotOn(frame);
      keyspdf->plotOn(frame);
      frame->SetMinimum(0.0);
      frame->SetTitle(Form(";%s;A.U.",x->GetName()));
      TLatex* txt1 = new TLatex(0.15,0.85,Form("%s",m_nick.c_str()));
      txt1->SetNDC();
      txt1->SetTextSize(0.04);
      frame->addObject(txt1);
      TLatex* txt2 = new TLatex(0.15,0.80,Form("%s",oname.c_str()));
      txt2->SetNDC();
      txt2->SetTextSize(0.035);
      frame->addObject(txt2);
      //float chi2 = frame->chiSquare();
      //TLatex* txt = new TLatex(0.5,0.25,Form("P(#Chi^{2}/df > %.3f) = %.3f",chi2,TMath::Prob(chi2*bins, bins)));
      //txt->SetNDC();
      //txt->SetTextSize(0.03);
      //frame->addObject(txt);
    }

    if(outfile){
      outfile->cd();
      if (frame) {
        frame->Write(Form("validation_%s",oname.c_str()));
        TCanvas* c = new TCanvas("c","c",800,800);
        frame->Draw();
        TString plotName = outfile->GetName();
        plotName.Replace(plotName.Last('/')+1,999,TString(Form("plots/%s_%s.eps",m_nick.c_str(),oname.c_str())));
        c->Print(plotName.Data());
      }
      h1->Write();
      delete keyspdf;   
    }
  }

  // 2D pdf
  TH2F *h2 = NULL;
  if (obsList.getSize() == 2){
    RooRealVar *y = (RooRealVar*) obsList.at(1); 
    arg2 = RooFit::YVar(*y, RooFit::Binning(y->getBinning()));

    Double_t rhoarr[2] = { m_rho[0], m_rho[1] };
    TVectorD rhovec(2,rhoarr);
    RooNDKeysPdf *keyspdf = new RooNDKeysPdf(Form("%s_RooNDKeysPdf", oname.c_str()), "keyspdf", treeobs, *ds, rhovec, (m_mirror=="NoMirror"?"a":"am"));

    h2 =(TH2F*) keyspdf->createHistogram(Form("%s", oname.c_str()), *x, RooFit::Binning(x->getBinning()), arg2);
    h2->Scale(ds->sumEntries()/h2->Integral());
    h2->SetName(oname.c_str());

    RooPlot* frameX=NULL;
    RooPlot* frameY=NULL;
    if (m_makeValidPlots){ 
      int bins = ds->numEntries()/40;
      if (bins>100) bins=100;
      if (bins<10) bins=10;
      while (bins>1){
        TH1F* htempX = (TH1F*)ds->createHistogram("htemp",*x,RooFit::Binning(bins));
        TH1F* htempY = (TH1F*)ds->createHistogram("htemp",*y,RooFit::Binning(bins));
        bool bad=Helper::TH1FHasEmptyBin(htempX) || Helper::TH1FHasEmptyBin(htempY);
        delete htempX;
        delete htempY;
        if (bad) --bins;
        else break;
      }
      frameX= x->frame(RooFit::Bins(bins));
      frameX->SetMinimum(0.);
      frameY= y->frame(RooFit::Bins(bins));
      frameY->SetMinimum(0.);
      ds->plotOn(frameX);
      ds->plotOn(frameY);
      keyspdf->plotOn(frameX);
      keyspdf->plotOn(frameY);

      float chi2X = frameX->chiSquare();
      float chi2Y = frameY->chiSquare();
      TLatex* txtX = new TLatex(0.5,0.25,Form("P(#Chi^{2}/df > %.3f) = %.3f",chi2X,TMath::Prob(chi2X*bins, bins)));
      TLatex* txtY = new TLatex(0.5,0.25,Form("P(#Chi^{2}/df > %.3f) = %.3f",chi2Y,TMath::Prob(chi2Y*bins, bins)));
      txtX->SetNDC();
      txtY->SetNDC();
      txtX->SetTextSize(0.05);
      txtY->SetTextSize(0.05);
      frameX->addObject(txtX);
      frameY->addObject(txtY);
    }


    if(outfile){
      outfile->cd();
      if (frameX || frameY) {
        if (frameX) frameX->Write(Form("validationX_%s",oname.c_str()));
        if (frameY) frameY->Write(Form("validationY_%s",oname.c_str()));
        TCanvas* c = new TCanvas("c","c",800,800);
        TString plotName = outfile->GetName();
        plotName.Replace(plotName.Last('/')+1,999,TString(Form("plots/%s_X_%s.eps",m_nick.c_str(),oname.c_str())));
        if (frameX) {
          frameX->Draw();
          c->Print(plotName.Data());
          c->Clear();
        }
        if (frameY){
          plotName.ReplaceAll("_X_","_Y_");
          frameY->Draw();
          c->Print(plotName.Data());
        }
      }
      h2->Write();
      delete keyspdf;   
    }
  }

  delete ds;
  delete tcut;

}


void Smoother::ReadMirror(RooKeysPdf::Mirror& mirror) const{
  if (m_mirror=="NoMirror") mirror=RooKeysPdf::NoMirror;
  else if (m_mirror=="MirrorLeft") mirror=RooKeysPdf::MirrorLeft;
  else if (m_mirror=="MirrorRight") mirror=RooKeysPdf::MirrorRight;
  else if (m_mirror=="MirrorBoth") mirror=RooKeysPdf::MirrorBoth;
  else if (m_mirror=="MirrorAsymLeft") mirror=RooKeysPdf::MirrorAsymLeft;
  else if (m_mirror=="MirrorAsymLeftRight") mirror=RooKeysPdf::MirrorAsymLeftRight;
  else if (m_mirror=="MirrorAsymRight") mirror=RooKeysPdf::MirrorAsymRight;
  else if (m_mirror=="MirrorLeftAsymRight") mirror=RooKeysPdf::MirrorLeftAsymRight;
  else if (m_mirror=="MirrorAsymBoth") mirror=RooKeysPdf::MirrorAsymBoth;
  else {
    std::cout<<"Warning: unrecognized mirror option received: "<<m_mirror<<", I am using NoMirror by default"<<std::endl;
  }
}
